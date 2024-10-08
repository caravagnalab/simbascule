#' Generate synthetic data and run the fit
#'
#' @param shared Signatures shared across groups
#' @param private Signatures private to one group
#' @param catalogue Reference catalogue
#' @param comb_matrix Combinations of parameters to test
#' @param py Either \code{NULL} or the \code{Python} package to use
#' @param fits_path Path to store the fits
#' @param data_path Path to store the data
#' @param seeds List of seeds to use
#' @param mut_range Number of mutation range for each patient
#' @param reference_catalogue Reference catalogue to use for the inference
#' @param input_catalogue Input catalogue to use for the fit
#' @param reg_weight Regularization weight
#' @param regularizer Type of regularization
#' @param CUDA Logical as whether run it on GPU or not
#' @param do.fits Logical as whether to run the fits
#' @param verbose Logical as whether show the progression
#' @param new_model Logical as whether to run the new model
#'
#' @return nothing
#' @export generate_and_run

generate_and_run = function(comb_matrix,
                            py,

                            fits_path = NULL,
                            data_path = NULL,
                            seeds = 1:30,

                            ## data generation
                            catalogue_sbs,
                            alpha_range = c(.15,0.2),
                            alpha_sigma = 0.1,
                            pi_conc = 1,
                            frac_rare = 1.,
                            n_muts_range = 100:5000,
                            shared_sbs = c("SBS1", "SBS5"),

                            ## inference
                            reference_catalogue = COSMIC_filt,
                            subset_reference = c("SBS1", "SBS5"),
                            keep_sigs = c("SBS1", "SBS5"),
                            hyperparameters = NULL,
                            lr = 0.005,
                            n_steps = 3000,
                            nonparametric = TRUE,
                            enforce_sparsity = TRUE,

                            reg_weight = 0.,
                            regularizer = "cosine",
                            regul_denovo = TRUE,

                            store_fits = FALSE,
                            seed_list = c(4,17,22),

                            CUDA = TRUE,
                            do.fits = FALSE,
                            cohort = "",

                            check_present = TRUE,
                            check_linear_comb = FALSE,

                            inference_type = c("flat","hier","clust"),
                            ...) {

  cat(paste("regularizer =", regularizer, "- reg_weight =", reg_weight, "- inference_type =", paste(inference_type, collapse=", "), "\n"))

  if (do.fits && !is.null(fits_path) && !dir.exists(fits_path))
    dir.create(fits_path, recursive=T)

  failed = paste0(fits_path, "failed_runs.txt")
  if (!file.exists(failed)) file.create(failed)

  for (i in 1:nrow(comb_matrix)) {

    for (j in seeds) {

      idd = paste0("N", comb_matrix$N_vals[i][[1]], ".G",
                   comb_matrix$n_groups_vals[i][[1]], ".s", j)
      cat(paste0(idd, "\n"))

      x = generate_data_aux(N=comb_matrix$N_vals[i][[1]],
                            G=comb_matrix$n_groups_vals[i][[1]],
                            catalogue_sbs=catalogue_sbs,
                            alpha_range=alpha_range,
                            alpha_sigma=alpha_sigma,
                            seed=j,
                            shared_sbs=shared_sbs,
                            pi_conc=pi_conc,
                            frac_rare=frac_rare,
                            n_muts_range=n_muts_range,
                            cohort=cohort, idd=idd,
                            out_path=data_path)

      if (!do.fits) next

      max_k = nrow(x$beta[[1]]) - length(subset_reference) + 2
      min_k = max(0, max_k - 4)
      k_list = min_k:max_k

      min_cl = max(comb_matrix$n_groups_vals[i][[1]] - 3, 1)
      max_cl = comb_matrix$n_groups_vals[i][[1]] + 3
      cluster_list = min_cl:max_cl

      cat(paste0(idd, "\n"))
      cat(paste("cluster_list =", paste(cluster_list, collapse=","), "\n"))
      cat(paste("k_list =", paste(k_list, collapse=","), "\n"))

      if (do.fits) {
        fname = paste0("N", comb_matrix$N_vals[i][[1]], ".G",
                       comb_matrix$n_groups_vals[i][[1]], ".s", seeds[j],
                       ".", cohort, ".Rds") %>% stringr::str_replace_all("\\.\\.", ".")

        run_model(x = x$counts[[1]],
                  k_list = k_list,
                  py = py,
                  reference_catalogue = reference_catalogue,
                  subset_reference = subset_reference,
                  keep_sigs = keep_sigs,
                  hyperparameters = hyperparameters,
                  n_steps = n_steps,
                  lr = lr,
                  enforce_sparsity = enforce_sparsity,
                  nonparametric = nonparametric,

                  reg_weight = reg_weight,
                  CUDA = CUDA,
                  regularizer = regularizer,
                  regul_denovo = regul_denovo,

                  store_fits = store_fits,
                  seed_list = seed_list,

                  filtered_catalogue = TRUE,
                  # groups = x$groups[[1]] - 1,
                  cluster_list = cluster_list,

                  error_file = failed,
                  idd = idd,
                  cohort = cohort,
                  path = fits_path,
                  out_name = fname,
                  check_present = check_present,
                  check_linear_comb = check_linear_comb,
                  inference_type = inference_type)
      }
    }
  }

  col.names = F
  if (!file.exists(paste0(data_path, "data_settings.csv"))) col.names = T

  write.table(comb_matrix %>%
                dplyr::mutate(dplyr::across(dplyr::where(is.list),
                                            function(x) paste0(x, collapse=","))),
              file=paste0(data_path, "data_settings.csv"), row.names=F, col.names=col.names,
              append=T, sep=",", quote=FALSE)

}

run_model = function(...,
                     k_list,
                     reference_catalogue = COSMIC_filt,
                     subset_reference = c("SBS1","SBS5"),
                     nonparametric = TRUE,
                     cluster_list = NULL,
                     error_file = NULL,
                     idd = "",
                     cohort = "",
                     path = NULL,
                     out_name = NULL,
                     regul_denovo =TRUE,
                     inference_type = c("flat","hier","clust"),
                     check_present = TRUE,
                     check_linear_comb = FALSE) {

  msg1 = paste0("fit.", idd, "\n")
  msg3 = paste0("fit_clust.", idd, "\n")

  expr_fit = (!is.null(path) && "flat" %in% inference_type)
  expr_fit_clust = (!is.null(path) && "clust" %in% inference_type)

  x.fit = x.fit.clust = NULL

  filename1 = paste0("fit.", idd, ".", cohort, ".Rds") %>%
                      stringr::str_replace_all("\\.\\.", ".")
  filename3 = paste0("fit_clust.", idd, ".", cohort, ".Rds") %>%
        stringr::str_replace_all("\\.\\.", ".")

  if (expr_fit) {
    cli::cli_process_start("Running non-hierarchical fit")

    x.fit = run_single_fit(..., pattern="fit.", path=path,
                           k_list = k_list, cluster_list=NULL,
                           reference_catalogue = reference_catalogue,
                           subset_reference = subset_reference,
                           cohort=cohort, error_file=error_file,
                           out_name=out_name,
                           nonparametric=FALSE, regul_denovo=regul_denovo,
                           check_present=check_present,
                           check_linear_comb=check_linear_comb, msg=msg1)

    cli::cli_process_done()

    cat(paste("Flat run class:", class(x.fit), "\n"))
    save_fit(x.fit, path, filename1, check_present=check_present,
             check_linear_comb=check_linear_comb)
  }

  if (expr_fit_clust) {
    cli::cli_process_start("Running clustering fit")

    x.fit.clust = run_single_fit(..., pattern="fit_clust.", path=path,
                                 out_name=out_name, k_list = k_list,
                                 cluster_list=cluster_list,
                                 reference_catalogue=reference_catalogue,
                                 subset_reference=subset_reference,
                                 cohort=cohort,
                                 error_file=error_file,
                                 nonparametric=nonparametric,
                                 regul_denovo=regul_denovo,
                                 check_present=check_present,
                                 check_linear_comb=check_linear_comb, msg=msg3)

    cli::cli_process_done()

    cat(paste("Clustering run class:", class(x.fit.clust), "\n"))
    save_fit(x.fit.clust, path, filename3, check_present=check_present,
             check_linear_comb=check_linear_comb)
  }
}


run_single_fit = function(...,
                          k_list,
                          pattern,
                          reference_catalogue = COSMIC_filt,
                          subset_reference = c("SBS1","SBS5"),
                          nonparametric = TRUE,
                          cluster_list = NULL,
                          error_file = NULL,
                          idd = "",
                          cohort = "",
                          path = NULL,
                          out_name = NULL,
                          regul_denovo =TRUE,
                          check_present = TRUE,
                          check_linear_comb = TRUE,
                          msg = "") {

  if (check_present && paste0(pattern, out_name) %in% list.files(path)) {
    if (!check_linear_comb) return(NULL)

    x.fit = readRDS(paste0(path, pattern, out_name))
  } else {
    x.fit = try_run(error_file,
                    expr =
                      fit(..., k = k_list,
                          reference_catalogue=reference_catalogue[subset_reference,],
                          nonparametric = nonparametric,
                          clusters = cluster_list,
                          cohort = cohort,
                          regul_denovo = regul_denovo),
                    msg = msg)
  }

  cat("AFTER READING FILE RDS\n")

  if (check_linear_comb && !is.null(get_denovo_signatures(x.fit))) {
    lc = filter_signatures_QP(sign1=get_denovo_signatures(x.fit),
                              sign2=reference_catalogue,
                              return_weights=FALSE,
                              filt_pi=0.2)
    new_sigs = unique(c(subset_reference, unlist(lc)))
    new_min_k = max(0, k_list[1] - length(setdiff(unlist(lc), subset_reference)))
    k_list = new_min_k:k_list[length(k_list)]

    if (any(!new_sigs %in% subset_reference)) {
      cat(paste(paste(new_sigs, collapse=","), "\n"))
      cat(paste(paste(k_list, collapse=","), "\n"))
      cat(paste(paste(cluster_list, collapse=","), "\n"))

      x.fit_new = try_run(error_file,
                  expr =
                    fit(..., k = k_list,
                        reference_catalogue=reference_catalogue[new_sigs, ],
                        nonparametric = nonparametric,
                        clusters = cluster_list,
                        cohort = cohort,
                        regul_denovo = regul_denovo),
                  msg = paste(msg, "checking linear comb"))
      x.fit$lc_check = x.fit_new
    }
  }

  return(x.fit)
}


try_run = function(error_file, expr, msg) {
  tryCatch(expr = expr,
           error = function(e) {
            print(paste(e))
            write(msg, file=error_file, append=T)
            write(paste(e), file=error_file, append=T)
            write(paste(reticulate::py_last_error()), file=error_file, append=T)
            return(NULL)
           })
}


save_fit = function(x.fit, path, filename, check_present=FALSE, check_linear_comb=FALSE) {
  if (is.null(path)) return()

  if (is.null(x.fit)) return()

  if (!dir.exists(path))
    dir.create(path, recursive=T)

  if (check_present && !check_linear_comb && filename %in% list.files(paste0(path)))
    return()

  saveRDS(x.fit, paste0(path, filename))
}

