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
                            # shared,
                            # private,
                            catalogue,
                            py,
                            shared = NULL,
                            private = NULL,
                            private_fracs = list("rare"=0.05, "common"=0.3),
                            fits_path = NULL,
                            data_path = NULL,
                            seeds = 1:30,
                            mut_range = 10:8000,
                            reference_catalogue = COSMIC_filt,
                            input_catalogue = NULL,
                            keep_sigs = c("SBS1", "SBS5"),
                            hyperparameters = NULL,
                            lr = 0.005,
                            n_steps = 1500,
                            nonparametric = FALSE,
                            enforce_sparsity = TRUE,

                            reg_weight = 1.,
                            regularizer = "cosine",
                            new_hier = FALSE,
                            regul_denovo = TRUE,

                            initializ_seed = FALSE,
                            initializ_pars_fit = TRUE,
                            save_runs_seed = TRUE,
                            seed_list = c(4,17,22),

                            CUDA = FALSE,
                            do.fits = TRUE,
                            verbose = FALSE,
                            new_model = TRUE,
                            cohort = "",

                            check_present = TRUE,
                            inference_type = c("flat","hier","clust"),
                            ...) {

  if ("shared" %in% colnames(comb_matrix))
    shared = comb_matrix$shared %>% unlist() %>% unique()

  cat(paste("regularizer =", regularizer, "- reg_weight =", reg_weight, "- inference_type =", paste(inference_type, collapse=", "), "\n"))

  if (do.fits && !is.null(fits_path) && !dir.exists(fits_path)) {
    dir.create(fits_path, recursive=T)
    failed = file(paste0(fits_path, "failed_runs.txt"), open="w")
  }

  shared_cat = catalogue[shared,]

  for (i in 1:nrow(comb_matrix)) {

    if ("rare" %in% colnames(comb_matrix)) {
      private_rare = comb_matrix[i,] %>% dplyr::pull(rare) %>% unlist()
      private_common = comb_matrix[i,] %>% dplyr::pull(common) %>% unlist()
    } else if (is.list(private)) {
      private_common = private$common
      private_rare = private$rare
    } else {
      private_common = sample(private, comb_matrix$n_priv_comm[i])
      tmp = setdiff(private, private_common)
      private_rare = sample(tmp, comb_matrix$n_priv_rare[i])
    }

    denovo_cat = catalogue[c(private_common, private_rare),]

    for (j in seeds) {

      x = single_dataset(
        N = comb_matrix$N_vals[i][[1]],
        n_groups = comb_matrix$n_groups_vals[i][[1]],
        samples_per_group = comb_matrix$samples_per_group[i][[1]],
        reference_cat = shared_cat,
        denovo_cat = denovo_cat,
        private_sigs = list("rare" = private_rare, "common" = private_common),
        private_fracs = private_fracs,
        mut_range = mut_range,
        seed = j,
        out_path = data_path,
        cohort_name = cohort
      )

      if (!do.fits) next

      if (!is.null(input_catalogue))
        input_sigs = nrow(input_catalogue) else
          input_sigs = length(keep_sigs)

      # min_k = max(0, length(shared) + nrow(denovo_cat) - input_sigs - 5)
      max_k = length(shared) + nrow(denovo_cat) - length(keep_sigs) + 2
      min_k = max(0, max_k - 4)
      k_list = min_k:max_k

      min_cl = max(comb_matrix$n_groups_vals[i][[1]] - 3, 1)
      max_cl = comb_matrix$n_groups_vals[i][[1]] + 3
      cluster_list = min_cl:max_cl

      idd = paste0("N", comb_matrix$N_vals[i][[1]], ".G", comb_matrix$n_groups_vals[i][[1]], ".s", j)

      cat(paste0(idd, "\n"))

      cat(paste("cluster_list =", paste(cluster_list, collapse=","), "\n"))
      cat(paste("k_list =", paste(k_list, collapse=","), "\n"))

      if (do.fits) {
        fname = paste0("N", comb_matrix$N_vals[i][[1]], ".G",
                       comb_matrix$n_groups_vals[i][[1]], ".s", seeds[j],
                       ".", cohort, ".Rds") %>% stringr::str_replace_all("\\.\\.", ".")

        run_model(x = x$x[[1]],
                  k = k_list,
                  py = py,
                  reference_catalogue = reference_catalogue,
                  input_catalogue = input_catalogue,
                  keep_sigs = keep_sigs,
                  hyperparameters = hyperparameters,
                  n_steps = n_steps,
                  lr = lr,
                  enforce_sparsity = enforce_sparsity,
                  nonparametric = nonparametric,

                  reg_weight = reg_weight,
                  CUDA = CUDA,
                  regularizer = regularizer,
                  new_hier = new_hier,
                  regul_denovo = regul_denovo,

                  initializ_seed = initializ_seed,
                  initializ_pars_fit = initializ_pars_fit,
                  save_runs_seed = save_runs_seed,
                  seed_list = seed_list,

                  filtered_cat = TRUE,
                  verbose = verbose,
                  groups = x$groups[[1]] - 1,
                  cluster_list = cluster_list,

                  new_model = new_model,
                  error_file = failed,
                  idd = idd,
                  cohort = cohort,
                  path = fits_path,
                  out_name = fname,
                  check_present = check_present,
                  inference_type = inference_type)
      }
    }
  }

  write.csv(comb_matrix %>%
              dplyr::mutate(dplyr::across(dplyr::where(is.list),
                                          function(x) paste0(x, collapse=","))),
            file=paste0(data_path, "data_settings.csv"), row.names=F)

  if (do.fits && !is.null(fits_path) && !dir.exists(fits_path))
    close(failed)
}


#' Function to generate a synthetic dataset
#'
#' @param N Number of samples
#' @param n_groups add
#' @param samples_per_group add
#' @param reference_cat add
#' @param denovo_cat add
#' @param private_sigs add
#' @param private_fracs add
#' @param cosine_limit add
#' @param seed add
#' @param reference_cosine add
#' @param denovo_cosine add
#' @param mut_range add
#' @param cohort_name add
#' @param out_path add
#'
#' @return The generated data
#' @export single_dataset

single_dataset = function(N, n_groups, samples_per_group,
                          reference_cat, denovo_cat,
                          private_sigs, private_fracs,
                          seed, cosine_limit=NULL,
                          reference_cosine=NULL, denovo_cosine=NULL,
                          mut_range=10:8000, cohort_name="",
                          out_path=NULL) {

  groups = sample(1:n_groups, N, replace=T)
  while (!all(lapply(1:n_groups, function(n) length(groups[groups==n]) %in% samples_per_group) %>% unlist()))
    groups = sample(1:n_groups, N, replace=T)

  idd = paste0("N", N, ".G", n_groups, ".s", seed)

  out_name = paste0("simul.N", N, ".G", n_groups, ".s", seed, ".", cohort_name, ".Rds") %>%
    stringr::str_replace_all("\\.\\.", ".")

  if (!is.null(out_path) && paste0(out_name) %in% list.files(paste0(out_path)))
    return(readRDS(paste0(out_path, out_name)))

  x = generate.data(
    reference_catalogue=reference_cat,
    denovo_catalogue=denovo_cat,
    reference_cosine=reference_cosine,
    denovo_cosine=denovo_cosine,
    targetX=-1,
    inputX=NULL,
    similarity_limit=cosine_limit,
    groups=groups,
    private_sigs=private_sigs,
    private_fracs=private_fracs,
    mut_range=mut_range,
    seed=seed)

  x$private_rare = list(private_sigs$rare)
  x$private_common = list(private_sigs$common)

  if (is.null(out_path)) return(x)

  if (!dir.exists(out_path))
    dir.create(out_path, recursive=T)

  if (cohort_name == "")
    saveRDS(x, paste0(out_path, "simul.", idd, ".Rds"))
  else
    saveRDS(x, paste0(out_path, "simul.", idd, ".", cohort_name, ".Rds"))

  return(x)
}


save_fit = function(x.fit, path, filename, check_present=FALSE) {
  if (is.null(path)) return()

  if (is.null(x.fit)) return()

  if (!dir.exists(path))
    dir.create(path, recursive=T)

  if (check_present && filename %in% list.files(paste0(path)))
    return()

  saveRDS(x.fit, paste0(path, filename))
}


run_model = function(...,
                     input_catalogue=NULL,
                     keep_sigs = c("SBS1","SBS5"),
                     filtered_cat=TRUE,
                     enforce_sparsity = TRUE,
                     nonparametric = FALSE,
                     groups=NULL,
                     cluster_list=NULL,
                     new_model=TRUE,
                     error_file=NULL,
                     idd="",
                     cohort="",
                     path=NULL,
                     out_name=NULL,
                     new_hier = FALSE,
                     regul_denovo =TRUE,
                     inference_type = c("flat","hier","clust"),
                     check_present = TRUE) {

  msg1 = paste0("fit.", idd, "\n")
  msg2 = paste0("fit_hier.", idd, "\n")
  msg3 = paste0("fit_clust.", idd, "\n")

  expr_fit = (!is.null(path) &&
              !(check_present && paste0("fit.", out_name) %in% list.files(path)) &&
              "flat" %in% inference_type)
  expr_fit_hier = (!is.null(path) &&
                   !(check_present && paste0("fit_hier.", out_name) %in% list.files(path)) &&
                   "hier" %in% inference_type)
  expr_fit_clust = (!is.null(path) &&
                    !(check_present && paste0("fit_clust.", out_name) %in% list.files(path)) &&
                    "clust" %in% inference_type)

  x.fit = x.fit.hier = x.fit.clust = NULL
  # if (!new_model) {
  #   if (expr_fit)
  #     x.fit = try_run(error_file,
  #                     expr =
  #                       fit(..., groups=NULL,
  #                           input_catalogue=input_catalogue),
  #                     msg = msg1)

  #   if (expr_fit_hier)
  #     x.fit.hier = try_run(error_file,
  #                          expr =
  #                            fit(..., groups=groups,
  #                                input_catalogue=input_catalogue),
  #                          msg = msg2)

  # } else {

  filename1 = paste0("fit.", idd, ".", cohort, ".Rds") %>%
                      stringr::str_replace_all("\\.\\.", ".")

  filename2 = paste0("fit_hier.", idd, ".", cohort, ".Rds") %>%
        stringr::str_replace_all("\\.\\.", ".")

  filename3 = paste0("fit_clust.", idd, ".", cohort, ".Rds") %>%
        stringr::str_replace_all("\\.\\.", ".")

  if (expr_fit) {
    cli::cli_process_start("Running non-hierarchical fit")
    x.fit = try_run(error_file,
                    expr =
                      two_steps_inference(..., cohort=cohort, keep_sigs=keep_sigs,
                                          groups=NULL, new_hier=FALSE, enforce_sparsity2 = enforce_sparsity,
                                          nonparametric=FALSE, regul_denovo=regul_denovo),
                    msg = msg1)
    cli::cli_process_done()

    cat(paste("Flat run class:", class(x.fit), "\n"))
    save_fit(x.fit, path, filename1, check_present=check_present)
  }

  if (expr_fit_hier) {
    cli::cli_process_start("Running hierarchical fit")
    x.fit.hier = try_run(error_file,
                          expr =
                            two_steps_inference(..., cohort=cohort, keep_sigs=keep_sigs,
                                                groups=groups, new_hier=new_hier, enforce_sparsity2=enforce_sparsity,
                                                nonparametric=FALSE, regul_denovo=regul_denovo),
                          msg = msg2)
    cli::cli_process_done()

    cat(paste("Hierarchical run class:", class(x.fit.hier), "\n"))
    save_fit(x.fit.hier, path, filename2, check_present=check_present)
  }

  if (expr_fit_clust) {
    cli::cli_process_start("Running clustering fit")
    x.fit.clust = try_run(error_file,
                            expr =
                              two_steps_inference(..., cohort=cohort, keep_sigs=keep_sigs,
                                                  groups=NULL, new_hier=new_hier, enforce_sparsity2 = enforce_sparsity,
                                                  clusters=cluster_list, nonparametric=nonparametric, 
                                                  regul_denovo=regul_denovo),
                            msg = msg3)
    cli::cli_process_done()

    cat(paste("Clustering run class:", class(x.fit.clust), "\n"))
    save_fit(x.fit.clust, path, filename3, check_present=check_present)
  }

  # return(list("fit1"=x.fit, "fit.hier"=x.fit.hier, "fit.clust"=x.fit.clust))
}


try_run = function(error_file, expr, msg) {
  tryCatch(expr = expr,
           error = function(e) {
             writeLines(msg, error_file)
             writeLines(paste(e))
             writeLines(paste(reticulate::py_last_error()))
             return(NULL)
           })
}

