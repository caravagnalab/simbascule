
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

generate_and_run = function(shared,
                            private,
                            catalogue,
                            comb_matrix,
                            py,
                            private_fracs = list("rare"=0.05, "common"=0.3),
                            fits_path = NULL,
                            data_path = NULL,
                            seeds = 1:30,
                            mut_range = 10:8000,
                            reference_catalogue = COSMIC_filt_merged,
                            input_catalogue = NULL,
                            keep_sigs = c("SBS1", "SBS5"),
                            hyperparameters = NULL,

                            reg_weight = 0.,
                            regularizer = "cosine",

                            initializ_seed = FALSE,
                            initializ_pars_fit = FALSE,
                            save_runs_seed = FALSE,
                            seed_list = c(4,17,22),

                            CUDA = FALSE,
                            do.fits = FALSE,
                            verbose = FALSE,
                            new_model = TRUE,
                            cohort = "",
                            ...) {

  if (!is.null(fits_path) && !dir.exists(fits_path))
    dir.create(fits_path, recursive=T)

  failed = file(paste0(fits_path, "failed_runs.txt"), open="w")

  shared_cat = catalogue[shared,]

  for (i in 1:nrow(comb_matrix)) {
    print(comb_matrix)

    if (is.list(private)) {
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

      if (!is.null(input_catalogue))
        input_sigs = nrow(input_catalogue) else
          input_sigs = length(keep_sigs)

      # min_k = max(0, length(shared) + nrow(denovo_cat) - input_sigs - 5)
      max_k = length(shared) + nrow(denovo_cat) + 5
      k_list = 0:max_k

      idd = paste0("N", comb_matrix$N_vals[i][[1]], ".G", comb_matrix$n_groups_vals[i][[1]], ".s", j)

      cat(paste0(idd, "\n"))

      if (do.fits) {
        fname = paste0("N", comb_matrix$N_vals[i][[1]], ".G",
                       comb_matrix$n_groups_vals[i][[1]], ".s", seeds[j],
                       ".", cohort, ".Rds") %>% stringr::str_replace_all("\\.\\.", ".")

        fits = run_model(x = x$x[[1]],
                         k = k_list,
                         py = py,
                         reference_catalogue = reference_catalogue,
                         input_catalogue = input_catalogue,
                         keep_sigs = keep_sigs,
                         hyperparameters = hyperparameters,

                         reg_weight = reg_weight,
                         CUDA = CUDA,
                         regularizer = regularizer,

                         initializ_seed = initializ_seed,
                         initializ_pars_fit = initializ_pars_fit,
                         save_runs_seed = save_runs_seed,
                         seed_list = seed_list,

                         filtered_cat = TRUE,
                         verbose = verbose,
                         groups = x$groups[[1]] - 1,
                         new_model = new_model,
                         error_file = failed,
                         idd = idd,
                         cohort = cohort,
                         path = fits_path,
                         out_name = fname)

        filename1 = paste0("fit.", idd, ".", cohort, ".Rds") %>%
                           stringr::str_replace_all("\\.\\.", ".")

        filename2 = paste0("fit.hier.", idd, ".", cohort, ".Rds") %>% 
		          stringr::str_replace_all("\\.\\.", ".")

        save_fit(fits$fit1, fits_path, filename1)
        save_fit(fits$fit.hier, fits_path, filename2)
      }
    }
  }

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
                          cosine_limit, seed,
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
  print(x)

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
                     # hyperparameters = NULL,
                     filtered_cat=TRUE,
                     groups=NULL,
                     new_model=FALSE,
                     error_file=NULL,
                     idd="",
                     path=NULL,
                     out_name=NULL) {

  msg1 = paste0("fit.", idd, "\n")
  msg2 = paste0("fit.hier.", idd, "\n")

  expr_fit = (is.null(path) || !paste0("fit.", out_name) %in% list.files(path))
  expr_fit_hier = (is.null(path) || !paste0("fit.hier.", out_name) %in% list.files(path))
  
  x.fit = x.fit.hier = NULL
  if (!new_model) {
    if (expr_fit)
      x.fit = try_run(error_file,
                      expr =
                        fit(..., groups=NULL,
                            input_catalogue=input_catalogue),
                      msg = msg1)

    if (expr_fit_hier)
      x.fit.hier = try_run(error_file,
                           expr =
                             fit(..., groups=groups,
                                 input_catalogue=input_catalogue),
                           msg = msg2)

  } else {
    if (expr_fit)
      x.fit = try_run(error_file,
                      expr =
                        two_steps_inference(..., keep_sigs=keep_sigs,
                                            groups=NULL)$tot,
                      msg = msg1)

    if (expr_fit_hier)
      x.fit.hier = try_run(error_file,
                           expr =
                             two_steps_inference(..., keep_sigs=keep_sigs,
                                                 groups=groups)$tot,
                           msg = msg2)
  }

  return(list("fit1"=x.fit, "fit.hier"=x.fit.hier))
}


try_run = function(error_file, expr, msg) {

  tryCatch(expr = expr,
           error = function(e) {
             writeLines(msg, error_file)
             writeLines(paste(e))
             writeLines(paste(reticulate::py_last_error()))
             return(c(paste(e), paste(reticulate::py_last_error()) ) )
           })
}

