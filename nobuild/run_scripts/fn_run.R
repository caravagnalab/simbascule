gen_run_aux = function(N, G, seed, private, shared, reference_cat, n_steps=2000,
                       path=NULL, data_path=NULL, run_fits=FALSE, run_name="", filter_dn=FALSE) {
  fname = paste0("simul_fit.N", N, ".G", G, ".s", seed, ".", run_name, ".Rds")

  simul_ng = x_ng.0 = x_ng.N = NULL
  if (!is.null(data_path)) {
    # if (!dir.exists(path)) dir.create(path, recursive=TRUE)

    fname_fpath = paste0(data_path, fname)
    if (file.exists(fname_fpath)) simul_ng = readRDS(fname_fpath)$dataset
    if (file.exists(fname_fpath)) {
      x_ng.0 = readRDS(fname_fpath)$fit.0
      x_ng.N = readRDS(fname_fpath)$fit.N
    }
  }

  if (!is.null(simul_ng) & !is.null(x_ng.0) & !is.null(x_ng.N))
    return(list("dataset"=simul_ng, "fit.0"=x_ng.0, "fit.N"=x_ng.N))

  if (is.null(simul_ng)) {
    seed_list = list("SBS"=seed,
                     "DBS"=seed*2)

    simul_ng = generate_simulation_dataset_matched(N=N, G=G,
                                                   private=private,
                                                   shared=shared,
                                                   alpha_range=c(0.16,0.25),
                                                   alpha_sigma=0.15,
                                                   py=py,
                                                   seed=seed_list) %>%
      create_bascule_obj_simul()
  }

  cat("Simulated dataset generated!\n\n")

  cli::cli_process_start(paste0("Fitting ", fname))
  counts_ng = get_input(simul_ng, matrix=TRUE)
  fixed_beta = list("SBS"=reference_cat[["SBS"]][shared$SBS,],
                    "DBS"=reference_cat[["DBS"]][shared$DBS,])
  dn_signames = lapply(names(fixed_beta), function(tid) setdiff(get_signames(simul_ng)[[tid]], rownames(fixed_beta[[tid]]))) %>%
    setNames(names(fixed_beta))
  min_K = max(min(sapply(dn_signames, length)) - 2, 0)
  max_K = max(max(sapply(dn_signames, length)) + 2, min_K+1)
  cat(paste0("min_K ", min_K, ", max_K ", max_K, "\n"))

  if (is.null(x_ng.0) & run_fits)
    x_ng.0 = fit(counts=counts_ng, k_list=min_K:max_K, cluster=G*2, n_steps=n_steps,
               reference_cat=fixed_beta,
               keep_sigs=unlist(shared),
               hyperparameters=list("penalty_scale"=0),
               seed_list=c(10,33,455), filter_dn=filter_dn, store_fits=TRUE,
               py=py)
  # if (is.null(x_ng.N) & run_fits)
  #   x_ng.N = fit(counts=counts_ng, k_list=min_K:max_K, cluster=G*2, n_steps=n_steps,
  #                reference_cat=fixed_beta,
  #                keep_sigs=unlist(shared),
  #                hyperparameters=list("penalty_scale"=N),
  #                seed_list=c(10,33,455), filter_dn=filter_dn, store_fits=TRUE,
  #                py=py)

  cli::cli_process_done()

  final_list = list("dataset"=simul_ng, "fit.0"=x_ng.0, "fit.N"=x_ng.N)

  if (!is.null(path) & !dir.exists(path)) dir.create(path, recursive=TRUE)
  if (!is.null(path)) saveRDS(final_list, paste0(path, fname))
  return(final_list)
}

