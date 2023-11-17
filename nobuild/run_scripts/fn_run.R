gen_run_aux = function(N, G, seed, private, shared, path,
                       run_fits=FALSE, run_name="") {
  fname = paste0("simul_fit.N", N, ".G", G, ".s", seed, ".", run_name, ".Rds")
  fname_fpath = paste0(path, fname)

  simul_ng = x_ng = NULL
  if (file.exists(fname_fpath)) simul_ng = readRDS(fname_fpath)$dataset
  if (file.exists(fname_fpath) && !run_fits) x_ng = readRDS(fname_fpath)$fit

  if (is.null(simul_ng)) {
    seed_list = list("SBS"=seed,
                     "DBS"=seed*2)

    simul_ng = generate_simulation_dataset_matched(N=N, G=G,
                                                   private=private,
                                                   shared=shared,
                                                   alpha_range=c(0.25,0.3),
                                                   alpha_sigma=0.2,
                                                   py=py,
                                                   seed=seed_list) %>%
      create_basilica_obj_simul()
  }

  cat("Simulated dataset generated!")

  if (is.null(x_ng)) {
    cli::cli_process_start(paste0("Fitting ", fname))
    counts_ng = get_input(simul_ng, matrix=TRUE)
    max_K = max(sapply(get_signames(simul_ng), length)) - 1
    min_K = max(0, max_K - 2)
    x_ng = fit(counts=counts_ng, k_list=min_K:max_K, cluster=NULL, n_steps=2000,
               reference_cat=list("SBS"=COSMIC_filt[shared$SBS,],
                                  "DBS"=COSMIC_dbs[shared$DBS,]),
               keep_sigs=unlist(shared),
               hyperparameters=list("penalty_scale"=penalty_scale),
               seed_list=c(10,33,455), filter_dn=TRUE, store_fits=TRUE,
               py=py)
    cli::cli_process_done()
  }

  saveRDS(list("dataset"=simul_ng, "fit"=x_ng), paste0(path, fname))
  return(list("dataset"=simul_ng, "fit"=x_ng))
}
