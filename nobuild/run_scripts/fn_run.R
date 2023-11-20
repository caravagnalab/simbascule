gen_run_aux = function(N, G, seed, private, shared, n_steps=2000,
                       path=NULL, run_fits=FALSE, run_name="") {
  fname = paste0("simul_fit.N", N, ".G", G, ".s", seed, ".", run_name, ".Rds")

  simul_ng = x_ng = NULL
  if (!is.null(path)) {
    fname_fpath = paste0(path, fname)

    if (file.exists(fname_fpath)) simul_ng = readRDS(fname_fpath)$dataset
    if (file.exists(fname_fpath) && !run_fits) x_ng = readRDS(fname_fpath)$fit
  }

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
      create_basilica_obj_simul()
  }

  cat("Simulated dataset generated!")

  if (is.null(x_ng)) {
    cli::cli_process_start(paste0("Fitting ", fname))
    counts_ng = get_input(simul_ng, matrix=TRUE)
    max_K = max(sapply(get_signames(simul_ng), length)) - 1
    min_K = max(0, max_K - 2)
    x_ng.0 = fit(counts=counts_ng, k_list=min_K:max_K, cluster=G*2, n_steps=n_steps,
               reference_cat=list("SBS"=COSMIC_filt[shared$SBS,],
                                  "DBS"=COSMIC_dbs[shared$DBS,]),
               keep_sigs=unlist(shared),
               hyperparameters=list("penalty_scale"=0),
               seed_list=c(10,33,455), filter_dn=TRUE, store_fits=TRUE,
               py=py)
    x_ng.N = fit(counts=counts_ng, k_list=min_K:max_K, cluster=G*2, n_steps=n_steps,
                 reference_cat=list("SBS"=COSMIC_filt[shared$SBS,],
                                    "DBS"=COSMIC_dbs[shared$DBS,]),
                 keep_sigs=unlist(shared),
                 hyperparameters=list("penalty_scale"=N),
                 seed_list=c(10,33,455), filter_dn=TRUE, store_fits=TRUE,
                 py=py)
    cli::cli_process_done()
  }

  final_list = list("dataset"=simul_ng, "fit.0"=x_ng.0, "fit.N"=x_ng.N)

  if (!is.null(path)) saveRDS(final_list, paste0(path, fname))
  return(final_list)
}
