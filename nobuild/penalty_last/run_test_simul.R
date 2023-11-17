devtools::load_all("~/GitHub/simbasilica/")
load_deps()

## Global info ####
path = "~/Dropbox/shared/2022. Basilica/simulations/matched_signals/fits/"
run_name = "no_omega.last_penalty"


## Generate and run ####
private = list(
  "SBS"=c("SBS17b", "SBS4", "SBS7c", "SBS13", "SBS20", "SBS22")
  # "DBS"=c("DBS1", "DBS2", "DBS7", "DBS11", "DBS4", "DBS10")
               )
shared = list(
  "SBS"=c("SBS1","SBS5")
  # "DBS"=c("DBS3","DBS5")
  )

N_list = c(150, 300)
G_list = c(2)
penalty_scales = c(0, 200, 400)
n_gs = expand.grid(N_list, G_list, penalty_scales) %>% dplyr::arrange(Var1)

easypar_pars = lapply(1:nrow(n_gs), function(i) {
  list(N=n_gs[i,"Var1"],
       G=n_gs[i,"Var2"],
       penalty_scale=n_gs[i,"Var3"],
       private=private,
       shared=shared)
})


easypar_fn = function(N, G, penalty_scale, private, shared, path,
                      run_fits=FALSE, run_name="") {
  fname = paste0("simul_fit_", N, ".", G, ".", penalty_scale, ".", run_name, ".Rds")
  fname_fpath = paste0(path, fname)

  simul_ng = x_ng = NULL
  if (file.exists(fname_fpath)) simul_ng = readRDS(fname_fpath)$dataset
  if (file.exists(fname_fpath) && !run_fits) x_ng = readRDS(fname_fpath)$fit

  devtools::load_all("~/GitHub/basilica/")
  devtools::load_all("~/GitHub/simbasilica/")
  reticulate::use_condaenv("basilica-env")
  library(ggplot2)
  py = reticulate::import_from_path("pybasilica", "~/GitHub/pybasilica/")

  if (is.null(simul_ng)) {
    seed_list = list("SBS"=N+G,
                     "DBS"=N+G*2)

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
               seed_list=c(10,33), filter_dn=TRUE, store_fits=TRUE,
               py=py)
    cli::cli_process_done()
  }

  saveRDS(list("dataset"=simul_ng, "fit"=x_ng), paste0(path, fname))
  return(list("dataset"=simul_ng, "fit"=x_ng))
}


## RUN #####
lapply(easypar_pars, function(i) {
  easypar_fn(N=i$N, G=i$G, penalty_scale=i$penalty_scale, private=i$private,
             shared=i$shared, path=path,
             run_fits=TRUE, run_name=run_name)
})

