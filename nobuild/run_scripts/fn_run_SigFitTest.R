gen_run_aux = function(fname, reference_cat, keep_sigs, n_steps=2000,
                       run_fits=TRUE, filter_dn=FALSE, save_path=NULL, py=NULL) {

  out_fname = strsplit(fname, "/")[[1]]
  out_fname = out_fname[length(out_fname)]
  if (!is.null(save_path) & !dir.exists(save_path)) dir.create(save_path, recursive=TRUE)
  if (file.exists(paste0(save_path, out_fname))) return(cli::cli_text("The file already exixts. Not performing inference."))

  simul_ng = readRDS(fname)$dataset
  G = simul_ng %>% get_n_groups()

  cli::cli_process_start("\nFitting")
  cat(paste("File", fname, "\n"))
  counts_ng = get_input(simul_ng, matrix=TRUE)

  # reorder column names - contexts - to match catalogue ones
  counts_ng = lapply(names(counts_ng), function(name_i) {
			     if (is.null(reference_cat[[name_i]])) return(counts_ng[[name_i]])
			     counts_ng[[name_i]][,colnames(reference_cat[[name_i]])]
		       } ) %>%
    setNames(names(counts_ng))

  fixed_beta = reference_cat

  dn_signames = lapply(names(fixed_beta), function(tid) {
			       if (is.null(fixed_beta[[tid]])) return(get_signames(simul_ng)[[tid]])
			       setdiff(get_signames(simul_ng)[[tid]], rownames(fixed_beta[[tid]]))
		       }) %>%
    setNames(names(fixed_beta))
  reference_sigs_max = lapply(fixed_beta, nrow) %>% unlist()
  if (is.null(reference_sigs_max)) reference_sigs_max = 0
  else reference_sigs_max = min(reference_sigs_max)

  min_K = max(min(sapply(dn_signames, length)) - reference_sigs_max - 5, 0)
  max_K = max(max(sapply(dn_signames, length)) + 5, min_K+5)
  cli::cli_text("min_K {min_K}, max_K {max_K}\n")

  # if (is.null(x_ng.0) & run_fitss
  if (run_fits) {
    # TIME = as.POSIXct(Sys.time(), format = "%H:%M:%S")
    start_time = Sys.time()
    
    x_ng.0 = fit(counts=counts_ng, k_list=min_K:max_K,
                 # cluster=max(G*2, 5),
                 n_steps=n_steps,
                 reference_cat=fixed_beta,
                 keep_sigs=unlist(keep_sigs),
                 # hyperparameters=list("penalty_scale"=0),
                 seed_list=c(10,33,455), filter_dn=filter_dn, store_fits=TRUE,
                 py=py)
    
    time_fit = Sys.time()

    cli::cli_text("Optimal number of signatures: {length(get_signames(x_ng.0)[['SBS']])}")
    
    # saveRDS(x_ng.0, paste0(save_path, "FIT.", out_fname))
    # x_ng.0 = readRDS(paste0(save_path, "FIT.", out_fname))

    if (length(get_denovo_signames(x_ng.0)[["SBS"]]) > 0 &
          length(get_signames(x_ng.0)[["SBS"]]) > 1) {
      x_ng_refined.0 = refine_denovo_signatures(x_ng.0)
      time_refinement = Sys.time()
      x_ng_refined.0 = fit_clustering(x_ng_refined.0, cluster=max(G*2, 5), py=py) %>% merge_clusters()
      time_refinement_clustering = Sys.time()
    }

    x_ng.0 = fit_clustering(x_ng.0, cluster=max(G*2, 5), py=py) %>% merge_clusters()
    time_clustering = Sys.time()

    x_ng.0$time = list("fit"=time_fit - start_time,
                       "clustering"=time_clustering - time_refinement_clustering)
    
    x_ng_refined.0$time = list("fit"=time_fit - start_time,
                               "refinement"=time_refinement - time_fit,
                               "clustering"=time_refinement_clustering - time_refinement)
  }


  cli::cli_process_done()

  final_list = list("dataset"=simul_ng, "fit.0"=x_ng.0, "fit_refined.0"=x_ng_refined.0) #, "fit.N"=x_ng.N)

  if (!is.null(save_path)) saveRDS(final_list, paste0(save_path, out_fname))
  return(final_list)
}

