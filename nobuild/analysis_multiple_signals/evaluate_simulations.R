devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")

source("~/GitHub/simbasilica/nobuild/analysis_multiple_signals/evaluate_aux_fns.R")

run_id = "clustering.matched.2011"
path = paste0("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.", run_id, "/")
files = list.files(path, full.names=T, pattern=".Rds")

runids = c("Autoguide", "ManualGuide")
fitnames = c("x.fit0.auto", "x.fit0.man")


merge_clusters = function(x, cutoff=0.9) {
  alpha_prior = get_centroids(x, matrix=T)
  if (nrow(alpha_prior) == 1) return(x)
  cosine_simil = lsa::cosine(t(alpha_prior)) %>% as.data.frame()
  cosine_simil[lower.tri(cosine_simil, diag=T)] = 0

  rownames(cosine_simil) = colnames(cosine_simil) = rownames(alpha_prior)

  merging = cosine_simil %>% tibble::rownames_to_column(var="gid1") %>%
    reshape2::melt(id="gid1", variable.name="gid2", value.name="cosine") %>%
    dplyr::mutate(gid1=as.character(gid1), gid2=as.character(gid2)) %>%
    dplyr::filter(cosine > cutoff) %>% dplyr::select(-cosine) %>%
    dplyr::group_by(gid1) %>%
    dplyr::summarise(cl_old=list(c(gid1,gid2) %>% unique())) %>% dplyr::ungroup() %>%
    dplyr::rename(cl_name=gid1)

  if (nrow(merging) == 0) return(x)

  grps = get_cluster_assignments(x)
  for (i in 1:nrow(merging)) {
    old_cl = dplyr::pull(merging[i,], cl_old)[[1]]
    new_cl = dplyr::pull(merging[i,], cl_name)
    grps = grps %>% dplyr::mutate(clusters=ifelse(clusters %in% old_cl,
                                                  new_cl, clusters))
  }

  x$clustering$clusters = grps
  return(x)
}


# all_stats = lapply(files, function(fname) {
#   stats_single_data(fname, names_fits=fitnames %>% setNames(runids))
# }) %>% dplyr::bind_rows()
# saveRDS(all_stats, paste0("~/Dropbox/shared/2022. Basilica/simulations/stats_", run_id, ".Rds"))


# stats_plots = lapply(runids, function(penalty_id) {
#   stats_p = all_stats %>% dplyr::filter(penalty==penalty_id)
#   make_plots_stats(stats_p)
# }) %>% setNames(runids)
# pdf(paste0("~/Dropbox/shared/2022. Basilica/simulations/stats_", run_id, ".pdf"), height=10, width=25)
# patchwork::wrap_plots(stats_plots) & patchwork::plot_annotation(title=paste(runids, collapse=" and "))
# dev.off()


lapply(files, function(fname) {
  tmp = strsplit(fname, split="/")[[1]]; fit_id = tmp[length(tmp)]
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  types = get_types(x.simul)

  plots = lapply(fitnames, function(fitname) {
    x.fit = simul_fit[[fitname]] %>%
      # rename_dn_expos() %>%
      merge_clusters()
    assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul)
    added = sapply(assigned_missing, function(i) i[["added_fp"]]) %>% unlist() %>% paste(collapse=",")
    missing = sapply(assigned_missing, function(i) i[["missing_fn"]]) %>% unlist() %>% paste(collapse=",")

    caption = paste("Added signatures (FP):", added, "- Missing signatures (FN):", missing)

    patchwork::wrap_plots(plot_fit(x.fit),
                          plot_fit(x.simul), ncol=1) &
      patchwork::plot_annotation(title=paste0(fit_id, " , fitname: ", fitname),
                                 subtitle="Fit (top) and simulated (bottom)",
                                 caption=caption)
  })

  QC = lapply(fitnames, function(fitname) {
    x.fit = simul_fit[[fitname]] %>%
      # rename_dn_expos() %>%
      merge_clusters()
    assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul)
    added = sapply(assigned_missing, function(i) i[["added_fp"]]) %>% unlist() %>% paste(collapse=",")
    missing = sapply(assigned_missing, function(i) i[["missing_fn"]]) %>% unlist() %>% paste(collapse=",")

    caption = paste("Added signatures (FP):", added, "- Missing signatures (FN):", missing)

    plot_QC(x.fit) & patchwork::plot_annotation(title=paste0(fit_id, " , fitname: ", fitname), caption=caption)
  })

  pdf(stringr::str_replace_all(fname, ".Rds", ".pdf") %>% stringr::str_replace_all("simul_fit","plots_fit"), width=20, height=16)
  print(plots)
  dev.off()

  pdf(stringr::str_replace_all(fname, ".Rds", ".pdf") %>% stringr::str_replace_all("simul_fit","plots_QC"), width=20, height=16)
  print(QC)
  dev.off()

})







