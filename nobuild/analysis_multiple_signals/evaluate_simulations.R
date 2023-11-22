devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")

source("~/GitHub/simbasilica/nobuild/analysis_multiple_signals/evaluate_aux_fns.R")

run_id = "matched.2011"
path = paste0("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.", run_id, "/")
files = list.files(path, full.names=T, pattern=".Rds")

all_stats = lapply(files, function(fname) {
  stats_single_data(fname)
}) %>% dplyr::bind_rows()

saveRDS(all_stats, "~/Dropbox/shared/2022. Basilica/simulations/stats_", run_id, ".Rds")

stats_plots = lapply(c("NoPenalty","PenaltyN"), function(penalty_id) {
  stats_p = all_stats %>% dplyr::filter(penalty==penalty_id)
  make_plots_stats(stats_p)
}) %>% setNames(c("NoPenalty","PenaltyN"))

pdf("~/Dropbox/shared/2022. Basilica/simulations/stats_", run_id, ".pdf", height=10, width=20)
patchwork::wrap_plots(stats_plots) & patchwork::plot_annotation(title="NoPenalty and PenaltyN")
dev.off()


lapply(files, function(fname) {
  tmp = strsplit(fname, split="/")[[1]]; fit_id = tmp[length(tmp)]
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  types = get_types(x.simul)

  plots = lapply(c("fit.0","fit.N"), function(fitname) {
    x.fit = simul_fit[[fitname]] %>% rename_dn_expos()
    assigned_missing = lapply(types, function(tid) {
      get_assigned_missing(x.fit=x.fit, x.simul=x.simul, type=tid)
    }) %>% setNames(types)

    added = sapply(assigned_missing, function(i) i[["added_fp"]]) %>% unlist() %>% paste(collapse=",")
    missing = sapply(assigned_missing, function(i) i[["missing_fn"]]) %>% unlist() %>% paste(collapse=",")

    caption = paste("Added signatures (FP):", added, "- Missing signatures (FN):", missing)

    patchwork::wrap_plots(plot_fit(x.fit),
                          plot_fit(x.simul), ncol=1) &
      patchwork::plot_annotation(title=paste0(fit_id, " , fitname: ", fitname),
                                 subtitle="Fit (top) and simulated (bottom)",
                                 caption=caption)
  })

  QC = lapply(c("fit.0","fit.N"), function(fitname) {
    x.fit = simul_fit[[fitname]] %>% rename_dn_expos()
    assigned_missing = lapply(types, function(tid) {
      get_assigned_missing(x.fit=x.fit, x.simul=x.simul, type=tid)
    }) %>% setNames(types)

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







