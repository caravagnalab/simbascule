devtools::load_all("~/GitHub/simbasilica/")
devtools::load_all("~/GitHub/basilica/")

source("~/GitHub/simbasilica/nobuild/analysis_multiple_signals/evaluate_aux_fns.R")

path = "~/Dropbox/shared/2022. Basilica/simulations/fits/last_model/fits_dn.matched.2011/"
files = list.files(path, full.names=T)

all_stats = lapply(files[1:10], function(fname) {
  stats_single_data(fname)
}) %>% dplyr::bind_rows()

saveRDS(all_stats, "~/Dropbox/shared/2022. Basilica/simulations/last_model_stats.Rds")

stats_plots = lapply(c("NoPenalty","PenaltyN"), function(penalty_id) {
  stats_p = all_stats %>% dplyr::filter(penalty==penalty_id)
  make_plots_stats(stats_p)
}) %>% setNames(c("NoPenalty","PenaltyN"))


pdf("~/Dropbox/shared/2022. Basilica/simulations/last_model_stats.pdf", height=10, width=20)
patchwork::wrap_plots(stats_plots) & patchwork::plot_annotation(title="NoPenalty and PenaltyN")
dev.off()


lapply(files, function(fname) {
  tmp = strsplit(fname, split="/")[[1]]; fit_id = tmp[length(tmp)]
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  x.fit.0 = simul_fit$fit.0
  x.fit.N = simul_fit$fit.N

  plots.0 = patchwork::wrap_plots(plot_fit(x.fit.0),
                                  plot_fit(x.simul), ncol=1) &
    patchwork::plot_annotation(title=fit_id,
                               subtitle="Fit no penalty (top) and Simulated (bottom)")

  })







