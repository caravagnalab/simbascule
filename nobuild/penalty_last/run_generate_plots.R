plot_devtools::load_all("~/GitHub/simbasilica/")
load_deps()

run_name = "no_omega"
path = "~/Dropbox/shared/2022. Basilica/simulations/matched_signals/"
files = list.files(paste0(path, "fits/"), pattern=run_name, full.names=TRUE) %>%
  stringr::str_sort(numeric=T)

pdf(paste0(path, "fits_QC.", run_name, ".pdf"), height=8, width=14)
lapply(files, function(fname) {
  fit_i = readRDS(fname)$fit
  tmp = strsplit(fname, split="/")[[1]]
  fit_id = tmp[length(tmp)]
  plot_QC(fit_i) & patchwork::plot_annotation(title=fit_id)
})
dev.off()


pdf(paste0(path, "fits.", run_name, ".pdf"), height=12, width=14)
lapply(files, function(fname) {
  simul_fit_i = readRDS(fname)
  fit_i = simul_fit_i$fit
  simul_i = simul_fit_i$dataset
  tmp = strsplit(fname, split="/")[[1]]
  fit_id = tmp[length(tmp)]
  patchwork::wrap_plots(plot_fit(fit_i),
                        plot_fit(simul_i), ncol=1) &
    patchwork::plot_annotation(title=fit_id,
                               subtitle="Fit (top) and Simulated (bottom)")
})
dev.off()



