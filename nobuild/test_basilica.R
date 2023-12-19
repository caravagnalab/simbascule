reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")


## Run fits on small test ####
fits_dir = "~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/"
fitsname = list.files(fits_dir, pattern="s12.matched.2011.Rds")


## fits ####
# fits_cl = lapply(fitsname, function(fname) {
#   G = strsplit(fname, "[.]")[[1]][3] %>% stringr::str_replace_all("G","") %>% as.integer()
#   x.nmf = readRDS(paste0(fits_dir, fname))$fit.0
#   x.cl = fit_clustering(x.nmf, cluster=G * 2,
#                         n_steps=3000, lr=0.005,
#                         nonparametric=TRUE,
#                         store_parameters=FALSE,
#                         seed_list=c(11,33,4392),
#                         py=py)
#   return(x.cl)
#   }) %>% setNames(fitsname)


## plots ####
save_path = "~/Dropbox/shared/2022. Basilica/simulations/matched_signals/fits_cl_test/"
fitsname = list.files(save_path, pattern=".Rds") %>% gtools::mixedsort()

## save objects
# lapply(fitsname, function(fname) {
#   simul_obj = readRDS(paste0(fits_dir, fname))
#   simul_obj$fit.0.cl = fits_cl[[fname]]
#   saveRDS(object=simul_obj, file=paste0(save_path, fname))
# })


lapply(fitsname, function(fname) {
  simul_fit = readRDS(paste0(save_path, fname))
  x = convert_dn_names(simul_fit$fit.0.cl %>% filter_denovo(), simul_fit$dataset)

  tt = table(simul_fit$dataset$clustering$clusters$clusters)
  true_pis = (tt / sum(tt)) %>% as.numeric()
  simul_fit$dataset$clustering[["pyro"]][["params"]][["infered_params"]][["pi"]] = true_pis

  fname_tmp = stringr::str_replace_all(fname, ".Rds", "")
  pdf(paste0(save_path, stringr::str_replace_all(fname_tmp, "simul_fit", "plots"), ".pdf"),
      width=14, height=14)
  print(plot_fit(x) & patchwork::plot_annotation(title="Fit"))
  print(plot_fit(simul_fit$dataset) & patchwork::plot_annotation(title="Simulated"))
  print(plot_posterior_probs(x))
  print(plot_QC(x))
  dev.off()
})


pdf(paste0(save_path, "fits_clustering.pdf"), width=18, height=18)
p_tot = lapply(fitsname, function(fname) {
  simul_fit = readRDS(paste0(save_path, fname))

  tt = table(simul_fit$dataset$clustering$clusters$clusters)
  true_pis = (tt / sum(tt)) %>% as.numeric()
  simul_fit$dataset$clustering[["pyro"]][["params"]][["infered_params"]][["pi"]] = true_pis

  x = convert_dn_names(simul_fit$fit.0.cl %>% filter_denovo(), simul_fit$dataset)

  p_fit = plot_fit(x)
  p_simul = plot_fit(simul_fit$dataset)

  patchwork::wrap_plots(p_fit, p_simul, ncol=1) &
    patchwork::plot_annotation(title=stringr::str_replace_all(fname,"simul_fit.|.matched.2011.Rds",""),
                               subtitle="Fitted (top) and simulated (bottom) data.")
  })
p_tot
dev.off()




## example ####
# input.simul = readRDS("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s9.matched.2011.Rds")
input_simul = readRDS("/Users/elenab/Dropbox/shared/2022. Basilica/simulations/matched_signals/fits_cl_test/simul_fit.N150.G3.s12.matched.2011.Rds")

x.simul = input_simul$dataset
x.cl = input_simul$fit.0.cl





## Real data ####
input.real = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N500_CRC.Rds")

counts = list("SBS"=input.real$counts[[1]][,colnames(COSMIC_filt)])
reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),])

x.real = fit(counts=counts,
             k_list=0:10, # n of denovo signatures
             # cluster=7,  # n of clusters
             reference_cat=reference_cat,
             n_steps=2000, lr=0.005,
             seed_list=c(19,33,2),
             py=py)

plot_fit(x.real)
plot_QC(x.real)
plot_similarity_reference(x.real, reference=COSMIC_filt)


