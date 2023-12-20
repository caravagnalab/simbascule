reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")


## Run fits on small test ####
fits_dir = "~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/"
fitsname = list.files(fits_dir, pattern="s12.matched.2011.Rds")


## fits ####
fits_cl = lapply(fitsname, function(fname) {
  G = strsplit(fname, "[.]")[[1]][3] %>% stringr::str_replace_all("G","") %>% as.integer()
  x.nmf = readRDS(paste0(fits_dir, fname))$fit.0
  x.cl = fit_clustering(x.nmf, cluster=G * 2,
                        n_steps=3000, lr=0.005,
                        nonparametric=TRUE,
                        store_parameters=FALSE,
                        store_fits=TRUE,
                        seed_list=c(11,33,4392),
                        py=py)
  return(x.cl)
  }) %>% setNames(fitsname)


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
input_simul = readRDS("/Users/elenab/Dropbox/shared/2022. Basilica/simulations/matched_signals/fits_cl_test/simul_fit.N150.G6.s12.matched.2011.Rds")

x.simul = input_simul$dataset
x.cl = input_simul$fit.0.cl %>% filter_denovo() %>% convert_dn_names(x.simul)
# x.cl$clustering = NULL

x.cl %>% plot_exposures()

counts = get_input(x.simul, matrix=T)
reference = get_fixed_signatures(x.cl, matrix=T)

x.cl_new = fit_clustering(x=x.cl, cluster=12, py=py)
x.cl_new %>% plot_exposures()



samples = get_input(x.cl, clusters="G2")$SBS$samples %>% unique()
samples = samples[which(get_exposure(x.cl, matrix=T)[["SBS"]][samples, "SBS18"] > 0.2)]
plot_exposures(x.cl, samples = samples, sample_name = T)
plot_centroids(x.cl)

centr = get_centroids(x.cl, matrix=T)
dbs = get_signames(x.cl)$DBS

exp1 = get_exposure(x.cl, matrix=T)$DBS["G4_3", ] %>% as.numeric()
centrg2 = centr[c("G2"), paste0("1_",dbs)] %>% as.numeric()
centrg3 = centr[c("G3"), paste0("1_",dbs)] %>% as.numeric()

gtools::ddirichlet(exp1, centrg2)
gtools::ddirichlet(exp1, centrg3)


shifter = function(par) {
  par = par[par >= 0.01]
  return(par / sum(par))
}

exp1_b = shifter(exp1)
centrg2_b = shifter(centrg2)
centrg3_b = shifter(centrg3)

gtools::ddirichlet(exp1_b, centrg2_b[1:3])
gtools::ddirichlet(exp1_b, centrg3_b[1:3])


get_pyro_stat(x.cl, what="clustering",statname="params")[[1]]$infered_params$post_probs[samples,]

x.cl %>%
  get_initial_object() %>%
  plot_exposures() %>%
  patchwork::wrap_plots(x.cl %>% plot_exposures(), guides="collect")

x.cl %>%
  get_initial_object() %>%
  plot_centroids() %>%
  patchwork::wrap_plots(x.cl %>% plot_centroids(), guides="collect")

x.cl %>% plot_gradient_norms()


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


