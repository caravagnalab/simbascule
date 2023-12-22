reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")


## Run fits on small test ####
fits_dir = "~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/"
save_path = "~/Dropbox/shared/2022. Basilica/simulations/matched_signals/fits_cl_test/"
fitsname = list.files(fits_dir, pattern="s12.matched.2011.Rds")

## fits ####
fits_cl_autoguide = lapply(fitsname, function(fname) {
  G = strsplit(fname, "[.]")[[1]][3] %>% stringr::str_replace_all("G","") %>% as.integer()
  x.nmf = readRDS(paste0(fits_dir, fname))$fit.0
  x.cl = fit_clustering(x.nmf, cluster=G * 2,
                        n_steps=3000, lr=0.005,
                        nonparametric=TRUE,
                        autoguide = TRUE,
                        store_parameters=FALSE,
                        store_fits=TRUE,
                        seed_list=c(11,33,4392),
                        py=py)
  return(x.cl)
  }) %>% setNames(fitsname)


fits_cl_autoguide$simul_fit.N1000.G3.s12.matched.2011.Rds %>% plot_exposures()


# fits_cl_manguide = lapply(fitsname, function(fname) {
#   G = strsplit(fname, "[.]")[[1]][3] %>% stringr::str_replace_all("G","") %>% as.integer()
#   x.nmf = readRDS(paste0(fits_dir, fname))$fit.0
#   x.cl = fit_clustering(x.nmf, cluster=G * 2,
#                         n_steps=3000, lr=0.005,
#                         nonparametric=TRUE,
#                         autoguide = FALSE,
#                         store_parameters=FALSE,
#                         store_fits=TRUE,
#                         seed_list=c(11,33,4392),
#                         py=py)
#   return(x.cl)
# }) %>% setNames(fitsname)


## save objects
lapply(fitsname, function(fname) {
  simul_obj = readRDS(paste0(fits_dir, fname))
  simul_obj$fit.0.cl = fits_cl_autoguide[[fname]]
  saveRDS(object=simul_obj, file=paste0(save_path, "autoguide_", fname))
})

lapply(fitsname, function(fname) {
  simul_obj = readRDS(paste0(fits_dir, fname))
  simul_obj$fit.0.cl = fits_cl_manguide[[fname]]
  saveRDS(object=simul_obj, file=paste0(save_path, "manguide_", fname))
})


## plots ####
fitsname = list.files(save_path, pattern=".Rds") %>% gtools::mixedsort()

lapply(fitsname, function(fname) {
  simul_fit = readRDS(paste0(save_path, fname))
  x = convert_dn_names(simul_fit$fit.0.cl %>% filter_denovo(), simul_fit$dataset)

  tt = table(simul_fit$dataset$clustering$clusters$clusters)
  true_pis = (tt / sum(tt)) %>% as.numeric()
  simul_fit$dataset$clustering[["pyro"]][["params"]][["infered_params"]][["pi"]] = true_pis

  fname_tmp = stringr::str_replace_all(fname, ".Rds", "")
  pdf(paste0(save_path,
             stringr::str_replace_all(fname_tmp, "simul_fit", "plots"),
             ".pdf"),
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
reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")

input_simul = readRDS("/Users/elenab/Dropbox/shared/2022. Basilica/simulations/matched_signals/fits_cl_test/autoguide_simul_fit.N150.G3.s12.matched.2011.Rds")

x.simul = input_simul$dataset
x.cl = input_simul$fit.0 # %>% filter_denovo() %>% convert_dn_names(x.simul)
# x.cl$clustering = NULL

x.cl.autoguide = fit_clustering(x=x.cl, cluster=6, lr=0.005,
                                autoguide=TRUE,
                                hyperparameters=list("z_tau"=10),
                                nonparametric=TRUE, py=py,
                                store_fits=TRUE,
                                seed_list=c(11,33,382,48329))

x.cl.manguide = fit_clustering(x=x.cl, cluster=6, lr=0.005,
                                autoguide=FALSE,
                                hyperparameters=list("z_tau"=10),
                                nonparametric=TRUE, py=py,
                                store_fits=TRUE,
                                seed_list=c(11,33,382,48329))

plot_exposures(x.cl.autoguide) %>% patchwork::wrap_plots(plot_exposures(x.cl.manguide))
plot_centroids(x.cl.autoguide) %>% patchwork::wrap_plots(plot_centroids(x.cl.manguide))
plot_mixture_weights(x.cl.autoguide) %>% patchwork::wrap_plots(plot_mixture_weights(x.cl.manguide))

plot_gradient_norms(x.cl.autoguide)
plot_gradient_norms(x.cl.manguide)

plot_exposures(x.cl.autoguide, add_centroid = T)
plot_posterior_probs(x.cl.autoguide)

plot_scores(x.cl.autoguide)
plot_scores(input_simul$fit.0.cl)


counts = get_input(x.simul, matrix=T)
reference = get_fixed_signatures(x.cl, matrix=T)


x = x.cl.autoguide
samples = get_input(x, clusters="G1")$SBS$samples %>% unique()
samples = samples[which(get_exposure(x, matrix=T)[["SBS"]][samples, "SBSD3"] > 0.2)]
plot_exposures(x, samples = samples, sample_name = T, add_centroid=T)
plot_posterior_probs(x, samples=samples)

sid = "G1_11"
centr = get_centroids(x, matrix=T)
sbs = get_signames(x)$SBS
dbs = get_signames(x)$DBS

exp1 = get_exposure(x, matrix=T)$SBS[sid, ] %>% as.numeric()
centrg0 = centr[c("G0"), paste0("0_",sbs)] %>% as.numeric()
centrg2 = centr[c("G1"), paste0("0_",sbs)] %>% as.numeric()

gtools::ddirichlet(exp1, centrg0)
gtools::ddirichlet(exp1, centrg2)


get_idxs = function(par) {
  return(which(par >= 0.01))
}

shifter = function(par) {
  par[par < 0.01] = 1e-10
  return(par / sum(par))
}

idxs = get_idxs(exp1)
exp1_b = exp1[idxs] / sum(exp1[idxs])
centrg0_b = centrg0[idxs] / sum(centrg0[idxs])
centrg2_b = centrg2[idxs] / sum(centrg2[idxs])

gtools::ddirichlet(exp1_b, centrg0_b*100)
gtools::ddirichlet(exp1_b, centrg2_b*100)


gtools::ddirichlet(shifter(exp1), shifter(centrg0)*100)
gtools::ddirichlet(shifter(exp1), shifter(centrg2)*100)


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


