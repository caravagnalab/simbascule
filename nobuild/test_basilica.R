py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")

## Simulated ####
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N150_G2_simul.Rds")
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s9.matched.2011.Rds")
x.simul = input.simul$fit.0
x.simul$clustering = NULL

# counts = get_input(input.simul$dataset, matrix=T, types="SBS")
# reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),])
# x.simul = fit(counts=counts,
#               k_list=3, # n of denovo signatures
#               # cluster=6,  # n of clusters
#               reference_cat=reference_cat,
#               n_steps=3000, lr=0.005,
#               # hyperparameters=list("scale_factor_centroid"=10000, "scale_factor_alpha"=10000),
#               seed_list=c(19,33,2),
#               py=py)

pars = expand.grid(list(scale_factor_alpha=c(1),
                        scale_factor_centroid=c(1),
                        pi_conc0=c(0.006, 0.6, 60)))
n_clusters = 1:6; fits = c()
# fits = c(fits, lapply(1:nrow(pars), function(rowid) {
#   x.cl = fit_clustering(x.simul, cluster=n_clusters,
#                         n_steps=1000, lr=0.005,
#                         # scale_factor_centroid = 1 and scale_factor_alpha = 1000 -> good centroids, bad clustering
#                         # scale_factor_centroid = 1 and scale_factor_alpha = 1 -> good centroids, bad clustering
#                         hyperparameters=list(
#                           "pi_conc0"=pars[rowid, "pi_conc0"],
#                           "scale_factor_alpha"=pars[rowid, "scale_factor_alpha"],
#                           "scale_factor_centroid"=pars[rowid, "scale_factor_centroid"],
#                           "tau"=0
#                           ),
#                         store_parameters=FALSE,
#                         seed_list=c(33),
#                         py=py)
#   return(x.cl)
# }) %>% setNames(paste0("row",1:nrow(pars),"_G", n_clusters))
# )

rowid = 2
x.cl = fit_clustering(x.simul, cluster=1:5,
                      nonparametric=TRUE,
                      n_steps=3000, lr=0.005, # optim_gamma=1e-10,
                      hyperparameters=list("scale_factor_centroid"=1000),
                      store_parameters=FALSE,
                      store_fits=TRUE,
                      seed_list=c(33),
                      py=py)

x.cl %>% plot_fit()
x.cl %>% plot_gradient_norms()
x.cl %>% plot_scores()
x.cl %>% plot_mixture_weights()
x.cl %>% plot_exposures()
x.cl %>% plot_centroids()
x.cl %>% plot_posterior_probs()

x.cl %>% get_initial_object() %>% plot_centroids()

alt_run = get_alternative_run(x.cl, G=3, seed=list("clustering"=33, "nmf"=get_seed(x.cl)[["nmf"]]))
alt_run %>% get_initial_object() %>% plot_centroids() %>%
  patchwork::wrap_plots(plot_centroids(alt_run)) %>% patchwork::wrap_plots(
    alt_run %>% get_initial_object() %>% plot_exposures() %>%
      patchwork::wrap_plots(plot_exposures(alt_run)), ncol=2
  )
alt_run %>% plot_mixture_weights()
alt_run %>% plot_posterior_probs()



x.cl %>% get_initial_object() %>% plot_mixture_weights() %>%
  patchwork::wrap_plots(x.cl %>% plot_mixture_weights(), ncol=1)
x.cl %>% get_initial_object() %>% plot_centroids() %>%
  patchwork::wrap_plots(x.cl %>% plot_centroids(), ncol=1)

x.cl %>% get_initial_object() %>% plot_fit()

pis = lapply(names(fits), function(fitname) {
  rowid = strsplit(fitname, split="_")[[1]][1] %>% stringr::str_replace_all("row","") %>% as.integer()
  titlee = paste0(colnames(pars), " ", pars[rowid,], collapse=" ")
  (plot_mixture_weights(fits[[fitname]] %>% get_initial_object()) + labs(title=titlee)) %>%
    patchwork::wrap_plots(fits[[fitname]] %>% plot_mixture_weights(), ncol=1, guides="collect")
}) %>% patchwork::wrap_plots()


centr = lapply(names(fits), function(fitname) {
  rowid = strsplit(fitname, split="_")[[1]][1] %>% stringr::str_replace_all("row","") %>% as.integer()
  titlee = paste0(colnames(pars), " ", pars[rowid,], collapse=" ")
  (plot_centroids(fits[[fitname]] %>% get_initial_object()) + labs(title=titlee)) %>%
    patchwork::wrap_plots(fits[[fitname]] %>% plot_centroids(), ncol=1, guides="collect")
}) %>% patchwork::wrap_plots(guides="collect")


expos = lapply(names(fits), function(fitname) {
  rowid = strsplit(fitname, split="_")[[1]][1] %>% stringr::str_replace_all("row","") %>% as.integer()
  titlee = paste0(colnames(pars), " ", pars[rowid,], collapse=" ")
  (plot_exposures(fits[[fitname]] %>% get_initial_object()) + labs(title=titlee)) %>%
    patchwork::wrap_plots(fits[[fitname]] %>% plot_exposures(), ncol=1, guides="collect")
}) %>% patchwork::wrap_plots() & theme(legend.position="none")

plot_fit(fits$row2_G6)



# x.cl %>% get_initial_object() %>% plot_mixture_weights() %>%
#   patchwork::wrap_plots(x.cl %>% plot_mixture_weights(), ncol=1, guides="collect")
# x.cl %>% get_initial_object() %>% plot_exposures() %>%
#   patchwork::wrap_plots(x.cl %>% plot_exposures(), ncol=1, guides="collect")
# x.cl %>% get_initial_object() %>% plot_centroids() %>%
#   patchwork::wrap_plots(x.cl %>% plot_centroids(), ncol=2, guides="collect")
# x.cl %>% plot_gradient_norms()
# x.cl %>% plot_posterior_probs()
# get_params(x.cl, what="clustering")
# fits$row1_G6


# pdf("~/Dropbox/shared/2022. Basilica/simulations/clustering_test_scale_factors.pdf", height=15, width=20)
# print(pis)
# print(centr)
# print(expos)
# dev.off()

plot_fit(x.simul)
plot_QC(x.simul)
plot_similarity_reference(x.simul, reference=get_signatures(input.simul, matrix=T)[["SBS"]])


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


