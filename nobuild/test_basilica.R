py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")

## Simulated ####
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N150_G2_simul.Rds")
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s9.matched.2011.Rds")
x.simul = input.simul$fit.0
x.simul$clustering = NULL

# counts = get_input(input.simul, matrix=T)
# # reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),])
# reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),],
#                      "DBS"=COSMIC_dbs[c("DBS2","DBS5"),])
# x.simul = fit(counts=counts,
#               k_list=3, # n of denovo signatures
#               # cluster=6,  # n of clusters
#               reference_cat=reference_cat,
#               n_steps=3000, lr=0.005,
#               # hyperparameters=list("scale_factor_centroid"=10000, "scale_factor_alpha"=10000),
#               seed_list=c(19,33,2),
#               py=py)

pi_list = c(0.6); n_clusters = 6; fits = c()
fits = c(fits, lapply(pi_list, function(pi_conc0) {
  x.cl = fit_clustering(x.simul, cluster=n_clusters,
                        n_steps=3000, lr=0.005,
                        # scale_factor_centroid = 1 and scale_factor_alpha = 1000 -> good centroids, bad clustering
                        # scale_factor_centroid = 1 and scale_factor_alpha = 1 -> good centroids, bad clustering
                        hyperparameters=list(
                          "pi_conc0"=pi_conc0,
                          "scale_factor_alpha"=1,
                          "scale_factor_centroid"=1
                          ),
                        seed_list=c(33),
                        py=py)
  return(x.cl)
}) %>% setNames(paste0("pi_",pi_list,"_G_", n_clusters))
)

x.cl = fits$pi_0.6_G_6

x.cl %>% get_initial_object() %>% plot_mixture_weights() %>%
  patchwork::wrap_plots(x.cl %>% plot_mixture_weights(), ncol=1, guides="collect")
x.cl %>% get_initial_object() %>% plot_exposures() %>%
  patchwork::wrap_plots(x.cl %>% plot_exposures(), ncol=1, guides="collect")
x.cl %>% get_initial_object() %>% plot_centroids() %>%
  patchwork::wrap_plots(x.cl %>% plot_centroids(), ncol=2, guides="collect")
x.cl %>% plot_gradient_norms()
x.cl %>% plot_posterior_probs()

get_params(x.cl, what="clustering")

# lapply(fits, function(fitid) {
#   fitid %>% get_initial_object() %>% plot_mixture_weights() %>%
#     patchwork::wrap_plots(fitid %>% plot_mixture_weights(), ncol=1, guides="collect") &
#     patchwork::plot_annotation(title=names(fitid))
# }) %>% patchwork::wrap_plots()

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


