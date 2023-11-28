py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")

## Simulated ####
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N150_G2_simul.Rds")
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/simulations/fits/fits_dn.matched.2011/simul_fit.N150.G3.s9.matched.2011.Rds")$dataset
input.simul %>% plot_fit()

counts = get_input(input.simul, matrix=T)
# reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),])
reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),],
                     "DBS"=COSMIC_dbs[c("DBS2","DBS5"),])
x.simul = fit(counts=counts,
              k_list=3, # n of denovo signatures
              cluster=6,  # n of clusters
              reference_cat=reference_cat,
              n_steps=1500, lr=0.005,
              hyperparameters=list("scale_factor_centroid"=5000),
              seed_list=c(19,33,2),
              py=py)

x.simul %>% get_initial_object() %>% plot_exposures()
x.simul %>% get_initial_object() %>% plot_centroids()

x.simul %>% plot_exposures()
x.simul %>% plot_centroids()

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


