py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")

## Simulated ####
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N150_G2_simul.Rds")
input.simul %>% plot_fit()

counts = get_input(input.simul, matrix=T)
reference_cat = list("SBS"=COSMIC_filt[c("SBS1","SBS5"),])
x.simul = fit(counts=counts,
              k_list=0:4, # n of denovo signatures
              cluster=5,  # n of clusters
              reference_cat=reference_cat,
              n_steps=2000, lr=0.005,
              seed_list=c(19,33,2),
              py=py)

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


