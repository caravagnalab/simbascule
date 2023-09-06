devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")
# load_deps(base_path="~/GitHub/")

## Simulated ####
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N500_G3_simul.Rds")
input.simul %>% create_basilica_obj_simul()  # creates a basilica obj with simulated data

x.simul = fit(x=input.simul$counts[[1]],
              k=2:6, # n of denovo signatures
              clusters=6,  # n of clusters
              n_steps=2000, lr=0.005,
              enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
              dirichlet_prior=FALSE,
              reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
              hyperparameters=list("alpha_sigma"=0.05, "scale_factor"=1000),  # change default values to hyperparameters
              nonparametric=TRUE, py=py)


## Real data ####
input.real = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N500_CRC.Rds")

x.real = fit(x=input.real$counts[[1]],
             k=15, # n of denovo signatures
             clusters=7,  # n of clusters
             n_steps=2000, lr=0.005,
             enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
             hyperparameters=list("alpha_sigma"=0.05, "scale_factor"=1000),  # change default values to hyperparameters
             nonparametric=TRUE)


## functions
get_obj_initial_params(x.simul)  # initial conditions
get_group(x.simul, groupIDs=c("1"), return_idx=T)  # sample names belonging to grp 1

convert_sigs_names(x.simul, reference_cat=COSMIC_filt)  # assigns denovo to catalogue

plot_exposures(x.simul)
plot_signatures(x.simul)
plot_mutations(x.simul)
plot_posterior_probs(x.simul)  # heatmap with posterior probs
plot_scores(x.simul)
plot_gradient_norms(x.simul)







