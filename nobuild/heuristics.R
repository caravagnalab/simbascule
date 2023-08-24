devtools::load_all()
devtools::load_all("~/GitHub/simbasilica/")
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_3107/"
fits_path = "~/Dropbox/shared/2022. Basilica/simulations/fits_dn.clust.nonparametric.sparsity.noreg.old_hier.0208/"
idd = "N150.G3.s30.1"
simul = readRDS(paste0(data_path, "simul.", idd, ".Rds")) %>% create_basilica_obj_simul()
x = readRDS(paste0(fits_path, "fit_clust.", idd, ".Rds")) %>% convert_sigs_names(simul)

## EM to recluster ####
x %>% fix_assignments() %>% plot_exposures()



## Linear combination ####
sign1 = get_denovo_signatures(x)
sign2 = COSMIC_filt

filter_signatures_QP(sign1, sign2, return_weights=TRUE)
# filter_signatures_QP(sign2, sign1, return_weights=TRUE)

x %>% plot_signatures(catalogue = get_signatures(simul))


## Compare with K-Means on exposures ####
x_init = get_obj_initial_params(x)


aricode::ARI(fix_assignments(x)$groups, simul$groups)
aricode::ARI(x_init$groups, simul$groups)
aricode::ARI(fix_assignments(x_init)$groups, simul$groups)




