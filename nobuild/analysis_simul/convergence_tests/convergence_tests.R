x_clust = readRDS("~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.noreg.old_hier/fit_clust.N1000.G4.s10.Rds")
x_hier = readRDS("~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.noreg.old_hier/fit_hier.N1000.G4.s10.Rds")
x_flat = readRDS("~/GitHub/simbasilica/nobuild/simulations/fits_dn_1606.noreg.old_hier/fit.N1000.G4.s10.Rds")
x_simul = readRDS("~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/simul.N1000.G4.s10.Rds") %>%
  create_basilica_obj_simul()

counts = get_data(x_simul)
true_k = length(get_signames(x_simul))
true_dn = length(get_dn_signames(x_simul))
min_k = max(0, true_dn-3); max_k = true_dn+3

min_gr = max(0, length(x_simul$groups %>% unique()) - 3); max_gr = length(x_simul$groups %>% unique()) + 3


load_deps()

## lower lr #####
x_clust.lr = two_steps_inference(x=counts,
                                 k=min_k:max_k,
                                 reference_catalogue=COSMIC_filt_merged,
                                 lr=0.005,
                                 py=py,
                                 filtered_catalogue=TRUE,
                                 clusters=min_gr:max_gr,

                                 reg_weight=0., seed_list=(10,))

