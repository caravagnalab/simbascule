devtools::load_all()
load_deps()
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_0507/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/fits_dn.flat_clust.noreg.new_hier.0507/"
save_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"
idd = "N150.G3.s1"

simul = readRDS(paste0(data_path, "simul.", idd, ".Rds")) %>% create_basilica_obj_simul()
fit = readRDS(paste0(fits_path, "fit.", idd, ".Rds")) %>% convert_sigs_names(simul)
fit_clust = readRDS(paste0(fits_path, "fit_clust.", idd, ".Rds")) %>% convert_sigs_names(simul)

fname = ".test_lr005.modsel.1207"

k_list = 4 # 2:6  # true = 4
clusters_list = 4 # 3:5  # true = 4

samples_rare = get_samples_with_sigs(simul, "SBS7c")
cls = get_color_palette(simul)

fit2 = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 0., lr=0.005,
                           k=fit$k_list, seed_list=c(10, 33, 92), n_steps=1500,
                           reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                           # hyperparameters = list("alpha_noise_var"=0.005),
                           py=py, new_hier=FALSE, nonparametric=FALSE)

fit_cl = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 0., lr=0.005,
                             k=k_list, clusters=clusters_list, seed_list=c(10, 33, 92), n_steps=1500,
                             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                             # hyperparameters = list("alpha_noise_var"=0.005),
                             py=py, new_hier=TRUE, nonparametric=FALSE)

fit_cl.nonpar = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 0., lr=0.005,
                                    k=k_list, clusters=max(clusters_list), seed_list=c(10, 33, 92), n_steps=1500,
                                    reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                    # hyperparameters = list("alpha_noise_var"=0.005),
                                    py=py, new_hier=TRUE, nonparametric=TRUE)

saveRDS(list("param"=fit_cl, "nonparam"=fit_cl.nonpar, "simul"=simul), paste0(save_path, idd, fname, ".Rds"))


fit_cl %>% plot_gradient_norms()
fit_cl %>% plot_posterior_probs()
fit_cl %>% convert_sigs_names(simul) %>% plot_exposures(cls=cls, sample_name=T)

fit_cl %>% convert_sigs_names(simul, cutoff=0.6) %>% plot_exposures(plot_noise = T)
fit_cl %>% convert_sigs_names(simul, cutoff=0.6) %>% plot_exposures(plot_noise = T, add_centroid = T, sampleIDs = 1)

# exposures_noise_centroid.fit_cl = lapply(unique(fit_cl$groups), function(gid) {
#   p1 = fit_cl %>% convert_sigs_names(simul) %>%
#     plot_exposures(sampleIDs=get_group(fit_cl, gid, return_idx=T), cls=cls, sample_name=T, add_centroid=T)
#   p2 = fit_cl %>% convert_sigs_names(simul) %>%
#     plot_exposures(sampleIDs=get_group(fit_cl, gid, return_idx=T), cls=cls, sample_name=T, plot_noise=T, add_centroid=T)
#
#   return(patchwork::wrap_plots(p1, p2, ncol=1, guides="collect"))
#   }
# ) %>% patchwork::wrap_plots(guides="collect")

expos_fit_cl = plot_exposures(fit_cl %>% convert_sigs_names(simul), cls=cls, add_centroid=T, sample_name=T) %>%
  patchwork::wrap_plots(plot_exposures(fit_cl %>% convert_sigs_names(simul), cls=cls, plot_noise=T, add_centroid=T, sample_name=T),
                        guides="collect", ncol=1)

expos_fit_cl.nonpar = plot_exposures(fit_cl.nonpar %>% convert_sigs_names(simul), cls=cls, add_centroid=T, sample_name=T) %>%
  patchwork::wrap_plots(plot_exposures(fit_cl.nonpar %>% convert_sigs_names(simul), cls=cls, plot_noise=T, add_centroid=T, sample_name=T),
                        guides="collect", ncol=1)



pdf(paste0(save_path, "plots.", idd, fname, ".pdf"), height=8, width=12)

expos_fit_cl &
  patchwork::plot_annotation(title="Parametric",
                             subtitle=paste("MSE expos =",
                                            compute.mse(get_exposure(fit_cl), get_exposure(simul)),
                                            "- Cosine expos =",
                                            compute.cosine(get_exposure(fit_cl), get_exposure(simul), what="expos",
                                                           assigned_missing=get_assigned_missing(fit_cl, simul)),
                                            "- ARI =", aricode::ARI(fit_cl$groups, get_groups_rare(simul, fit_cl)),
                                            "- NMI =", aricode::NMI(fit_cl$groups, get_groups_rare(simul, fit_cl))))
expos_fit_cl.nonpar &
  patchwork::plot_annotation(title="Non-parametric",
                             subtitle=paste("MSE expos =",
                                            compute.mse(get_exposure(fit_cl.nonpar), get_exposure(simul)),
                                            "- Cosine expos =",
                                            compute.cosine(get_exposure(fit_cl.nonpar), get_exposure(simul), what="expos",
                                                           assigned_missing=get_assigned_missing(fit_cl.nonpar, simul)),
                                            "- ARI =", aricode::ARI(fit_cl.nonpar$groups, get_groups_rare(simul, fit_cl.nonpar)),
                                            "- NMI =", aricode::NMI(fit_cl.nonpar$groups, get_groups_rare(simul, fit_cl.nonpar))))

plot_fit(fit_cl %>% convert_sigs_names(simul), simul, cls=cls) & patchwork::plot_annotation(title="Parametric")
plot_fit(fit_cl.nonpar %>% convert_sigs_names(simul), simul, cls=cls) & patchwork::plot_annotation(title="Non-parametric")

dev.off()



fit_cl %>% plot_scores()


