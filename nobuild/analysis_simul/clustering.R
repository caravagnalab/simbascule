devtools::load_all()
load_deps()
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_0507/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/fits_dn.flat_clust.noreg.new_hier.0507/"
save_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"
idd = "N150.G3.s1"

simul = readRDS(paste0(data_path, "simul.", idd, ".Rds")) %>% create_basilica_obj_simul()
fit = readRDS(paste0(fits_path, "fit.", idd, ".Rds")) %>% convert_sigs_names(simul)
fit_clust = readRDS(paste0(fits_path, "fit_clust.", idd, ".Rds")) %>% convert_sigs_names(simul)

fname = ".test_lr005.modsel.1907"

k_list = 3:5  # true = 4
# k_list = 0
clusters_list = 1:5  # true = 2

samples_rare = get_samples_with_sigs(simul, "SBS7c")
cls = get_color_palette(simul)

counts = get_group(simul, c(2,3))
catal = get_signatures(simul)

# fits new_hier enforce_sparsity
alpha_sigma = 0.05
# alpha_p_conc0 = alpha_p_conc1 = 0.5
fit_cl.old.spars = two_steps_inference(x=counts, reg_weight = 0., lr=0.005,
                                      k=1:4,
                                      reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                      clusters=1:3, n_steps=2000,
                                      seed_list=c(10, 33, 92), # reference_catalogue=catal,
                                      # hyperparameters = list("alpha_sigma"=alpha_sigma,
                                      #                        "alpha_p_conc1"=0.6,
                                      #                        "alpha_p_conc0"=0.6),
                                      enforce_sparsity2 = TRUE, py=py, new_hier=FALSE,
                                      nonparametric=FALSE)

fit_cl.old.spars.nonparam = two_steps_inference(x=counts, reg_weight = 0., lr=0.005,
                                       k=1:4,
                                       reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                       clusters=1:5, n_steps=2000,
                                       seed_list=c(10, 33, 92), # reference_catalogue=catal,
                                       hyperparameters = list("alpha_sigma"=alpha_sigma,
                                                              "alpha_p_conc1"=0.6,
                                                              "alpha_p_conc0"=0.6),
                                       enforce_sparsity2 = TRUE, py=py, new_hier=FALSE,
                                       nonparametric=TRUE)

tmp = fit_cl.old.spars.nonparam %>% convert_sigs_names(simul)
(tmp$fit$params$alpha_prior / rowSums(tmp$fit$params$alpha_prior)) %>% as.data.frame() %>%
  tibble::rownames_to_column(var="gid") %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=gid, y=value, fill=variable), stat="identity") +
  scale_fill_manual(values=gen_palette(5))

fit_cl.old.spars.all = two_steps_inference(x=get_data(simul), reg_weight = 0., lr=0.005,
                                       k=k_list,
                                       # reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                       clusters=2:5, n_steps=2000,
                                       seed_list=c(10, 33, 92), reference_catalogue=catal,
                                       hyperparameters = list("alpha_sigma"=alpha_sigma,
                                                              "alpha_p_conc1"=0.6,
                                                              "alpha_p_conc0"=0.6),
                                       enforce_sparsity2 = TRUE, py=py, new_hier=FALSE,
                                       nonparametric=FALSE)

fit_cl.old.spars.all.nonparam = two_steps_inference(x=get_data(simul), reg_weight = 0., lr=0.005,
                                           k=k_list,
                                           # reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                           clusters=2:5, n_steps=2000,
                                           seed_list=c(10, 33, 92), reference_catalogue=catal,
                                           hyperparameters = list("alpha_sigma"=alpha_sigma,
                                                                  "alpha_p_conc1"=0.6,
                                                                  "alpha_p_conc0"=0.6),
                                           enforce_sparsity2 = TRUE, py=py, new_hier=FALSE,
                                           nonparametric=TRUE)

fit_cl.old.nospars = two_steps_inference(x=counts, reg_weight = 0., lr=0.005,
                                         k=k_list,
                                         # reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                         clusters=clusters_list, n_steps=2000,
                                         seed_list=c(10, 33, 92), reference_catalogue=catal,
                                         hyperparameters = list("alpha_sigma"=alpha_sigma),
                                         enforce_sparsity2 = FALSE, py=py, new_hier=FALSE,
                                         nonparametric=FALSE)

fits = readRDS(paste0(save_path, idd, fname, ".Rds"))


## sparsity
aricode::ARI(fit_cl.old.spars$groups, simul$groups[as.integer(get_group(simul, groupIDs = c(2,3), return_idx = T))])
aricode::NMI(fit_cl.old.spars$groups, simul$groups[as.integer(get_group(simul, groupIDs = c(2,3), return_idx = T))])

compute.cosine(m1=get_exposure(fit_cl.old.spars),
               m2=get_exposure(simul)[rownames(counts),],
               assigned_missing=get_assigned_missing(fit_cl.old.spars, simul),
               what="expos", subset_cols=get_sigs_group(simul, groupID=c(2,3)))

compute.mse(m_inf=get_exposure(fit_cl.old.spars),
            m_true=get_exposure(simul)[rownames(counts),],
            subset_cols=get_sigs_group(simul, groupID=c(2,3)))


## sparsity non parametric
aricode::ARI(fit_cl.old.spars.nonparam$groups, simul$groups[as.integer(get_group(simul, groupIDs = c(2,3), return_idx = T))])
aricode::NMI(fit_cl.old.spars.nonparam$groups, simul$groups[as.integer(get_group(simul, groupIDs = c(2,3), return_idx = T))])

compute.cosine(m1=get_exposure(fit_cl.old.spars.nonparam),
               m2=get_exposure(simul)[rownames(counts),],
               assigned_missing=get_assigned_missing(fit_cl.old.spars.nonparam, simul),
               what="expos", subset_cols=get_sigs_group(simul, groupID=c(2,3)))

compute.mse(m_inf=get_exposure(fit_cl.old.spars.nonparam),
            m_true=get_exposure(simul)[rownames(counts),],
            subset_cols=get_sigs_group(simul, groupID=c(2,3)))


## no sparsity
aricode::ARI(fit_cl.old.nospars$groups, simul$groups[as.integer(get_group(simul, groupIDs = c(2,3), return_idx = T))])
aricode::NMI(fit_cl.old.nospars$groups, simul$groups[as.integer(get_group(simul, groupIDs = c(2,3), return_idx = T))])

compute.cosine(m1=get_exposure(fit_cl.old.nospars),
               m2=get_exposure(simul)[rownames(counts),],
               assigned_missing=get_assigned_missing(fit_cl.old.nospars, simul),
               what="expos", subset_cols=get_sigs_group(simul, groupID=c(2,3)))

compute.mse(m_inf=get_exposure(fit_cl.old.nospars),
            m_true=get_exposure(simul)[rownames(counts),],
            subset_cols=get_sigs_group(simul, groupID=c(2,3)))


## sparsity all
aricode::ARI(fit_cl.old.spars.all$groups, simul$groups)
aricode::NMI(fit_cl.old.spars.all$groups, simul$groups)

aricode::ARI(fit_cl.old.spars.all$groups, get_groups_rare(simul))
aricode::NMI(fit_cl.old.spars.all$groups, get_groups_rare(simul))

compute.cosine(m1=get_exposure(fit_cl.old.spars.all), m2=get_exposure(simul),
               assigned_missing=get_assigned_missing(fit_cl.old.spars.all, simul),
               what="expos")

compute.mse(m_inf=get_exposure(fit_cl.old.spars.all), m_true=get_exposure(simul))


## sparsity all non parametric
aricode::ARI(fit_cl.old.spars.all.nonparam$groups, simul$groups)
aricode::NMI(fit_cl.old.spars.all.nonparam$groups, simul$groups)

aricode::ARI(fit_cl.old.spars.all.nonparam$groups, get_groups_rare(simul))
aricode::NMI(fit_cl.old.spars.all.nonparam$groups, get_groups_rare(simul))

compute.cosine(m1=get_exposure(fit_cl.old.spars.all.nonparam),
               m2=get_exposure(simul),
               assigned_missing=get_assigned_missing(fit_cl.old.spars.all.nonparam, simul),
               what="expos")

compute.mse(m_inf=get_exposure(fit_cl.old.spars.all.nonparam),
            m_true=get_exposure(simul))


saveRDS(list("fit_cl.old.spars"=fit_cl.old.spars,
             "fit_cl.old.spars.nonparam"=fit_cl.old.spars.nonparam,
             "fit_cl.old.spars.all"=fit_cl.old.spars.all,
             "fit_cl.old.spars.all.nonparam"=fit_cl.old.spars.all.nonparam,
             "fit_cl.old.nospars"=fit_cl.old.nospars,
             "simul"=simul), paste0(save_path, idd, fname, ".Rds"))



fit_cl2 %>% plot_gradient_norms()
fit_cl2 %>% plot_posterior_probs()
fit_cl2 %>% convert_sigs_names(simul) %>% plot_exposures()

fit_cl %>% convert_sigs_names(simul) %>% plot_exposures(plot_noise = T)
fit_cl %>% convert_sigs_names(simul) %>% plot_exposures(add_centroid = T)


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


