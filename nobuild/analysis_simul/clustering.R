devtools::load_all()
load_deps()
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_0507/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/fits_dn.flat_clust.noreg.new_hier.0507/"
idd = "N150.G3.s1"

simul = readRDS(paste0(data_path, "simul.", idd, ".Rds")) %>% create_basilica_obj_simul()
fit = readRDS(paste0(fits_path, "fit.", idd, ".Rds"))
fit_clust = readRDS(paste0(fits_path, "fit_clust.", idd, ".Rds"))

k_list = fit_clust$k_list
clusters_list = 1:5

samples_rare = get_samples_with_sigs(simul, "SBS7c")
cls = get_color_palette(simul)

fit_cl = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 0.,
                             k=k_list, clusters=clusters_list, seed_list=c(10, 33, 92), n_steps=500,
                             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                             py=py, new_hier=TRUE)

fit_cl$fit$post_probs

pheatmap::pheatmap(fit_cl$fit$post_probs, cluster_rows=T, cluster_cols=F)
fit_cl$fit$gradient_norms %>% as.data.frame() %>%
  tibble::rownames_to_column(var="step") %>%
  reshape2::melt(id="step", variable.name="parameter", value.name="value") %>%
  ggplot() + geom_line(aes(x=as.integer(step), y=value)) + facet_wrap(~parameter, scales="free_y") + theme_bw()
fit_cl$fit$gradient_norms$pi_param %>% plot()
fit_cl$fit$gradient_norms$alpha_t_param %>% plot()


fit_cl %>% convert_sigs_names(simul) %>% plot_exposures(sort_by="SBS7c", cls=cls, sample_name=T)
fit_cl %>% convert_sigs_names(simul) %>% plot_exposures(sort_by="SBS7c", sampleIDs=samples_rare, cls=cls,
                                                        sample_name=T)
fit_cl %>% convert_sigs_names(simul) %>% plot_exposures(sampleIDs=rownames(get_group(fit_cl, 0)), cls=cls,
                                                        sample_name=T)

simul %>% plot_exposures(sort_by="SBS7c", cls=cls, sampleIDs=samples_rare, sample_name=T)


# fit.cl %>% convert_sigs_names(simul) %>%
#   plot_fit(x.true=simul,
#            cls=gen_palette(n=6) %>%
#              setNames(get_signames(simul)),
#            sample_name=T, sampleIDs=5)

