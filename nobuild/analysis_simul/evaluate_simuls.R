devtools::load_all()
load_deps()

main_path = "~/Github/simbasilica/nobuild/simulations/"

save_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"

data_path = paste0(main_path, "synthetic_datasets_2507/")
fits_path = c(paste0(main_path, "fits_dn.flat_clust.parametric.noreg.old_hier.2507/"),
               paste0(main_path, "fits_dn.flat_clust.nonparametric.noreg.old_hier.2507/"))

run_id = c("noreg.old_hier.param", "noreg.old_hier.nonparam") %>% setNames(fits_path)


cutoff = 0.8; df_id = "2507"
stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=cutoff, fits_pattern=c("fit_clust.")) %>%
  dplyr::mutate(run_id=run_id[fits_path],
                unique_id=paste(inf_type, run_id, sep=".")) %>%

  dplyr::add_row(
    get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=cutoff, fits_pattern=c("fit.")) %>%
        dplyr::mutate(run_id=run_id[fits_path])
  ) %>%
  dplyr::mutate(clust_type=dplyr::case_when(
    grepl(".nonparam", run_id) ~ "non-parametric",
    grepl(".param", run_id) ~ "parametric",
    .default="flat")
  )
saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))


stats_df %>%
  dplyr::mutate(rare_ratio=n_rare_found/n_rare) %>%
  ggplot() +
  geom_jitter(aes(x=rare_freq, y=rare_ratio, color=inf_type), height=0.01, width=0.01) +
  ylim(-0.01,1+0.01) + xlim(-0.01,1+0.01) +
  facet_grid(G~clust_type) +
  theme_bw()



pdf(paste0(save_path, "stats_report.sim", cutoff*100, ".", df_id, ".pdf"), height=6, width=10)

plot_sigs_clusters_found(stats_df, what="all", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
# plot_sigs_clusters_found(stats_df, what="common", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
plot_sigs_clusters_found(stats_df, what="rare", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))

# plot_sigs_clusters_found(stats_df, what="all", facet=T, ratio=F, scales="free_y")
# plot_sigs_clusters_found(stats_df, what="common", facet=T, ratio=F, scales="free_y")
# plot_sigs_clusters_found(stats_df, what="rare", facet=T, ratio=F, scales="free_y")

plot_mse_cosine(stats_df, "mse_counts", scales="free_y", ylim=c(0,NA))
plot_mse_cosine(stats_df, "mse_expos", scales="free_y", ylim=c(0,NA))
plot_mse_cosine(stats_df, "mse_expos_rare", scales="free_y", ylim=c(0,NA))

plot_mse_cosine(stats_df, "cosine_sigs", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df, "cosine_expos", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df, "cosine_expos_rare", scales="free_y", ylim=c(0,1))

plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"), what="clusters", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
# plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"), what="clusters", facet=T, ratio=F, scales="free_y")

plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "nmi", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "nmi_rare", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "ari", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "ari_rare", scales="free_y", ylim=c(0,1))

dev.off()


idd = (stats_df %>% dplyr::filter(G==1, N==150) %>% dplyr::arrange(nmi_rare) %>% dplyr::pull(idd))[1]
fits_path = (stats_df %>% dplyr::filter(G==1, N==150) %>% dplyr::arrange(nmi_rare) %>%
               dplyr::pull(fits_path))[1]
simul = readRDS(paste0(data_path, "simul.", idd, ".Rds")) %>% create_basilica_obj_simul()
fit_cl = readRDS(paste0(fits_path, "fit_clust.", idd, ".Rds")) %>% convert_sigs_names(simul)

aricode::ARI(fit_cl$groups, get_groups_rare(simul))

lsa::cosine(t(fit_cl$fit$params$alpha_prior)) %>% pheatmap::pheatmap(display_numbers=T)



## Example #####
cls = x.simul$color_palette

simul = readRDS("nobuild/simulations/synthetic_datasets_2507/simul.N150.G1.s1.1.Rds") %>% create_basilica_obj_simul()
fit_cl = readRDS("nobuild/simulations/fits_dn.flat_clust.nonparametric.noreg.old_hier.2507/fit_clust.N150.G1.s1.1.Rds")

fit_cl2 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                              new_hier=TRUE, do_initial_fit=TRUE,
                              save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)

fit_cl3 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                              new_hier=FALSE, save_runs_seed=TRUE,
                              # do_initial_fit=TRUE,
                              save_all_fits=TRUE, nonparametric=TRUE)

fit_cl4 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                              new_hier=FALSE, save_runs_seed=TRUE,
                              do_initial_fit=TRUE,
                              save_all_fits=TRUE, nonparametric=TRUE)

fit2 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=NULL,
                           enforce_sparsity2=F, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                           lr=0.005, n_steps=2000, py=py,
                           hyperparameters=list("alpha_sigma"=0.05),
                           seed_list=c(10, 33, 92),
                           reg_weight=0., regularizer="noreg", new_hier=F,
                           save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=T)

fit3 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=NULL,
                           enforce_sparsity2=F, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                           lr=0.005, n_steps=2000, py=py,
                           hyperparameters=list("alpha_sigma"=1), seed_list=c(10, 33, 92),
                           reg_weight=0., regularizer="noreg", new_hier=F,
                           save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=T)

x.fit_cl %>% convert_sigs_names(x.simul) %>% plot_exposures(cls=cls) %>%
  patchwork::wrap_plots(x.simul %>% plot_exposures(cls=cls), ncol=1)





## Sigprofiler #####
stats2 = stats_df %>% dplyr::select(N, G, mse_counts, idd, inf_type) %>%
  # dplyr::mutate(N=paste0("N",N), G=paste0("G",G)) %>%
  dplyr::add_row(
    df_rerror %>% dplyr::rename(mse_counts=error, N=nsamples, G=groups, idd=id) %>%
      dplyr::mutate(inf_type="sigprofiler",
                    N=stringr::str_replace_all(N,"N","") %>% as.integer(),
                    G=stringr::str_replace_all(G,"G","") %>% as.integer(),
                    mse_counts = mse_counts/100)
  )

plot_mse_cosine(stats2, colname="mse_counts", runs_names=stats2$inf_type %>% unique())



df_rerror = sigprofiler_error
df_rerror$groups = factor(df_rerror$groups, levels=c('G2', 'G4', 'G6', 'G10'))
df_rerror$nsamples = factor(df_rerror$nsamples, levels=c('N50', 'N300', 'N1000', 'N5000'))
df_plot = df_rerror[!is.na(df_rerror$error),]

stats2 %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(N), y=mse_counts, color=inf_type, alpha=0) ) +
  facet_grid(~G, scales = "free_x") +
  theme_bw() +
  labs(title="Frobenius error between true and reconstructed counts (sigprofiler)") +
  ylim(0,0.2)





## example #####

fitname = list.files(fits_path, pattern="fit.")[1]
x.simul = readRDS(paste0(data_path, fitname %>% stringr::str_replace_all("fit.hier.","simul."))) %>%
  create_basilica_obj_simul()
x.fit = readRDS(paste0(fits_path, fitname)) %>% convert_sigs_names(x.simul)






generate_and_run(catalogue = COSMIC_filt,
                 comb_matrix = comb_i,
                 py = py,
                 inference_type = inference_type,
                 private_fracs = list("rare"=0.01, "common"=0.3),
                 fits_path = fits_path,
                 data_path = data_path,

                 seeds = 1,
                 mut_range = 10:8000,
                 reference_catalogue = COSMIC_filt[c("SBS1","SBS5"), ],
                 input_catalogue = NULL,
                 keep_sigs = c("SBS1", "SBS5"),
                 hyperparameters = list("alpha_sigma"=0.05),
                 lr = 0.005,
                 n_steps = 2000,

                 enforce_sparsity = T,
                 nonparametric = F,

                 reg_weight = 0.,
                 regularizer = "noreg",
                 new_hier = F,

                 initializ_seed = FALSE,
                 initializ_pars_fit = FALSE,
                 save_runs_seed = TRUE,
                 save_all_fits = TRUE,
                 #   seed_list = c(10),
                 seed_list = c(10, 33, 92),

                 do.fits = TRUE,
                 verbose = FALSE,
                 check_present = FALSE,
                 cohort = "test")



