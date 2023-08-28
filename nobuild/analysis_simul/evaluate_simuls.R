devtools::load_all()
load_deps()

main_path = "~/Dropbox/shared/2022. Basilica/simulations/simuls_lc/"
save_path = paste0(main_path, "../stats_dataframes/")
data_path = paste0(main_path, "../synthetic_datasets_3107/")

fits_path = c(paste0(main_path, "fits_dn.clust.nonparametric.sparsity.noreg.old_hier.0208/"),
              paste0(main_path, "fits_dn.clust.nonparametric.nonsparsity.noreg.old_hier.0208/"),
              paste0(main_path, "fits_dn.flat.nonparametric.sparsity.noreg.old_hier.0208/"),
              paste0(main_path, "fits_dn.flat.nonparametric.nonsparsity.noreg.old_hier.0208/"))
run_id = c("nparam.spars", "nparam.nspars", "flat.spars", "flat.nspars") %>% setNames(fits_path)

cutoff = 0.8; min_expos=0.; df_id = "lc.2408"
stats_df = get_stats_df(data_path=data_path, fits_path=fits_path,
                        cutoff=cutoff, fits_pattern=c("fit.", "fit_clust."),
                        run_id=run_id,
                        min_exposure=min_expos, save_plots=FALSE) %>%

  dplyr::mutate(clust_type=dplyr::case_when(
    grepl(".nonparam", fits_path) ~ "non-parametric",
    grepl(".param", fits_path) ~ "parametric",
    .default="flat")
  )
# saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))
# fname = paste0(cutoff*100, ".", df_id)
# report_stats(stats_df=stats_df, fname=fname, save_path=save_path, fill="run_id")

stats_df = readRDS(paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))


## Example ####
# idd_good = stats_df %>%
#   dplyr::filter(clust_type=="non-parametric", nmi_rare==max(nmi_rare, na.rm=T)) %>%
#   dplyr::pull(idd)
# idd_bad = stats_df %>%
#   dplyr::filter(clust_type=="non-parametric", nmi_rare<0.4) %>%
#   dplyr::pull(idd)
#
# simul_good = readRDS(paste0(data_path, "simul.", idd_good[3], ".Rds")) %>% create_basilica_obj_simul()
# fit_cl_good = readRDS(paste0(fits_path[2], "fit_clust.", idd_good[3], ".Rds")) %>% convert_sigs_names(simul_good)
# fit_cl_good %>% plot_fit(simul_good, add_centroid = T)
# fit_cl_good %>% plot_exposures(centroids = T)
#
# simul_bad = readRDS(paste0(data_path, "simul.", idd_bad[3], ".Rds")) %>% create_basilica_obj_simul()
# simul_bad$groups = get_groups_rare(simul_bad)
# fit_cl_bad = readRDS(paste0(fits_path[2], "fit_clust.", idd_bad[3], ".Rds")) %>% convert_sigs_names(simul_bad)
# fit_cl_bad %>% plot_fit(simul_bad, add_centroid = T)
# plot_exposures()
#
# for (gid in fit_cl_good$groups %>% unique()) {
#   samples_gid = get_group(fit_cl_good, groupIDs=gid, return_idx=T)
#   exposures_gid = get_exposure(fit_cl_good)[samples_gid,]
#   centroids_gid = (fit_cl_good$fit$params$alpha_prior /
#                      rowSums(fit_cl_good$fit$params$alpha_prior))[gid+1,]
#   lapply(1:nrow(exposures_gid), function(i) exposures_gid[i,] - centroids)
#   exposures_gid - centroids
# }


## Example #####

# simul = readRDS("nobuild/simulations/synthetic_datasets_2507/simul.N150.G1.s1.1.Rds") %>% create_basilica_obj_simul()
# fit_cl = readRDS("nobuild/simulations/fits_dn.flat_clust.nonparametric.noreg.old_hier.2507/fit_clust.N150.G1.s1.1.Rds")
# cls = simul$color_palette
#
# fit_cl2 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
#                               enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
#                               lr=0.005, n_steps=2000, py=py,
#                               hyperparameters=list("alpha_sigma"=0.05),
#                               seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
#                               new_hier=TRUE, do_initial_fit=TRUE,
#                               save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)
#
#
# ## Real data ####
# crc = read.csv("~/GitHub/simbasilica/nobuild/processed_data/SBS_v2.03/catalogues/GEL/catalogues_Colorectal_SBS.tsv",
#                sep="\t") %>% t() %>% as.data.frame()
#
# set.seed(1234)
# crc_subsample = crc[sample(1:nrow(crc), size=150), ]
# fit_crc = two_steps_inference(x=crc_subsample, k=0:10, clusters=1:5,
#                               enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
#                               lr=0.005, n_steps=2000, py=py,
#                               hyperparameters=list("alpha_sigma"=0.05),
#                               seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
#                               new_hier=FALSE, save_runs_seed=TRUE,
#                               do_initial_fit=TRUE,
#                               save_all_fits=TRUE, nonparametric=TRUE)
#
# fit_crc_cat = two_steps_inference(x=crc_subsample, k=0:10, clusters=1:5,
#                                   enforce_sparsity2=T, reference_catalogue=COSMIC_filt,
#                                   lr=0.005, n_steps=2000, py=py,
#                                   hyperparameters=list("alpha_sigma"=0.05),
#                                   seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
#                                   new_hier=FALSE, save_runs_seed=TRUE,
#                                   do_initial_fit=TRUE,
#                                   save_all_fits=TRUE, nonparametric=TRUE)
#
#
# fit_crc %>% convert_sigs_names(reference_cat = COSMIC_filt, cutoff = 0.75) %>%
#   plot_fit(add_centroid = T)
# get_assigned_missing(fit_crc, reference_cat = COSMIC_filt_merged, cutoff = 0.75)


## Sigprofiler #####
# stats2 = stats_df %>% dplyr::select(N, G, mse_counts, idd, inf_type) %>%
#   # dplyr::mutate(N=paste0("N",N), G=paste0("G",G)) %>%
#   dplyr::add_row(
#     df_rerror %>% dplyr::rename(mse_counts=error, N=nsamples, G=groups, idd=id) %>%
#       dplyr::mutate(inf_type="sigprofiler",
#                     N=stringr::str_replace_all(N,"N","") %>% as.integer(),
#                     G=stringr::str_replace_all(G,"G","") %>% as.integer(),
#                     mse_counts = mse_counts/100)
#   )
#
# plot_mse_cosine(stats2, colname="mse_counts", runs_names=stats2$inf_type %>% unique())
#
#
#
# df_rerror = sigprofiler_error
# df_rerror$groups = factor(df_rerror$groups, levels=c('G2', 'G4', 'G6', 'G10'))
# df_rerror$nsamples = factor(df_rerror$nsamples, levels=c('N50', 'N300', 'N1000', 'N5000'))
# df_plot = df_rerror[!is.na(df_rerror$error),]
#
# stats2 %>%
#   ggplot() +
#   geom_boxplot(aes(x=as.factor(N), y=mse_counts, color=inf_type, alpha=0) ) +
#   facet_grid(~G, scales = "free_x") +
#   theme_bw() +
#   labs(title="Frobenius error between true and reconstructed counts (sigprofiler)") +
#   ylim(0,0.2)
