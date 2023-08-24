devtools::load_all()
load_deps()

main_path = "~/Dropbox/shared_folders/2022. Basilica/simulations/"
save_path = paste0(main_path, "stats_dataframes/")


data_path = paste0(main_path, "synthetic_datasets_3107/")
# fits_path = c(paste0(main_path, "fits_dn.clust.nonparametric.nonsparsity.noreg.old_hier.3107/"),
#                paste0(main_path, "fits_dn.clust.nonparametric.sparsity.noreg.old_hier.3107/"))
# run_id = c("noreg.nonparam.nonsparsity", "noreg.nonparam.sparsity") %>% setNames(fits_path)

fits_path = c(paste0(main_path, "fits_dn.clust.nonparametric.sparsity.noreg.old_hier.0208/"),
              paste0(main_path, "fits_dn.clust.nonparametric.nonsparsity.noreg.old_hier.0208/"))
run_id = c("nparam.spars", "nparam.nspars") %>% setNames(fits_path)

# fits_path = paste0(main_path, "fits_dn.clust.nonparametric.nonsparsity.noreg.old_hier.cauchy_0208/")
# run_id = c("noreg.nonparam.nonsparsity.cauchy") %>% setNames(fits_path)

cutoff = 0.8; min_expos=0.; df_id = "0208"
# stats_df = get_stats_df(data_path=data_path, fits_path=fits_path,
#                         cutoff=cutoff, fits_pattern=c("fit_clust."),
#                         run_id=run_id,
#                         min_exposure=min_expos, save_plots=FALSE) %>%
#   # dplyr::mutate(run_id=run_id[fits_path],
#   #               unique_id=paste(fits_pattern, run_id, sep=".")) %>%
#
#   dplyr::mutate(clust_type=dplyr::case_when(
#     grepl(".nonparam", fits_path) ~ "non-parametric",
#     grepl(".param", fits_path) ~ "parametric",
#     .default="flat")
#   )
# saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))
stats_df = readRDS(paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))

# fname = paste0(cutoff*100, ".", df_id)
# report_stats(stats_df=stats_df, fname=fname, save_path=save_path, fill="run_id")




x.simul = readRDS(paste0(data_path, "simul.N150.G3.s4.1.Rds")) %>%
  create_basilica_obj_simul()

x = two_steps_inference(x=get_data(x.simul), k=2:6,
                        reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5","SBS90","SBS10b"),],
                        keep_sigs=c("SBS1","SBS5"), n_steps = 2000, nonparametric = T,
                        clusters=6, py = py, reg_weight = 0., hyperparameters=list("alpha_sigma"=0.1))

x2 = fix_assignments(xbis)
pp = make_plots_compare(xbis %>% convert_sigs_names(x.simul), x.simul)
pp2 = make_plots_compare(x2 %>% convert_sigs_names(x.simul), x.simul)

pdf("./tmp.pdf", height = 12, width = 12)
print(pp)
dev.off()

pdf("./tmp2.pdf", height = 12, width = 12)
print(pp2)
dev.off()


x.fit = readRDS(paste0(fits_path[1], "fit_clust.N150.G3.s4.1.Rds")) %>% convert_sigs_names(x.simul)
fit_new = x.fit %>% fix_assignments(cutoff=0.8, max_iters=20)
pl = make_plots_compare(x.fit, fit_new, "fit", "fit fixed")
pl$expos_centr


sigs = get_signatures(x.fit)
substitutions = get_contexts(x.fit) %>% dplyr::pull(subs) %>% unique()

similar = data.frame()
for (subs in substitutions) {
  sigs_s = sigs[, grepl(subs, colnames(sigs))]
  # sigs_s = sigs_s / rowSums(sigs_s)

  cosine = lsa::cosine(t(sigs_s))
  cosine[is.nan(cosine)] = 0  # lower.tri(cosine, diag=T) |

  similar = similar %>% dplyr::bind_rows(
    cosine %>% as.data.frame() %>% tibble::rownames_to_column(var="sigs1") %>%
      reshape2::melt(variable.name="sigs2", value.name="cosine") %>%
      dplyr::filter(sigs1 %in% get_fixed_signames(x.fit),
                    sigs2 %in% get_dn_signames(x.fit),
                    cosine > cutoff) %>%
      dplyr::mutate(subs=subs)
  )
}






new_fit = fix_assignments(x.fit)

plots = make_plots_compare(x.fit, new_fit, "fit init", "fit fixed")

cls = gen_palette(n=unique(c(get_signames(x.simul), get_signames(x.fit))) %>% length()) %>%
  setNames(unique(c(get_signames(x.simul), get_signames(x.fit))))
new_fit %>% convert_sigs_names(x.simul) %>% plot_exposures(add_centroid=TRUE, cls=cls) %>%
  patchwork::wrap_plots(x.fit %>% convert_sigs_names(x.simul) %>%
                          plot_exposures(add_centroid=TRUE, cls=cls),
                        x.simul %>% plot_exposures(add_centroid=TRUE, cls=cls),
                        ncol=1, guides="collect")


x.fit = get_simul_fit(stats_df, return_fit=T, condition="grepl('N150.G1.s1',idd)")$x.fit
x.simul = get_simul_fit(stats_df, return_fit=T, condition="grepl('N150.G1.s1',idd)")$x.simul

x.fit %>% plot_exposures()




nspars.g1 = get_simul_fit(stats_df, condition="idd==spars.g1$idd & grepl('.nonsparsity', run_id)",
                       return_fit=T)

spars.g1$x.fit %>% convert_sigs_names(spars.g1$x.simul) %>%
  get_centroids(normalize=T) %>% tibble::rownames_to_column() %>%
  reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=rowname, y=value, fill=variable), stat="identity") +
  scale_fill_manual(values=gen_palette(spars.g1$x.fit$color_palette %>% length()))

nspars.g1$x.fit %>% convert_sigs_names(spars.g1$x.simul) %>%
  get_centroids(normalize=T) %>% tibble::rownames_to_column() %>%
  reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=rowname, y=value, fill=variable), stat="identity") +
  scale_fill_manual(values=gen_palette(nspars.g1$x.fit$color_palette %>% length()))


simul = spars.g1$x.simul
xnew = recompute_centroids(spars.g1$x.fit) %>% convert_sigs_names(simul)
xnew %>% plot_exposures(centroids = T)
xnew %>% plot_exposures(add_centroid = T)

xnew %>% merge_clusters() %>% plot_exposures(add_centroid = T)



x = nspars.g1$x.fit %>% convert_sigs_names(spars.g1$x.simul)
simul = nspars.g1$x.simul

x_bis = two_steps_inference(x=get_data(simul), k=x$k_list, enforce_sparsity2=F,
                    clusters=get_centroids(x) %>% nrow(), n_steps=2000,
                    nonparametric=TRUE, reg_weight=0., do_initial_fit=TRUE,
                    save_all_fits=TRUE, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),])

x_bis.normal = two_steps_inference(x=get_data(simul), k=x$k_list, enforce_sparsity2=F,
                            clusters=get_centroids(x) %>% nrow(), n_steps=2000,
                            hyperparameters = list("alpha_sigma"=0.1),
                            nonparametric=TRUE, reg_weight=0., do_initial_fit=TRUE,
                            save_all_fits=TRUE, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),])

x_bis.cat = two_steps_inference(x=get_data(simul), k=x$k_list, enforce_sparsity2=F,
                            clusters=get_centroids(x) %>% nrow(), n_steps=2000,
                            nonparametric=TRUE, reg_weight=0., do_initial_fit=TRUE,
                            save_all_fits=TRUE, reference_catalogue=COSMIC_filt)

x_bis.spars = two_steps_inference(x=get_data(simul), k=x$k_list, enforce_sparsity2=T,
                            clusters=get_centroids(x) %>% nrow(), n_steps=2000,
                            nonparametric=TRUE, reg_weight=0., do_initial_fit=TRUE,
                            save_all_fits=TRUE, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),])


x_bis %>% convert_sigs_names(x.simul=simul) %>%
  plot_exposures(add_centroid=TRUE, cls=simul$color_palette) %>%
  patchwork::wrap_plots(simul %>% plot_exposures(add_centroid = T),
                        guides="collect", ncol=1)

x_bis %>% convert_sigs_names(x.simul=simul) %>% plot_exposures() %>%
  patchwork::wrap_plots(
    x_bis %>% convert_sigs_names(x.simul=simul) %>% plot_exposures(centroids=TRUE),
    widths=c(3,1), guides="collect"
  )


x_bis.normal %>% convert_sigs_names(x.simul=simul) %>%
  plot_exposures(add_centroid=TRUE, cls=simul$color_palette) %>%
  patchwork::wrap_plots(simul %>% plot_exposures(add_centroid = T),
                        guides="collect", ncol=1)


x_bis.cat %>% convert_sigs_names(x.simul=simul) %>%
  plot_exposures(add_centroid=TRUE, cls=simul$color_palette) %>%
  patchwork::wrap_plots(simul %>% plot_exposures(add_centroid = T),
                        guides="collect", ncol=1)


x_bis.spars %>% convert_sigs_names(x.simul=simul) %>%
  plot_exposures(add_centroid=TRUE, cls=simul$color_palette) %>%
  patchwork::wrap_plots(simul %>% plot_exposures(add_centroid = T),
                        guides="collect", ncol=1)


x_bis.spars %>% convert_sigs_names(x.simul=simul) %>% plot_exposures() %>%
  patchwork::wrap_plots(
    x_bis.spars %>% convert_sigs_names(x.simul=simul) %>% plot_exposures(centroids=TRUE),
    widths=c(3,1), guides="collect"
  )

# alpha_prior = get_centroids(x, normalize = T)[3,"SBS7c"]
#
# alpha_sigma = 0.1
# q_01 = alpha_prior - alpha_sigma # * alpha_prior
# q_99 = alpha_prior + alpha_sigma # * alpha_prior
# alpha_sigma_corr = (q_99 - q_01) / (2 * qnorm(p=0.99, mean=alpha_prior, sd=1))
#
# (q_99 - q_01) / (2 * qnorm(p=0.99, mean=alpha_prior, sd=1))
# (q_99 - q_01) / (2 * quantile(rcauchy(1000, location=alpha_prior, scale=1), 0.95))
#
# rcauchy(40, location=alpha_prior, scale=alpha_sigma_corr) %>% hist(breaks=100)
# rnorm(40, mean=alpha_prior, sd=alpha_sigma_corr) %>% hist(breaks=100)


spars.g6 = get_simul_fit(stats_df, return_fit=T,
                      condition="N==500 & G==3 & !grepl('.nonsparsity', run_id)")
nspars.g6 = get_simul_fit(stats_df, condition="idd==spars.g6$idd & grepl('.nonsparsity', run_id)",
                       return_fit=T)
spar

p1 = spars.g6$filtered$plot_centroids[[1]]
p2 = spars.g6$filtered$plot_expos[[1]] & theme(legend.position="none")
patchwork::wrap_plots(p1, p2, widths=c(1,3))
spars.g6$filtered$plot_expos %>% patchwork::wrap_plots(spars.g6$filtered$plot_centroids[[1]], widths=c(3, 1))
spars.g6$filtered$plot_centroids
spars.g6$filtered$plot_sigs

nspars.g6$filtered$plot_expos



spars_umap = umap::umap(get_exposure(spars$x.fit))




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
