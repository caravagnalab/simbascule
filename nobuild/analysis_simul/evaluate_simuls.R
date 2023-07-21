devtools::load_all()
load_deps()

save_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"

# data_path ="~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/"
# fits_path1 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn_1606.cosine.old_hier/",
#                "~/GitHub/simbasilica/nobuild/simulations/fits_dn_1606.noreg.old_hier/")
#
# fits_path2 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.cosine.new_hier/",
#                "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.cosine.old_hier/",
#                "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.noreg.new_hier/",
#                "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.noreg.old_hier/",
#                "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.KL.new_hier/",
#                "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.KL.old_hier/")

# run_id1 = c("cosine.nohier", "noreg.nohier") %>% setNames(fits_path1)
# run_id2 = c("cosine.new", "cosine.old", "noreg.new", "noreg.old", "KL.new", "KL.old") %>%
#   setNames(fits_path2)


data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_0507/"
fits_path1 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn.flat_clust.noreg.old_hier.0507/")
run_id1 = c("noreg.old_hier") %>% setNames(fits_path1)

fits_path2 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn.clust.nonparametric.noreg.old_hier.2007/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.clust.parametric.noreg.old_hier.2007/")
run_id2 = c("noreg.old_hier.nonparam", "noreg.old_hier.param") %>% setNames(fits_path2)

cutoff = 0.8; df_id = "2007"
stats_df = get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_clust.")) %>%
  dplyr::mutate(run_id=run_id2[fits_path],
                unique_id=paste(inf_type, run_id, sep=".")) %>%

  dplyr::add_row(
    get_stats_df(data_path=data_path, fits_path=fits_path1, cutoff=cutoff, fits_pattern=c("fit.")) %>%
        dplyr::mutate(run_id=run_id1[fits_path])
  ) %>%
  dplyr::mutate(clust_type=dplyr::case_when(
    grepl(".nonparam", run_id) ~ "non-parametric",
    grepl(".param", run_id) ~ "parametric",
    .default="flat")
  )
saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))

# cutoff = 0.8; df_id = "1807"
# stats_df = readRDS(paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))

# stats_df = get_stats_df(data_path=data_path, fits_path=fits_path1, cutoff=cutoff, fits_pattern=c("fit.")) %>%
#   dplyr::mutate(run_id=run_id1[fits_path]) %>%
#
#   # dplyr::add_row(
#   #   get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_hier.")) %>%
#   #     dplyr::mutate(run_id=run_id2[fits_path])
#   # ) %>%
#
#   dplyr::add_row(
#     get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_clust.")) %>%
#       dplyr::mutate(run_id=run_id2[fits_path])
#   ) %>%
#   dplyr::mutate(unique_id=paste(inf_type, run_id, sep=".")) %>%
#   dplyr::mutate(regularizer=dplyr::case_when(
#     grepl("cosine", run_id) ~ "cosine",
#     grepl("KL", run_id) ~ "KL",
#     grepl("noreg", run_id) ~ "noreg"
#   ), model=dplyr::case_when(
#     grepl("old", run_id) ~ "old_hier",
#     grepl("new", run_id) ~ "new_hier",
#     grepl("nohier", run_id) ~ "no_hier",
#   ))
#
#
# saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))



pdf(paste0(save_path, "stats_report.sim", cutoff*100, ".", df_id, ".pdf"), height=6, width=10)

plot_sigs_clusters_found(stats_df, what="all", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
plot_sigs_clusters_found(stats_df, what="common", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
plot_sigs_clusters_found(stats_df, what="rare", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))

plot_sigs_clusters_found(stats_df, what="all", facet=T, ratio=F, scales="free_y")
plot_sigs_clusters_found(stats_df, what="common", facet=T, ratio=F, scales="free_y")
plot_sigs_clusters_found(stats_df, what="rare", facet=T, ratio=F, scales="free_y")

plot_mse_cosine(stats_df, "mse_counts", scales="free_y", ylim=c(0,NA))
plot_mse_cosine(stats_df, "mse_expos", scales="free_y", ylim=c(0,NA))
plot_mse_cosine(stats_df, "mse_expos_rare", scales="free_y", ylim=c(0,NA))

plot_mse_cosine(stats_df, "cosine_sigs", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df, "cosine_expos", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df, "cosine_expos_rare", scales="free_y", ylim=c(0,1))

plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"), what="clusters", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"), what="clusters", facet=T, ratio=F, scales="free_y")

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


