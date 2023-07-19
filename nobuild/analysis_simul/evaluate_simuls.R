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
# fits_path1 = fits_path2 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn.flat_clust.noreg.old_hier.0507/",
#                             "~/GitHub/simbasilica/nobuild/simulations/fits_dn.flat_clust.noreg.new_hier.0507/")
# run_id1 = run_id2 = c("noreg.old_hier", "noreg.new_hier") %>% setNames(fits_path1)

fits_path2 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn.clust.noreg.old_hier.no_sparsity/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.clust.noreg.old_hier.sparsity/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.clust.noreg.new_hier.no_sparsity/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.clust.noreg.new_hier.sparsity/")
run_id2 = c("noreg.old_hier.no_sparsity", "noreg.old_hier.sparsity",
            "noreg.new_hier.no_sparsity", "noreg.new_hier.sparsity") %>% setNames(fits_path2)

stats_df = get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_clust.")) %>%
  dplyr::mutate(run_id=run_id2[fits_path])

cutoff = 0.8; df_id = "1807"
# stats_df = readRDS(paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))

stats_df = get_stats_df(data_path=data_path, fits_path=fits_path1, cutoff=cutoff, fits_pattern=c("fit.")) %>%
  dplyr::mutate(run_id=run_id1[fits_path]) %>%

  # dplyr::add_row(
  #   get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_hier.")) %>%
  #     dplyr::mutate(run_id=run_id2[fits_path])
  # ) %>%

  dplyr::add_row(
    get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_clust.")) %>%
      dplyr::mutate(run_id=run_id2[fits_path])
  ) %>%
  dplyr::mutate(unique_id=paste(inf_type, run_id, sep=".")) %>%
  dplyr::mutate(regularizer=dplyr::case_when(
    grepl("cosine", run_id) ~ "cosine",
    grepl("KL", run_id) ~ "KL",
    grepl("noreg", run_id) ~ "noreg"
  ), model=dplyr::case_when(
    grepl("old", run_id) ~ "old_hier",
    grepl("new", run_id) ~ "new_hier",
    grepl("nohier", run_id) ~ "no_hier",
  ))
#
#
# saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))



pdf(paste0(save_path, "stats_report.sim", cutoff*100, ".", df_id, ".pdf"), height=6, width=10)
plot_sigs_found(stats_df, which="all", facet=T, ratio=T, scales="free_y")

plot_sigs_found(stats_df, which="common", facet=T, ratio=T, scales="free_y")

plot_sigs_found(stats_df, which="rare", facet=T, ratio=T, scales="free_y")

plot_mse_cosine(stats_df, "mse_counts", scales="free_y")

plot_mse_cosine(stats_df, "mse_expos", scales="free_y")

plot_mse_cosine(stats_df, "mse_expos_rare", scales="free_y")

plot_mse_cosine(stats_df, "cosine_sigs", scales="free_y")

plot_mse_cosine(stats_df, "cosine_expos", scales="free_y")

plot_mse_cosine(stats_df, "cosine_expos_rare", scales="free_y")

plot_mse_cosine(stats_df, "nmi", scales="free_y")

plot_mse_cosine(stats_df, "nmi_rare", scales="free_y")

plot_mse_cosine(stats_df, "ari", scales="free_y")

plot_mse_cosine(stats_df, "ari_rare", scales="free_y")

dev.off()


# stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8, fits_pattern=fits_pattern) %>%
#   dplyr::mutate(run_id=run_id[fits_path]) %>%
#   dplyr::mutate(unique_id=paste(inf_type, run_id, sep="."))
#
# stats_df2 = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.6, fits_pattern=fits_pattern) %>%
#   dplyr::mutate(run_id=run_id[fits_path]) %>%
#   dplyr::mutate(unique_id=paste(inf_type, run_id, sep="."))

plot_sigs_stats(stats_df, wrap=T)

stats_df %>%
  dplyr::filter(unique_id != "fit_hier.cosine") %>%
  dplyr::filter(run_id %in% c("cosine","cosine_new")) %>%
  plot_sigs_stats(wrap=T, ratio=T, facet_groups=T)

stats_df %>%
  # dplyr::filter(unique_id != "fit_hier.cosine") %>%
  # dplyr::filter(run_id %in% c("cosine","cosine_new")) %>%
  plot_metrics_stats(wrap=T)

# plot_mse_cosine(stats_df, colname="cosine_expos_rare")  # probably error
plot_mse_cosine(stats_df, colname="mse_expos_rare")
plot_mse_cosine(stats_df, colname="mse_counts")


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


