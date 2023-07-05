devtools::load_all()
load_deps()

data_path ="~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/"
save_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"

fits_path1 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn_1606.cosine.old_hier/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn_1606.noreg.old_hier/")

fits_path2 = c("~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.cosine.new_hier/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.cosine.old_hier/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.noreg.new_hier/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.noreg.old_hier/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.KL.new_hier/",
               "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.KL.old_hier/")

run_id1 = c("cosine.nohier", "noreg.nohier") %>% setNames(fits_path1)
run_id2 = c("cosine.new", "cosine.old", "noreg.new", "noreg.old", "KL.new", "KL.old") %>%
  setNames(fits_path2)

cutoff = 0.6
stats_df = get_stats_df(data_path=data_path, fits_path=fits_path1, cutoff=cutoff, fits_pattern=c("fit.")) %>%
  dplyr::mutate(run_id=run_id1[fits_path]) %>%

  dplyr::add_row(
    get_stats_df(data_path=data_path, fits_path=fits_path2, cutoff=cutoff, fits_pattern=c("fit_hier.")) %>%
      dplyr::mutate(run_id=run_id2[fits_path])
  ) %>%

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

saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".0407.Rds"))



pdf(paste0(save_path, "stats_report.sim", cutoff*100, ".pdf"), height=6, width=10)
stats_df %>%
  dplyr::mutate(ratio = n_sigs_found / n_sigs) %>%
  ggplot() +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5, color="grey") +
  geom_violin(aes(x=as.factor(N), y=ratio, fill=inf_type),
              alpha=0.5, position=position_dodge(width=0.7), color=NA) +

  geom_jitter(aes(x=as.factor(N), y=ratio, color=inf_type),
              position=position_jitterdodge(jitter.width=0.05,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("Ratio") +
  labs(title=paste0("Ratio between signatures found and true")) +
  # ylim(0-0.01, NA) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")


stats_df %>%
  dplyr::mutate(ratio = n_common_found / n_common) %>%
  ggplot() +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5, color="grey") +
  geom_violin(aes(x=as.factor(N), y=ratio, fill=inf_type),
              alpha=0.5, position=position_dodge(width=0.7), color=NA) +

  geom_jitter(aes(x=as.factor(N), y=ratio, color=inf_type),
              position=position_jitterdodge(jitter.width=0.05,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("Ratio") +
  labs(title=paste0("Ratio between common signatures found and true")) +
  # ylim(0-0.01, NA) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")


stats_df %>%
  dplyr::mutate(ratio = n_rare_found / n_rare) %>%
  ggplot() +
  geom_hline(yintercept=1, linetype="dashed", linewidth=0.5, color="grey") +
  geom_violin(aes(x=as.factor(N), y=ratio, fill=inf_type),
              alpha=0.5, position=position_dodge(width=0.7), color=NA) +

  geom_jitter(aes(x=as.factor(N), y=ratio, color=inf_type),
              position=position_jitterdodge(jitter.width=0.05,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("Ratio") +
  labs(title=paste0("Ratio between rare signatures found and true")) +
  # ylim(0-0.1, NA) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")


stats_df %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(N), y=mse_counts, color=inf_type),
               alpha=0.5, position=position_dodge(width=0.7), outlier.colour=NA) +

  geom_jitter(aes(x=as.factor(N), y=mse_counts, color=inf_type),
              position=position_jitterdodge(jitter.width=0.01,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("MSE") +
  labs(title=paste0("MSE counts")) +
  # ylim(-0.01, 1) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")

stats_df %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(N), y=mse_expos, color=inf_type),
               alpha=0.5, position=position_dodge(width=0.7), outlier.colour=NA) +

  geom_jitter(aes(x=as.factor(N), y=mse_expos, color=inf_type),
              position=position_jitterdodge(jitter.width=0.01,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("MSE") +
  labs(title=paste0("MSE exposures")) +
  # ylim(-0.01, 1) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")


stats_df %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(N), y=cosine_sigs, color=inf_type),
               alpha=0.5, position=position_dodge(width=0.7), outlier.colour=NA) +

  geom_jitter(aes(x=as.factor(N), y=cosine_sigs, color=inf_type),
              position=position_jitterdodge(jitter.width=0.01,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("Cosine") +
  labs(title=paste0("Cosine signatures")) +
  # ylim(-0.01, 1) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")


stats_df %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(N), y=cosine_expos, color=inf_type),
               alpha=0.5, position=position_dodge(width=0.7), outlier.colour=NA) +

  geom_jitter(aes(x=as.factor(N), y=cosine_expos, color=inf_type),
              position=position_jitterdodge(jitter.width=0.01,
                                            jitter.height=0,
                                            dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("Cosine") +
  labs(title=paste0("Cosine exposures")) +
  # ylim(-0.01, 1) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")



stats_df %>%
  dplyr::filter(inf_type == "fit_clust") %>%
  dplyr::mutate(ratio = n_groups_found / G) %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(N), y=ratio),
               alpha=0.5, position=position_dodge(width=0.7), outlier.colour=NA) +

  geom_jitter(aes(x=as.factor(N), y=ratio), height = 0,
              # position=position_jitterdodge(jitter.width=0.01,
              #                               jitter.height=0.01,
              #                               dodge.width=0.7),
              size=0.8) +
  # scale_color_manual(values=cols, name="") +
  # scale_fill_manual(values=cols, name="") +
  theme_bw() + xlab("# samples") + ylab("Ratio") +
  labs(title=paste0("N groups found")) +
  # ylim(-0.01, 1) +
  ggh4x::facet_nested(G ~ regularizer + model, scales="free_y")


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


