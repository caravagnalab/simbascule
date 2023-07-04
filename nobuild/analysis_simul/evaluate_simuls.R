devtools::load_all()
load_deps()

data_path ="~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/"
# fits_path = "~/GitHub/simbasilica/nobuild/poster_bits/run_new_model_wholeCat_1606/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/fits_dn.hier_clust.cosine.new_hier/"


out_name = "dn_2606"

stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)

plot_sigs_stats(stats_df, wrap=T)
plot_sigs_stats(stats_df, wrap=T, ratio=T)
plot_metrics_stats(stats_df, wrap=T)
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


