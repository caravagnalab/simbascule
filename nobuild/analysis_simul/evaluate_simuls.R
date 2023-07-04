devtools::load_all()
load_deps()

data_path ="~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/"
# fits_path = "~/GitHub/simbasilica/nobuild/poster_bits/run_new_model_wholeCat_1606/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_newmodel_1606_dn/"


out_name = "dn_2606"

stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)

plot_sigs_stats(stats_df, wrap=T)
plot_sigs_stats(stats_df, wrap=T, ratio=T)
plot_metrics_stats(stats_df, wrap=T)
# plot_mse_cosine(stats_df, colname="cosine_expos_rare")  # probably error
plot_mse_cosine(stats_df, colname="mse_expos_rare")



## example

fitname = list.files(fits_path, pattern="fit.")[1]
x.simul = readRDS(paste0(data_path, fitname %>% stringr::str_replace_all("fit.hier.","simul."))) %>%
  create_basilica_obj_simul()
x.fit = readRDS(paste0(fits_path, fitname)) %>% convert_sigs_names(x.simul)


