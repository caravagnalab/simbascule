devtools::load_all()
load_deps()

data_path = "~/GitHub/simbasilica/nobuild/poster_bits/synthetic_datasets_1606/"
# fits_path = "~/GitHub/simbasilica/nobuild/poster_bits/run_new_model_wholeCat_1606/"
fits_path = "~/GitHub/simbasilica/nobuild/poster_bits/run_newmodel_1606_dn/"


out_name = "newmodel_dn"

stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)

plot_sigs_stats(stats_df, wrap=T)
plot_sigs_stats(stats_df, wrap=T, ratio=T)
plot_metrics_stats(stats_df, wrap=T)
# plot_mse_cosine(stats_df, colname="cosine_expos_rare")  # probably error
plot_mse_cosine(stats_df, colname="mse_expos_rare")


