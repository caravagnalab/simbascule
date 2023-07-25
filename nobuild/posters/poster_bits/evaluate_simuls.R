devtools::load_all()
load_deps()

main_path = "~/GitHub/simbasilica/nobuild/poster_bits/"
save_path = paste0(main_path, "results_simuls/")
data_path = paste0(main_path, "synthetic_datasets_1606/")

images_path = paste0("~/Google Drive/My Drive/work/signatures/bits/images/")


## Eval on denovo runs ####
fits_path = paste0(main_path, "run_newmodel_1606_dn/"); out_name = "dn"

# stats_dn = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)
# saveRDS(stats_dn, paste0(save_path, "stats.", out_name, ".Rds"))
stats_dn = readRDS(paste0(save_path, "stats.", out_name, ".Rds"))

stats_plots(stats_dn, save_path, out_name)




## Eval on denovo runs with no regularization ####
fits_path = paste0(main_path, "run_newmodel_1706_dn_noreg/"); out_name = "dn_noreg"

# stats_dn_nr = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)
# saveRDS(stats_dn_nr, paste0(save_path, "stats.", out_name, ".Rds"))
stats_dn_nr = readRDS(paste0(save_path, "stats.", out_name, ".Rds"))

stats_plots(stats_dn_nr, save_path, out_name)



## plots poster #####
fits_path = paste0(main_path, "run_newmodel_1606_wholeCat/"); out_name = "wholeCat"
stats_wholeCat = readRDS(paste0(save_path, "stats.", out_name, ".Rds"))
stats_nr_tmp = stats_dn_nr %>%
  dplyr::filter(G!=10, inf_type!="Hierarchical") %>%
  dplyr::mutate(inf_type="Denovo") %>%
  dplyr::add_row(
    stats_wholeCat %>%
      dplyr::filter(G!=10, inf_type!="Hierarchical") %>%
      dplyr::mutate(inf_type="Catalogue")
  )

runs_names1 = c("Denovo", "Catalogue")

stats_nr_tmp = stats_nr_tmp %>%
  dplyr::filter(inf_type=="Denovo") %>%
  dplyr::mutate(inf_type="Non-hierarchical")
sigs_found = plot_sigs_found(stats_nr_tmp, which="all", ratio=T,
                             compare_runs=F, facet_groups=F)
mse_expos = plot_mse_cosine(stats_nr_tmp, colname="mse_expos",
                            facet_groups=F, compare_runs=F) + ylab("RMSE")
cosine_sigs = plot_mse_cosine(stats_nr_tmp, colname="cosine_sigs",
                              facet_groups=F, compare_runs=F)
mse_rare = plot_mse_cosine(stats_nr_tmp, colname="mse_expos_rare",
                              facet_groups=F, compare_runs=F) + ylab("RMSE")
mse_counts = plot_mse_cosine(stats_nr_tmp, colname="mse_counts",
                             facet_groups=F, compare_runs=F) + ylab("RMSE")

patchwork::wrap_plots(sigs_found, cosine_sigs, mse_expos, mse_rare, nrow=1)
ggsave(paste0(images_path, out_name, ".simul_stats.pdf"), width = 18, height = 4)




## Eval on runs with whole catalogue ####
fits_path = paste0(main_path, "run_newmodel_1606_wholeCat/"); out_name = "wholeCat"

# stats_wholeCat = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)
# saveRDS(stats_wholeCat, paste0(save_path, "stats.", out_name, ".Rds"))
stats_wholeCat = readRDS(paste0(save_path, "stats.", out_name, ".Rds"))

stats_plots(stats_wholeCat, save_path, out_name)


