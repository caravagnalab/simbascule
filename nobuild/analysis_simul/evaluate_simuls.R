devtools::load_all()
load_deps()

main_path = "~/Dropbox/shared/2022. Basilica/simulations/"
save_path = paste0(main_path, "stats_dataframes/")
data_path = paste0(main_path, "synthetic_datasets_3107/")

cutoff = 0.8; min_expos=0.; df_id = "sf5000vlearning.1809"

fits_path = c(
  paste0(main_path, "fits_dn.clust.sf5000.dmm.1809/"),
  paste0(main_path, "fits_dn.clust.sf_learning.dmm.1809/")
  )
run_id = c("sf5000.DMM", "sfLearning.DMM") %>% setNames(fits_path)
stats_df = get_stats_df(data_path=data_path, fits_path=fits_path,
                        cutoff=cutoff, fits_pattern=c("fit.", "fit_clust."),
                        run_id=run_id,
                        min_exposure=min_expos,
                        save_plots=FALSE, check_plots=FALSE) %>%

  dplyr::mutate(clust_type=dplyr::case_when(
    grepl(".nonparam", fits_path) ~ "non-parametric",
    grepl(".param", fits_path) ~ "parametric",
    .default="flat"),
    clust_type=ifelse(grepl("clust",fits_pattern), clust_type, "non-clustering")
  )
saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))
fname = paste0(cutoff*100, ".", df_id)
report_stats(stats_df=stats_df, fname=paste(fname,"LC",sep="."),
             save_path=save_path, fill="run_id", suffix_name="LC")
report_stats(stats_df=stats_df, fname=paste(fname,"noLC",sep="."),
             save_path=save_path, fill="run_id", suffix_name="noLC")
report_stats(stats_df=stats_df, fname=paste(fname,"LCnomerge",sep="."),
             save_path=save_path, fill="run_id", suffix_name="LCnomerge")

# stats_df = readRDS(paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))



## Plots ####
stats_df_spars = stats_df %>%
  dplyr::filter(run_id == "sf5000.DMM") %>% dplyr::mutate(run_id="DMM")

stats_df_spars = stats_df

suffix_name = "LCnomerge"
figure = make_figure(stats_df_spars, suffix_name=suffix_name)

fig_id = paste(suffix_name,df_id,sep=".")
ggsave(plot=figure, filename=paste0(save_path,"fig_simulations.",fig_id,".pdf"), height=8, width=12)

