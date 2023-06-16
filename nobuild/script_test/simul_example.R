devtools::load_all()
load_deps()

fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_new_model_0806/"
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/"

out_name = "new_model_0806"

## Get stats df ####
# stats_df = get_stats_df(data_path, fits_path, cutoff=0.8)
stats_df = readRDS(paste0("~/GitHub/simbasilica/nobuild/simulations/stats.",
                          out_name, ".Rds"))

stats_df %>%
  dplyr::mutate(rare_ratio=n_priv_rare_found/n_priv_rare,
                common_ratio=n_priv_common_found/n_priv_common) %>%
  dplyr::filter(N==300) %>%
  dplyr::select(idd, inf_type, rare_ratio, common_ratio) %>%
  dplyr::arrange(rare_ratio, common_ratio)



## Compare previous and new runs ####
idd = "N300.G4.s6"
cutoff = 0.8
save_path = "~/GitHub/simbasilica/nobuild/simulations/simul_example/"

simul = readRDS(paste0(data_path, "simul.", idd, ".Rds")) %>%
  create_basilica_obj_simul()
fit0 = readRDS(paste0(fits_path, "fit.", idd, ".Rds"))%>%
  convert_sigs_names(simul, cutoff=cutoff)
hier0 = readRDS(paste0(fits_path, "fit.hier.", idd, ".Rds")) %>%
  convert_sigs_names(simul, cutoff=cutoff)


fit1 = paste0(save_path, "fit.", idd, ".1.Rds")
hier1 = paste0(save_path, "fit.hier.", idd, ".1.Rds")
# hier1b = paste0(save_path, "fit.hier.", idd, ".1b.Rds")

fit2 = paste0(save_path, "fit.", idd, ".2.Rds")
hier2 = paste0(save_path, "fit.hier.", idd, ".2.Rds")
# hier2b = paste0(save_path, "fit.hier.", idd, ".2b.Rds")


pll0 = make_plots(fname=idd, data_path=data_path, fits_path=fits_path,
                  save_path=save_path, idd_name=".0")
pll1 = make_plots(fname=idd, data_path=data_path, fits_path=save_path,
                  save_path=save_path, idd_name=".1")
pll2 = make_plots(fname=idd, data_path=data_path, fits_path=save_path,
                  save_path=save_path, idd_name=".2")




## Run new fits ####
# seed_list = c(20,55,76,90)
# hyperparameters = list("alpha_var"=0.5, "alpha_prior_var"=1)
#
# fit1 = two_steps_inference(x=get_data(simul), k=0:12,
#                            reference_catalogue=COSMIC_filt_merged[get_fixed_signames(simul),],
#                            keep_sigs=c("SBS1","SBS5"), filtered_catalogue=TRUE,
#                            hyperparameters=NULL, initializ_seed=FALSE,
#                            seed_list=seed_list, save_runs_seed=TRUE)$tot
# fit1 = fit1 %>%
#   convert_sigs_names(simul, cutoff=cutoff) %>%
#   add_groups(simul$groups)
#
# hier1 = two_steps_inference(x=get_data(simul), k=0:12,
#                             reference_catalogue=COSMIC_filt_merged[get_fixed_signames(simul),],
#                             keep_sigs=c("SBS1","SBS5"), filtered_catalogue=TRUE,
#                             groups=simul$groups-1,
#                             hyperparameters=NULL,
#                             initializ_seed=FALSE, seed_list=seed_list,
#                             save_runs_seed=TRUE)$tot %>%
#   convert_sigs_names(simul, cutoff=cutoff)
#
# hier1b = two_steps_inference(x=get_data(simul), k=0:12,
#                             reference_catalogue=COSMIC_filt_merged[get_fixed_signames(simul),],
#                             keep_sigs=c("SBS1","SBS5"), filtered_catalogue=TRUE,
#                             groups=simul$groups-1,
#                             hyperparameters=hyperparameters,
#                             initializ_seed=FALSE, seed_list=seed_list,
#                             save_runs_seed=TRUE)$tot %>%
#   convert_sigs_names(simul, cutoff=cutoff)
#
# saveRDS(fit1, paste0(save_path, "fit.", idd, ".1.Rds"))
# saveRDS(hier1, paste0(save_path, "fit.hier.", idd, ".1.Rds"))
# saveRDS(hier1b, paste0(save_path, "fit.hier.", idd, ".1b.Rds"))




## Run with regularization #####
# fit2 = two_steps_inference(x=get_data(simul), k=0:12, reg_weight=1.,
#                            reference_catalogue=COSMIC_filt_merged[get_fixed_signames(simul),],
#                            keep_sigs=c("SBS1","SBS5"), filtered_catalogue=TRUE,
#                            hyperparameters=NULL, initializ_seed=FALSE,
#                            seed_list=seed_list, save_runs_seed=TRUE)$tot %>%
#   convert_sigs_names(simul, cutoff=cutoff) %>%
#   add_groups(simul$groups)
#
# hier2 = two_steps_inference(x=get_data(simul), k=0:12, reg_weight=1.,
#                             reference_catalogue=COSMIC_filt_merged[get_fixed_signames(simul),],
#                             keep_sigs=c("SBS1","SBS5"), filtered_catalogue=TRUE,
#                             groups=simul$groups-1,
#                             hyperparameters=NULL,
#                             initializ_seed=FALSE, seed_list=seed_list,
#                             save_runs_seed=TRUE)$tot %>%
#   convert_sigs_names(simul, cutoff=cutoff)
#
# hier2b = two_steps_inference(x=get_data(simul), k=0:12, reg_weight=1.,
#                              reference_catalogue=COSMIC_filt_merged[get_fixed_signames(simul),],
#                              keep_sigs=c("SBS1","SBS5"), filtered_catalogue=TRUE,
#                              groups=simul$groups-1,
#                              hyperparameters=hyperparameters,
#                              initializ_seed=FALSE, seed_list=seed_list,
#                              save_runs_seed=TRUE)$tot %>%
#   convert_sigs_names(simul, cutoff=cutoff)
#
# saveRDS(fit2, paste0(save_path, "fit.", idd, ".2.Rds"))
# saveRDS(hier2, paste0(save_path, "fit.hier.", idd, ".2.Rds"))








