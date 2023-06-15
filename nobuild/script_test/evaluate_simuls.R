load_deps()

data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/"

fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_oldmodel_regul_0806/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_new_model_0806/"
fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_new_model_wholeCat_1206/"

out_name = "old_model_0806"
out_name = "new_model_0806"
out_name = "wholeCat"



## Stats df and plots ####
# stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=0.8)
# saveRDS(stats_df, paste0("~/GitHub/simbasilica/nobuild/simulations/stats.", out_name, ".Rds"))

stats_df = readRDS(paste0("~/GitHub/simbasilica/nobuild/simulations/stats.", out_name, ".Rds"))

pdf(file = paste0("~/GitHub/simbasilica/nobuild/simulations/plot.", out_name, ".pdf"),
    height = 8, width = 14)

plot_sigs_stats(stats_df, wrap=T)
plot_sigs_stats(stats_df, wrap=T, ratio=T)
plot_metrics_stats(stats_df, wrap=T)

dev.off()



## Example good ####
simul_id1 = stats %>% filter(mse_expos == min(stats$mse_expos)) %>% dplyr::pull(idd)

x.simul1 = readRDS(paste0(data_path, "simul.", simul_id1, ".Rds")) %>%
  create_basilica_obj_simul()
x.fit1 = readRDS(paste0(fits_path, "fit.", simul_id1, ".Rds"))
x.fit.hier1 = readRDS(paste0(fits_path, "fit.hier.", simul_id1, ".Rds"))

assigned.fit1 = compare_sigs_inf_gt(get_signatures(x.fit1), get_signatures(x.simul1), cutoff=0.8)
unassigned.fit1 = c(setdiff(rownames(get_signatures(x.fit1)), assigned.fit1),
                    setdiff(rownames(get_signatures(x.simul1)), names(assigned.fit1)))


x.fit1b = x.fit1$fit$secondBest %>% create_basilica_obj()
plot_similarity_reference(x.fit1, reference = x.simul1 %>% get_signatures())
plot_similarity_reference(x.fit1b, reference = x.simul1 %>% get_signatures())

plot1 = plot_fit(x.fit1)
plot2 = plot_fit(x.fit.hier1)
plot_similarity_reference(x.fit1, reference = get_signatures(x.simul1))



plot_mutations(x.simul1)
plot_exposures(x.simul1)
plot_signatures(x.simul1)

plot_signatures(x.fit1)
plot_signatures(x.fit.hier1)


## Example bad ####
simul_id = stats %>% filter(mse_expos == max(stats$mse_expos)) %>% dplyr::pull(idd)
# simul_id = "N50.G2.s4"

x.simul = readRDS(paste0(data_path, "simul.", simul_id, ".Rds")) %>%
  create_basilica_obj_simul()
x.fit = readRDS(paste0(fits_path, "fit.", simul_id, ".Rds"))
x.fit.hier = readRDS(paste0(fits_path, "fit.", simul_id, ".Rds"))


plot_fit(x.simul)
plot_mutations(x.simul)
plot_exposures(x.simul)
plot_signatures(x.simul)

plot_signatures(x.fit)
plot_signatures(x.fit.hier)





#########
x.fit.hier = readRDS("~/GitHub/simbasilica/nobuild/simulations/run_oldmodel_regul_0806/fit.hier.N1000.G6.s1.Rds")
x.fit = readRDS("~/GitHub/simbasilica/nobuild/simulations/run_oldmodel_regul_0806/fit.N1000.G6.s1.Rds")
x.simul = readRDS("~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/simul.N1000.G6.s1.Rds") %>%
  create_basilica_obj_simul()
rare_common_sigs(x.simul)

plot_fit(x.simul)
plot_exposures(x.simul, sort_by = "SBS28")

plot_exposures(x.fit)
plot_exposures(x.fit.hier)

assigned = compare_sigs_inf_gt(get_signatures(x.fit),
                               get_signatures(x.simul), cutoff=0.8)
assigned2 = compare_sigs_inf_gt(get_signatures(x.fit.hier),
                               get_signatures(x.simul), cutoff=0.8)

unassigned = c(setdiff(rownames(get_signatures(x.fit)), assigned),
               setdiff(rownames(get_signatures(x.simul)), names(assigned)))
unassigned2 = c(setdiff(rownames(get_signatures(x.fit.hier)), assigned2),
               setdiff(rownames(get_signatures(x.simul)), names(assigned2)))

compute.mse(get_exposure(x.fit), get_exposure(x.simul), assigned)
compute.mse(get_exposure(x.fit.hier), get_exposure(x.simul), assigned2)

compute.cosine(get_exposure(x.fit), get_exposure(x.simul), assigned, unassigned, what="expos")
compute.mse(get_exposure(x.fit.hier), get_exposure(x.simul), assigned2)



# min_k = max(1, nrow(x.simul$exp_fixed[[1]]) + nrow(x.simul$exp_denovo[[1]]) - 5)
# max_k = min_k + 10
# k_list = min_k:max_k
# py = reticulate::import_from_path("pybasilica", "~/GitHub/pybasilica/")
#
# x.fit2 = run_model(x = x.simul$x[[1]],
#                    k = k_list,
#                    py = py,
#                    reference_catalogue = COSMIC_filt_merged,
#                    input_catalogue = COSMIC_filt_merged[c("SBS1","SBS5"),],
#                    # keep_sigs = keep_sigs,
#                    reg_weight = 1.,
#                    CUDA = FALSE,
#                    regularizer = "cosine",
#                    filtered_cat = TRUE,
#                    verbose = TRUE,
#                    groups = x.simul$groups[[1]] - 1,
#                    new_model = FALSE,
#                    error_file = failed,
#                    idd = idd)
# close(failed)


x.fit2$fit1 %>% plot_similarity_reference(reference=get_signatures(x.simul %>% create_basilica_obj_simul()))
x.fit2$fit.hier %>% plot_similarity_reference(reference=get_signatures(x.simul %>% create_basilica_obj_simul()))
x.fit.hier %>% plot_similarity_reference(reference=get_signatures(x.simul %>% create_basilica_obj_simul()))
x.simul %>% create_basilica_obj_simul() %>% plot_mutations()





