devtools::load_all("~/GitHub/simbasilica/")
load_deps()

private = list("SBS"=c("SBS17b", "SBS4", "SBS7c", "SBS13"),
               "DBS"=c("DBS1", "DBS2", "DBS7", "DBS11"))
shared = list("SBS"=c("SBS1","SBS5"),
              "DBS"=c("DBS3","DBS5"))

N_list = c(150, 300, 1000)
G_list = c(2,4)
simul_objects = lapply(G_list, function(gid) {
  lapply(N_list, function(nid) {
    seed_list = list("SBS"=sample(1:100, size=1),
                     "DBS"=sample(1:100, size=1))
    simul_ng = generate_simulation_dataset_matched(N=nid, G=gid,
                                                    private=private,
                                                    shared=shared,
                                                    py=py,
                                                    seed=seed_list) %>%
      create_basilica_obj_simul()

    counts_ng = get_input(simul_obj, matrix=TRUE)
    max_K = sapply(get_signames(simul_obj), length) %>% max
    x_ng = fit(counts=counts, k_list=0:max_K, cluster=NULL, n_steps=3000,
               reference_cat=list("SBS"=COSMIC_filt[shared$SBS,],
                                  "DBS"=COSMIC_dbs[shared$DBS,]),
               keep_sigs=unlist(shared),
               hyperparameters=list("scale_factor_centroid"=5000,
                                    "scale_factor_alpha"=5000, "tau"=0),
               seed_list=c(10,33,4), filter_dn=TRUE, store_fits=TRUE)
    return(list("dataset"=simul_ng, "fit"=x_ng))
  }) %>% setNames(N_list)
}) %>% setNames(G_list)




simul_obj = generate_simulation_dataset_matched(N=150, G=2, private=private,
                                                shared=shared,
                                                py=py) %>%
  create_basilica_obj_simul()

counts = get_input(simul_obj, matrix=TRUE)


x = fit(counts=counts, k_list=0:2, cluster=NULL, n_steps=1000,
        reference_cat=list("SBS"=COSMIC_filt[c("SBS1","SBS5"),],
                           "DBS"=COSMIC_dbs[c("DBS3","DBS5"),]),
        keep_sigs=c("SBS1","SBS5","DBS3","DBS5"),
        hyperparameters=list("scale_factor_centroid"=5000,
                             "scale_factor_alpha"=5000, "tau"=0),
        seed_list=c(10,33,4), filter_dn=TRUE, store_fits=TRUE)
x$description = "Using Dirichlet as beta prior, centroid computed from fixed cumulative, adding penalty in pyro.factor, try learning whole omega"
saveRDS(x, "~/GitHub/simbasilica/nobuild/analysis_multiple_signals/ex_nmf3.Rds")

plot_scores(x)
plot_signatures(x)
# alt_sols = get_alternatives(x, what="nmf", types="SBS")
get_params(x, what="nmf", type="SBS")[[1]]$beta_w

beta_star = get_params(x, what="nmf", type="SBS")[[1]]$beta_star$numpy() %>%
  as.data.frame()
colnames(beta_star) = colnames(COSMIC_filt)
p_beta_star = beta_star %>% wide_to_long(what="beta") %>%
  reformat_contexts(what="SBS") %>%
  dplyr::mutate(type="SBS") %>%
  plot_signatures_aux()

alpha_star = get_params(x, what="nmf", type="SBS")[[1]]$alpha_star %>%
  as.data.frame()
p_alpha_star = alpha_star %>% wide_to_long(what="exposures") %>%
  dplyr::mutate(type="SBS") %>%
  plot_exposures_aux()

alpha_star_DBS = get_params(x, what="nmf", type="DBS")[[1]]$alpha_star %>%
  as.data.frame()
p_alpha_star_DBS = alpha_star_DBS %>% wide_to_long(what="exposures") %>%
  dplyr::mutate(type="DBS") %>%
  plot_exposures_aux()

beta_weights = get_params(x, what="nmf", type="SBS")[[1]]$beta_w$numpy()
colnames(beta_weights) = c(get_fixed_signames(x, types="SBS")$SBS, "DN")
rownames(beta_weights) = get_denovo_signames(x, types="SBS")$SBS

patchwork::wrap_plots(plot_signatures(simul_obj, types="SBS"), p_beta_star, ncol=1)
patchwork::wrap_plots(plot_signatures(x, types="SBS"), p_beta_star, ncol=1)

patchwork::wrap_plots(plot_exposures(simul_obj, types="SBS"), p_alpha_star, ncol=1)
patchwork::wrap_plots(plot_exposures(simul_obj, types="DBS"), p_alpha_star_DBS, ncol=1)

patchwork::wrap_plots(plot_exposures(simul_obj, types="SBS"), plot_exposures(x, types="SBS"), ncol=1)
heatmap(beta_weights, cluster_rows=F, cluster_cols=F)
plot_signatures(x)
# saveRDS(x, "~/GitHub/simbasilica/nobuild/analysis_multiple_signals/ex.Gamma.Rds")

patchwork::wrap_plots(plot_signatures(simul_obj, types="DBS"),
                      plot_signatures(x, types="DBS"), ncol=1)



# test = fit(counts=counts, k_list=4, cluster=NULL, n_steps=3000,
#            reference_cat=list("SBS"=COSMIC_filt[c("SBS1","SBS5"),], "DBS"=COSMIC_dbs["DBS4",]),
#            keep_sigs=c("SBS1","SBS5","DBS4"),
#            hyperparameters=list("scale_factor_centroid"=5000,
#                                 "scale_factor_alpha"=5000, "tau"=0),
#            seed_list=c(10,33,4), filter_dn=TRUE, store_fits=TRUE)

# get_alternative_run(x, params=list("seed"=10))






## real data ####
data_path = "~/Dropbox/shared/2022. Basilica/real_data/processed_data/"
tissues = c("Colorectal")
N = 500
counts_sbs = readRDS(paste0(data_path, "catalogues/counts_sbs.Rds")) %>%
  dplyr::filter(organ %in% tissues)
counts_dbs = readRDS(paste0(data_path, "catalogues/counts_dbs.Rds")) %>%
  dplyr::filter(organ %in% tissues)

common_samples = intersect(rownames(counts_sbs), rownames(counts_dbs))

set.seed(13)
snames = sample(common_samples, N, replace=F)

counts_sbs_n = counts_sbs[snames,] %>% dplyr::select(-organ, -cohort)
counts_dbs_n = counts_dbs[snames,] %>% dplyr::select(-organ, -cohort)

counts = list("SBS"=counts_sbs_n, "DBS"=counts_dbs_n)

x_crc_Unif = fit(counts=counts, k_list=5:10, cluster=6, n_steps=3000,
            reference_cat=list("SBS"=COSMIC_filt[c("SBS1","SBS5"),],
                               "DBS"=COSMIC_dbs[c("DBS2","DBS4"),]),
            keep_sigs=c("SBS1","SBS5"),
            hyperparameters=list("scale_factor_centroid"=5000,
                                 "scale_factor_alpha"=5000, "tau"=0),
            seed_list=c(10,33,4,57), filter_dn=TRUE, store_fits=TRUE)

saveRDS(x_crc_Unif, "~/GitHub/simbasilica/nobuild/analysis_multiple_signals/crc_sbs_dbs.Unif.Rds")


x_crc_Gamma = fit(counts=counts, k_list=5:10, cluster=6, n_steps=3000,
                 reference_cat=list("SBS"=COSMIC_filt[c("SBS1","SBS5"),],
                                    "DBS"=COSMIC_dbs[c("DBS2","DBS4"),]),
                 keep_sigs=c("SBS1","SBS5"),
                 hyperparameters=list("scale_factor_centroid"=5000,
                                      "scale_factor_alpha"=5000, "tau"=0),
                 seed_list=c(10,33,4,57), filter_dn=TRUE, store_fits=TRUE)

saveRDS(x_crc_Gamma, "~/GitHub/simbasilica/nobuild/analysis_multiple_signals/crc_sbs_dbs.Gamma.Rds")







