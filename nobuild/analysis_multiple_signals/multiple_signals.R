devtools::load_all("~/GitHub/simbasilica/")
load_deps()


## Generayte and run ####

private = list("SBS"=c("SBS17b", "SBS4", "SBS7c", "SBS13", "SBS20", "SBS22"),
               "DBS"=c("DBS1", "DBS2", "DBS7", "DBS11", "DBS4", "DBS10"))
shared = list("SBS"=c("SBS1","SBS5"),
              "DBS"=c("DBS3","DBS5"))

N_list = c(150,300,1000)
G_list = c(2,4,6)
n_gs = expand.grid(N_list, G_list)

easypar_pars = lapply(1:nrow(n_gs), function(i) {
  list(N=n_gs[i,"Var1"], G=n_gs[i,"Var2"],
       private=private,
       shared=shared)
})


easypar_fn = function(N, G, private, shared, run_fits=FALSE) {
  fname = paste0("simul_fit_", N, ".", G, "_macpro.Rds")
  fname_fpath = paste0("~/Dropbox/shared/2022. Basilica/simulations/matched_signals/",
                       fname)

  simul_ng = NULL
  if (file.exists(fname_fpath))
    simul_ng = readRDS(fname_fpath)$dataset

  devtools::load_all("~/GitHub/basilica/")
  devtools::load_all("~/GitHub/simbasilica/")
  reticulate::use_condaenv("basilica-env")
  library(ggplot2)
  py = reticulate::import_from_path("pybasilica", "~/GitHub/pybasilica/")

  if (is.null(simul_ng)) {
    seed_list = list("SBS"=N+G,
                     "DBS"=N+G*2)

    simul_ng = generate_simulation_dataset_matched(N=N, G=G,
                                                   private=private,
                                                   shared=shared,
                                                   alpha_range=c(0.25,0.3),
                                                   alpha_sigma=0.2,
                                                   py=py,
                                                   seed=seed_list) %>%
      create_basilica_obj_simul()
  }

  cat("Simulated dataset generated!")

  x_ng = NULL

  if (run_fits) {
    counts_ng = get_input(simul_ng, matrix=TRUE)
    max_K = max(sapply(get_signames(simul_ng), length)) + 1
    x_ng = fit(counts=counts_ng, k_list=0:max_K, cluster=NULL, n_steps=3000,
               reference_cat=list("SBS"=COSMIC_filt[shared$SBS,],
                                  "DBS"=COSMIC_dbs[shared$DBS,]),
               keep_sigs=unlist(shared),
               hyperparameters=list("scale_factor_centroid"=5000,
                                    "scale_factor_alpha"=5000,
                                    "tau"=0),
               seed_list=c(10,33,4), filter_dn=TRUE, store_fits=TRUE)
  }

  saveRDS(list("dataset"=simul_ng, "fit"=x_ng),
          paste0("~/Dropbox/shared/2022. Basilica/simulations/matched_signals/", fname))
  return(list("dataset"=simul_ng, "fit"=x_ng))
}

# easypar_output = easypar::run(FUN=easypar_fn,
#                               PARAMS=easypar_pars,
#                               parallel=FALSE,
#                               silent=FALSE,
#                               filter_errors=FALSE,
#                               outfile="~/Dropbox/shared/2022. Basilica/simulations/matched_signals/easypar.log")

lapply(easypar_pars[2:length(easypar_pars)], function(i) {
  print(i)
  easypar_fn(N=i$N, G=i$G, private=i$private, shared=i$shared, run_fits=TRUE)
})



## Example ####

path = "~/Dropbox/shared/2022. Basilica/simulations/matched_signals/"
simul_fit_i = readRDS(paste0(path, "simul_fit_150.2_macpro.Rds"))
simul_i = simul_fit_i$dataset
fit_i = simul_fit_i$fit

counts_ng = list("SBS"=get_input(fit_i)$SBS %>% long_to_wide(what="counts"))
x_ng = fit(counts=counts_ng, k_list=0:3, cluster=NULL, n_steps=3000,
           reference_cat=list("SBS"=COSMIC_filt[shared$SBS,]),
           keep_sigs=unlist(shared),
           hyperparameters=list("scale_factor_centroid"=5000,
                                "scale_factor_alpha"=5000,
                                "tau"=0, "pi_conc0"=1e-3),
           seed_list=c(10), filter_dn=FALSE, store_fits=TRUE, store_parameters=FALSE)

x_ng2 = fit(counts=counts_ng, k_list=0:3, cluster=NULL, n_steps=100,
           reference_cat=list("SBS"=COSMIC_filt[shared$SBS,]),
           keep_sigs=unlist(shared),
           hyperparameters=list("scale_factor_centroid"=5000,
                                "scale_factor_alpha"=5000,
                                "tau"=0, "pi_conc0"=0.6),
           seed_list=c(10), filter_dn=FALSE, store_fits=TRUE, store_parameters=TRUE)

# muts = plot_data(x_ng, types = "SBS", reconstructed=T)
muts_true = plot_data(x_ng, types = "SBS", reconstructed=F)
sigs = plot_signatures(x_ng)
weights = plot_beta_weights(x_ng)
beta_star = plot_beta_star(x_ng)
alpha_star = plot_alpha_star(x_ng)
patchwork::wrap_plots(sigs, muts_true, design="AABB
                                               CCBB")

patchwork::wrap_plots(alpha_star, beta_star, weights, design="AABB
                                                              CCBB")
plot_exposures(x_ng)
plot_gradient_norms(x_ng)

x_ng3 = fit(counts=counts_ng, k_list=3, cluster=NULL, n_steps=10,
            reference_cat=list("SBS"=COSMIC_filt[shared$SBS,]),
            keep_sigs=unlist(shared),
            hyperparameters=list("scale_factor_centroid"=5000,
                                 "scale_factor_alpha"=5000,
                                 "tau"=0, "pi_conc0"=1e-3),
            seed_list=c(10), filter_dn=TRUE, store_fits=TRUE, store_parameters=TRUE)

alts = get_alternatives(x_ng, what="nmf", types="SBS")[["SBS"]]$fits$`k_denovo:2`$`seed:10`[[1]]
tmp = x_ng
tmp$nmf$SBS$exposure = alts$exposure
tmp$nmf$SBS$beta_denovo = alts$beta_denovo
tmp$nmf$SBS$pyro = alts

tmp %>% plot_signatures()

p1 = plot_beta_weights(x_ng)
p2 = plot_beta_weights(x_ng2)
p3 = plot_beta_weights(x_ng3)

patchwork::wrap_plots(p1, p2, p3, guides="collect")


## Data analysis #####
path = "~/Dropbox/shared/2022. Basilica/simulations/matched_signals/"
simul_fit_i = readRDS(paste0(path, "simul_fit_150.2_macpro.Rds"))
simul_i = simul_fit_i$dataset
fit_i = simul_fit_i$fit

get_params(fit_i, what="nmf", type="SBS")[[1]]$beta_w


plots = lapply(list.files(path, pattern="^simul_fit", full.names=T), function(fname) {
  simul_fit_i = readRDS(fname)
  simul_i = simul_fit_i$dataset
  fit_i = simul_fit_i$fit
  colpalette = gen_palette(simul_i)
  colpalette_fit = gen_palette(fit_i)
  design = "AACC
            BBCC
            BBCC"
  design2 = "AACC
             BBCC
             BBCC"
  simul_plots = patchwork::wrap_plots(plot_data(simul_i, reconstructed=F),
                                      plot_exposures(simul_i, cls=colpalette) +
                                        theme(legend.position="bottom"),
                                      plot_signatures(simul_i, cls=colpalette),
                                      design=design)
  fit_plots = patchwork::wrap_plots(plot_data(fit_i, reconstructed=F),
                                    plot_exposures(fit_i, cls=colpalette_fit) +
                                      theme(legend.position="bottom"),
                                    plot_signatures(fit_i, cls=colpalette_fit),
                                    plot_beta_weights(fit_i),
                                    design=design2)
  return(
    list("simul"=simul_plots,
         "fit"=fit_plots)
  )
}) %>% setNames(list.files(path, pattern="^simul_fit", full.names=F))

pdf(paste0(path, "datasets.pdf"), height=10, width=12)
lapply(names(plots), function(i)
  plots[[i]] + labs(title=stringr::str_replace_all(i, pattern="_macpro.Rds", "")))
dev.off()


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




## Plotting functions ####

plot_beta_weights = function(x) {
  params = get_QC(x, what="nmf", types="SBS")[["SBS"]]$train_params
  tmp = params %>% dplyr::filter(paramname=="beta_w") %>% dplyr::select(-paramname) %>%
    dplyr::mutate(iteration=iteration) # %>% dplyr::filter(iteration %in% c(10, 100, 300, 500))
  tmp %>% ggplot(aes(x=iteration, y=rowname, colour=columnname, size=value)) +
    geom_point(position=position_dodge(width=0.5))
}





