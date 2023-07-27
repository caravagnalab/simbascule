devtools::load_all()
load_deps()

main_path = "~/Github/simbasilica/nobuild/simulations/"

save_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"

data_path = paste0(main_path, "synthetic_datasets_2507/")
fits_path = c(paste0(main_path, "fits_dn.flat_clust.parametric.noreg.old_hier.2507/"),
               paste0(main_path, "fits_dn.flat_clust.nonparametric.noreg.old_hier.2507/"))

run_id = c("noreg.old_hier.param", "noreg.old_hier.nonparam") %>% setNames(fits_path)


cutoff = 0.8; df_id = "2507"
# stats_df = get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=cutoff, fits_pattern=c("fit_clust.")) %>%
#   dplyr::mutate(run_id=run_id[fits_path],
#                 unique_id=paste(inf_type, run_id, sep=".")) %>%
#
#   dplyr::add_row(
#     get_stats_df(data_path=data_path, fits_path=fits_path, cutoff=cutoff, fits_pattern=c("fit.")) %>%
#         dplyr::mutate(run_id=run_id[fits_path])
#   ) %>%
#   dplyr::mutate(clust_type=dplyr::case_when(
#     grepl(".nonparam", run_id) ~ "non-parametric",
#     grepl(".param", run_id) ~ "parametric",
#     .default="flat")
#   )
# saveRDS(stats_df, paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))
stats_df = readRDS(paste0(save_path, "stats_df.sim", cutoff*100, ".", df_id, ".Rds"))


stats_df %>%
  dplyr::mutate(rare_ratio=n_rare_found/n_rare) %>%
  ggplot() +
  geom_jitter(aes(x=rare_freq, y=rare_ratio, color=inf_type), height=0.01, width=0.01) +
  ylim(-0.01,1+0.01) + xlim(-0.01,1+0.01) +
  facet_grid(G~clust_type) +
  theme_bw()



# pdf(paste0(save_path, "stats_report.sim", cutoff*100, ".", df_id, ".pdf"), height=6, width=10)

plot_sigs_clusters_found(stats_df, what="all", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
# plot_sigs_clusters_found(stats_df, what="common", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
plot_sigs_clusters_found(stats_df, what="rare", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))

plot_mse_cosine(stats_df, "mse_counts", scales="free_y", ylim=c(0,NA))
plot_mse_cosine(stats_df, "mse_expos", scales="free_y", ylim=c(0,NA))
plot_mse_cosine(stats_df, "mse_expos_rare", scales="free_y", ylim=c(0,NA))

plot_mse_cosine(stats_df, "cosine_sigs", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df, "cosine_expos", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df, "cosine_expos_rare", scales="free_y", ylim=c(0,1))

plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"), what="clusters", facet=T, ratio=T, scales="free_y", ylim=c(0,NA))
# plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"), what="clusters", facet=T, ratio=F, scales="free_y")

plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "nmi", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "nmi_rare", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "ari", scales="free_y", ylim=c(0,1))
plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), "ari_rare", scales="free_y", ylim=c(0,1))

# dev.off()

G = 2
N = 100
K = c(3, 4)
set.seed(1234)
pi_unn = sample(1:100, size=G)
pi = pi_unn / sum(pi_unn)

privs = c("SBS4","SBS6","SBS13")
beta = COSMIC_filt[c("SBS1","SBS5",privs),]

counts = tibble::tibble()
alpha = tibble::tibble()
alpha_prior = tibble::tibble()

for (gid in 1:G) {
  n_g = round(N * pi[gid])
  k_g = K[gid]
  n_muts_g = reticulate::r_to_py(sample(1000:5000, size=n_g))

  beta_g = beta[c("SBS1","SBS5", sample(privs, k_g-2)), ]
  alpha_prior_g = sample(1:100, size=k_g) %>% setNames(rownames(beta_g))
  privs = setdiff(privs, rownames(beta_g))

  data_g = py$generate_model(reticulate::r_to_py(alpha_prior_g),
                             beta_g, n_muts_g, N=n_g,
                             seed=as.integer(15), use_normal=T, alpha_sigma=0.4)
  # data_g$data$group = gid
  counts = dplyr::bind_rows(counts, data_g$data %>%
                              tibble::rownames_to_column(var="sample") %>%
                              dplyr::mutate(sample=paste0("G",gid,"_",sample)))
  alpha = dplyr::bind_rows(alpha, data_g$alpha %>%
                             tibble::rownames_to_column(var="sample") %>%
                             dplyr::mutate(sample=paste0("G",gid,"_",sample))) %>%
    replace(is.na(.), 0)
  alpha_prior = dplyr::bind_rows(alpha_prior, alpha_prior_g/sum(alpha_prior_g)) %>%
    replace(is.na(.), 0)
}

alpha %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=sample, y=value, fill=variable), stat="identity")

counts %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=sample, y=value, fill=variable), stat="identity")

alpha_prior = c(0.6,0.3,0.1)
beta = reticulate::r_to_py(COSMIC_filt[c("SBS1","SBS5","SBS17b"),])
n_muts = reticulate::r_to_py(sample(1000:5000,N))
data_alpha1 = py$generate_model(reticulate::r_to_py(alpha_prior), beta, n_muts, N=N, alpha_sigma=0.4,
                                seed=as.integer(15), use_normal=F)
data_alpha1$alpha %>% tibble::rownames_to_column() %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=rowname, y=value, fill=variable), stat="identity") +
  geom_hline(yintercept=1-cumsum(reticulate::py_to_r(alpha_prior)))


data_alpha2 = py$generate_model(reticulate::r_to_py(rev(alpha_prior)), beta, n_muts, N=N,
                                seed=as.integer(15), use_normal=T, alpha_sigma=0.4)

data_alpha1$alpha %>% dplyr::add_row(data_alpha2$alpha) %>%
  tibble::rownames_to_column() %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=as.integer(rowname), y=value, fill=variable), stat="identity") +
  geom_hline(yintercept=1-cumsum(alpha_prior))

counts = data_alpha1$data %>% dplyr::add_row(data_alpha2$data)

test = two_steps_inference(x=counts %>% dplyr::select(-sample), k=0:6, clusters=1:5,
                    enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                    lr=0.005, n_steps=2000, py=py,
                    hyperparameters=list("alpha_sigma"=0.05),
                    seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                    new_hier=FALSE, do_initial_fit=FALSE,
                    save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)




rnorm(1000, mean=reticulate::py_to_r(alpha_prior)[1],
        sd=data_alpha$alpha_sigma$numpy()[1]) %>% hist(breaks = 100)

rcauchy(1000, location=reticulate::py_to_r(alpha_prior)[1],
      scale=data_alpha$alpha_sigma$numpy()[1]) %>% hist(breaks = 100)



## Example ####
idd_good = stats_df %>%
  dplyr::filter(clust_type=="non-parametric", nmi_rare==max(nmi_rare, na.rm=T)) %>%
  dplyr::pull(idd)
idd_bad = stats_df %>%
  dplyr::filter(clust_type=="non-parametric", nmi_rare<0.4) %>%
  dplyr::pull(idd)

simul_good = readRDS(paste0(data_path, "simul.", idd_good[3], ".Rds")) %>% create_basilica_obj_simul()
fit_cl_good = readRDS(paste0(fits_path[2], "fit_clust.", idd_good[3], ".Rds")) %>% convert_sigs_names(simul_good)
fit_cl_good %>% plot_fit(simul_good, add_centroid = T)
fit_cl_good %>% plot_exposures(centroids = T)

simul_bad = readRDS(paste0(data_path, "simul.", idd_bad[3], ".Rds")) %>% create_basilica_obj_simul()
simul_bad$groups = get_groups_rare(simul_bad)
fit_cl_bad = readRDS(paste0(fits_path[2], "fit_clust.", idd_bad[3], ".Rds")) %>% convert_sigs_names(simul_bad)
fit_cl_bad %>% plot_fit(simul_bad, add_centroid = T)
plot_exposures()

for (gid in fit_cl_good$groups %>% unique()) {
  samples_gid = get_group(fit_cl_good, groupIDs=gid, return_idx=T)
  exposures_gid = get_exposure(fit_cl_good)[samples_gid,]
  centroids_gid = (fit_cl_good$fit$params$alpha_prior /
                     rowSums(fit_cl_good$fit$params$alpha_prior))[gid+1,]
  lapply(1:nrow(exposures_gid), function(i) exposures_gid[i,] - centroids)
  exposures_gid - centroids
}


## Example #####

simul = readRDS("nobuild/simulations/synthetic_datasets_2507/simul.N150.G1.s1.1.Rds") %>% create_basilica_obj_simul()
fit_cl = readRDS("nobuild/simulations/fits_dn.flat_clust.nonparametric.noreg.old_hier.2507/fit_clust.N150.G1.s1.1.Rds")
cls = simul$color_palette

fit_cl2 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                              new_hier=TRUE, do_initial_fit=TRUE,
                              save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)

# fit_cl3 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
#                               enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
#                               lr=0.005, n_steps=2000, py=py,
#                               hyperparameters=list("alpha_sigma"=0.05),
#                               seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
#                               new_hier=FALSE, save_runs_seed=TRUE,
#                               # do_initial_fit=TRUE,
#                               save_all_fits=TRUE, nonparametric=TRUE)

fit_cl4 = two_steps_inference(x=get_data(simul), k=fit_cl$k_list, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                              new_hier=FALSE, save_runs_seed=TRUE,
                              do_initial_fit=TRUE,
                              save_all_fits=TRUE, nonparametric=TRUE)

x.fit_cl %>% convert_sigs_names(x.simul) %>% plot_exposures(cls=cls) %>%
  patchwork::wrap_plots(x.simul %>% plot_exposures(cls=cls), ncol=1)



## Real data ####
crc = read.csv("~/GitHub/simbasilica/nobuild/processed_data/SBS_v2.03/catalogues/GEL/catalogues_Colorectal_SBS.tsv",
               sep="\t") %>% t() %>% as.data.frame()

set.seed(1234)
crc_subsample = crc[sample(1:nrow(crc), size=150), ]
fit_crc = two_steps_inference(x=crc_subsample, k=0:10, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                              new_hier=FALSE, save_runs_seed=TRUE,
                              do_initial_fit=TRUE,
                              save_all_fits=TRUE, nonparametric=TRUE)

fit_crc_cat = two_steps_inference(x=crc_subsample, k=0:10, clusters=1:5,
                                  enforce_sparsity2=T, reference_catalogue=COSMIC_filt,
                                  lr=0.005, n_steps=2000, py=py,
                                  hyperparameters=list("alpha_sigma"=0.05),
                                  seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                                  new_hier=FALSE, save_runs_seed=TRUE,
                                  do_initial_fit=TRUE,
                                  save_all_fits=TRUE, nonparametric=TRUE)

fit_crc_cat2 = two_steps_inference(x=crc_subsample, k=5:10, clusters=1:5,
                              enforce_sparsity2=T, reference_catalogue=COSMIC_filt,
                              lr=0.005, n_steps=2000, py=py,
                              hyperparameters=list("alpha_sigma"=0.05),
                              seed_list=c(10, 33, 92), reg_weight=1., regularizer="cosine",
                              new_hier=FALSE, save_runs_seed=TRUE,
                              do_initial_fit=TRUE,
                              save_all_fits=TRUE, nonparametric=TRUE)

fit_crc_cat3 = two_steps_inference(x=crc_subsample, k=5:10, clusters=1:5,
                                   enforce_sparsity2=T, reference_catalogue=COSMIC_filt,
                                   lr=0.005, n_steps=2000, py=py,
                                   hyperparameters=list("alpha_sigma"=0.05),
                                   seed_list=c(10, 33, 92), reg_weight=1., regularizer="cosine",
                                   new_hier=FALSE, save_runs_seed=TRUE,
                                   do_initial_fit=TRUE, regul_denovo = FALSE,
                                   save_all_fits=TRUE, nonparametric=TRUE)

fit_crc %>% convert_sigs_names(reference_cat = COSMIC_filt, cutoff = 0.75) %>%
  plot_fit(add_centroid = T)
get_assigned_missing(fit_crc, reference_cat = COSMIC_filt_merged, cutoff = 0.75)


## Sigprofiler #####
# stats2 = stats_df %>% dplyr::select(N, G, mse_counts, idd, inf_type) %>%
#   # dplyr::mutate(N=paste0("N",N), G=paste0("G",G)) %>%
#   dplyr::add_row(
#     df_rerror %>% dplyr::rename(mse_counts=error, N=nsamples, G=groups, idd=id) %>%
#       dplyr::mutate(inf_type="sigprofiler",
#                     N=stringr::str_replace_all(N,"N","") %>% as.integer(),
#                     G=stringr::str_replace_all(G,"G","") %>% as.integer(),
#                     mse_counts = mse_counts/100)
#   )
#
# plot_mse_cosine(stats2, colname="mse_counts", runs_names=stats2$inf_type %>% unique())
#
#
#
# df_rerror = sigprofiler_error
# df_rerror$groups = factor(df_rerror$groups, levels=c('G2', 'G4', 'G6', 'G10'))
# df_rerror$nsamples = factor(df_rerror$nsamples, levels=c('N50', 'N300', 'N1000', 'N5000'))
# df_plot = df_rerror[!is.na(df_rerror$error),]
#
# stats2 %>%
#   ggplot() +
#   geom_boxplot(aes(x=as.factor(N), y=mse_counts, color=inf_type, alpha=0) ) +
#   facet_grid(~G, scales = "free_x") +
#   theme_bw() +
#   labs(title="Frobenius error between true and reconstructed counts (sigprofiler)") +
#   ylim(0,0.2)





## example #####







generate_and_run(catalogue = COSMIC_filt,
                 comb_matrix = comb_i,
                 py = py,
                 inference_type = inference_type,
                 private_fracs = list("rare"=0.01, "common"=0.3),
                 fits_path = fits_path,
                 data_path = data_path,

                 seeds = 1,
                 mut_range = 10:8000,
                 reference_catalogue = COSMIC_filt[c("SBS1","SBS5"), ],
                 input_catalogue = NULL,
                 keep_sigs = c("SBS1", "SBS5"),
                 hyperparameters = list("alpha_sigma"=0.05),
                 lr = 0.005,
                 n_steps = 2000,

                 enforce_sparsity = T,
                 nonparametric = F,

                 reg_weight = 0.,
                 regularizer = "noreg",
                 new_hier = F,

                 initializ_seed = FALSE,
                 initializ_pars_fit = FALSE,
                 save_runs_seed = TRUE,
                 save_all_fits = TRUE,
                 #   seed_list = c(10),
                 seed_list = c(10, 33, 92),

                 do.fits = TRUE,
                 verbose = FALSE,
                 check_present = FALSE,
                 cohort = "test")



