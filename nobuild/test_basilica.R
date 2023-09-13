devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica/")
# load_deps(base_path="~/GitHub/")

## Simulated ####
input.simul = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N500_G3_simul.Rds")
input.simul %>% create_basilica_obj_simul()  # creates a basilica obj with simulated data
obj_simul = input.simul %>% create_basilica_obj_simul()

x.simul = fit(x=input.simul$counts[[1]],
              k=0:2, # n of denovo signatures
              clusters=6,  # n of clusters
              n_steps=3000, lr=0.005,
              enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
              dirichlet_prior=TRUE,
              reference_catalogue=COSMIC_filt[lc,],
              hyperparameters=list("scale_factor_alpha"=100, "scale_factor_centroid"=1000),  # change default values to hyperparameters
              nonparametric=TRUE, py=py, store_parameters=TRUE)

x.simul.100 = fit(x=input.simul$counts[[1]],
              k=2:6, # n of denovo signatures
              clusters=6,  # n of clusters
              n_steps=3000, lr=0.005,
              enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
              dirichlet_prior=TRUE,
              reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
              hyperparameters=list("scale_factor_alpha"=100, "scale_factor_centroid"=100),  # change default values to hyperparameters
              nonparametric=TRUE, py=py, store_parameters=TRUE)

saveRDS(list("sf_centr100"=x.simul.100,
             "sf_centr1000"=x.simul),
        "~/Dropbox/shared/2022. Basilica/real_data/results/fit.N500_G3_simul.onestep.Rds")

saveRDS(x.simul,
        "~/Dropbox/shared/2022. Basilica/real_data/results/fit.N500_G3_simul.dmm.lc.Rds")

lc = filter_signatures_QP(get_signatures(x.simul), COSMIC_filt) %>% unlist() %>% unique()

x.simul %>% convert_sigs_names(input.simul %>% create_basilica_obj_simul()) %>%
  get_obj_initial_params() %>% plot_exposures() %>%
  patchwork::wrap_plots(x.simul %>% convert_sigs_names(input.simul %>% create_basilica_obj_simul()) %>%
                          plot_exposures(sample_name = T))



## Real data ####
input.real = readRDS("~/Dropbox/shared/2022. Basilica/datasets/input.N500_CRC.Rds")

x.real = fit(x=input.real$counts[[1]],
             k=8:10, # n of denovo signatures
             clusters=6,  # n of clusters
             n_steps=3000, lr=0.005,
             enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
             dirichlet_prior=TRUE,
             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
             hyperparameters=list("scale_factor_alpha"=100, "scale_factor_centroid"=1000),  # change default values to hyperparameters
             nonparametric=TRUE, py=py, store_parameters=TRUE)

x.real.100 = fit(x=input.real$counts[[1]],
             k=8:10, # n of denovo signatures
             clusters=6,  # n of clusters
             n_steps=3000, lr=0.005,
             enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
             dirichlet_prior=TRUE,
             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
             hyperparameters=list("scale_factor_alpha"=100, "scale_factor_centroid"=100),  # change default values to hyperparameters
             nonparametric=TRUE, py=py, store_parameters=TRUE)

lc = filter_signatures_QP(get_signatures(x.real), COSMIC_filt, filt_pi=0.2) %>% unlist() %>% unique()

saveRDS(list("sf_centr100"=x.real.100,
             "sf_centr1000"=x.real),
        "~/Dropbox/shared/2022. Basilica/real_data/results/fit.N500_CRC.onestep.Rds")


## plots ####
fit100 = obj_simul
fit1000 = x.simul
idd = "N500_G3_simul.dmm.lc"

fit1 = fit1000 %>% convert_sigs_names(reference_cat=COSMIC_filt, cutoff=.6)
fit2 = fit100 %>% convert_sigs_names(reference_cat=COSMIC_filt, cutoff=.6)

name1 = "One step - scale centr 1000"
name2 = "Real"
ref_cls = COSMIC_color_palette()[c(get_signames(fit1),get_signames(fit2))] %>%
  purrr::discard(is.na)
dn_cls = c(get_color_palette(fit1),
           get_color_palette(fit2))[grep("D",c(get_signames(fit1),
                                               get_signames(fit2)) %>% unique())]
cls = c(ref_cls, dn_cls)
pp = make_plots_compare(fit1=fit1, fit2=fit2,
                        name1=name1, name2=name1,
                        min_exposure=.05, cls=cls)

centr_iters1 = lapply(seq(1,length(fit1$fit$train_params$iteration %>% unique()),by=10), function(iter) {
  fit1$fit$train_params %>%
    dplyr::filter(paramname=="centroid", iteration==iter) %>%
    dplyr::rename(sample=rowname, Signature=columnname, alpha=value) %>%
    dplyr::mutate(sample=factor(sample, levels=sort(unique(sample)))) %>%
    plot_exposures_aux(cls=fit1000$color_palette, titlee=paste0("Iter ", (iter-1)*50), sample_name=T)
}) %>% patchwork::wrap_plots(guides="collect")

centr_iters2 = lapply(seq(1,length(fit2$fit$train_params$iteration %>% unique()),by=10), function(iter) {
  fit2$fit$train_params %>%
    dplyr::filter(paramname=="centroid", iteration==iter) %>%
    dplyr::rename(sample=rowname, Signature=columnname, alpha=value) %>%
    dplyr::mutate(sample=factor(sample, levels=sort(unique(sample)))) %>%
    plot_exposures_aux(cls=fit100$color_palette, titlee=paste0("Iter ", (iter-1)*50), sample_name=T)
}) %>% patchwork::wrap_plots(guides="collect")


centr1 = lapply(unique(fit1$groups), function(gid) {
  idxs = get_group(fit1, groupIDs=gid, return_idx=TRUE)
  if (length(idxs) == 0) next
  plot_exposures(fit1, sampleIDs=idxs, cls=cls) + labs(title="")
})
centr1[["centr"]] = plot_exposures(fit1, centroids=T, cls=cls) + labs(title="")

centr2 = lapply(unique(fit2$groups), function(gid) {
  idxs = get_group(fit2, groupIDs=gid, return_idx=TRUE)
  if (length(idxs) == 0) next
  plot_exposures(fit2, sampleIDs=idxs, cls=cls) + labs(title="")
})
centr2[["centr"]] = plot_exposures(fit2, centroids=T, cls=cls) + labs(title="")


pdf(paste0("~/Dropbox/shared/2022. Basilica/real_data/results/plots.", idd, ".pdf"), height=12, width=16)
plot_fit(fit1, fit2, cls=cls, name1=name1, name2=name2) %>% print()
pp$expos_centr %>% print()
patchwork::wrap_plots(centr1, guides="collect") %>% print()
patchwork::wrap_plots(centr2, guides="collect") %>% print()
centr_iters1 %>% print()
# centr_iters2 %>% print()
patchwork::wrap_plots(pp$umap,
                      plot_gradient_norms(fit1),
                      # plot_gradient_norms(fit2),
                      ncol=2) %>% print()
plot_posterior_probs(fit1)
# plot_posterior_probs(fit2)
dev.off()





## functions
get_obj_initial_params(x.simul)  # initial conditions
get_group(x.simul, groupIDs=c("1"), return_idx=T)  # sample names belonging to grp 1

convert_sigs_names(x.simul, reference_cat=COSMIC_filt)  # assigns denovo to catalogue

plot_exposures(x.simul)
plot_signatures(x.simul)
plot_mutations(x.simul)
plot_posterior_probs(x.simul)  # heatmap with posterior probs
plot_scores(x.simul)
plot_gradient_norms(x.simul)







