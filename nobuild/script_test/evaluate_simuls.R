devtools::load_all("~/GitHub/basilica/")
devtools::load_all("~/GitHub/simbasilica/")
library(ggplot2)

fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_new_model_0806/"
# fits_path = "~/GitHub/simbasilica/nobuild/simulations/run_oldmodel_regul_0806/"
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/"

fits = list.files(path=fits_path,
                   pattern="fit.")

stats = lapply(fits, function(fitname) {
  print(fitname)
  compare_single_fit(fitname, fits_path, data_path, cutoff=0.8)
  } ) %>% do.call(what=rbind, args=.) %>%
  dplyr::mutate(inf_type=ifelse(is_hierarchical, "Hierarchical", "Non-hierarchical"))


mse_counts = stats %>%
  ggplot() +
  geom_jitter(aes(x=as.factor(N), y=mse_counts, color=inf_type), size=.1) +
  geom_boxplot(aes(x=as.factor(N), y=mse_counts, color=inf_type), alpha=0) +
  facet_grid(~G, scales="free_x") +
  theme_bw() + labs(title="MSE between true and reconstructed counts")


mse_expos = stats %>%
  ggplot() +
  geom_jitter(aes(x=as.factor(N), y=mse_expos, color=inf_type), size=.1) +
  geom_boxplot(aes(x=as.factor(N), y=mse_expos, color=inf_type), alpha=0) +
  facet_grid(~G, scales="free_x") +
  theme_bw() + labs(title="MSE between true and estimated exposures")


cos_sigs = stats %>%
  ggplot() +
  geom_jitter(aes(x=factor(N), y=cosine_sigs, color=inf_type), size=.1) +
  geom_boxplot(aes(x=factor(N), y=cosine_sigs, color=inf_type), alpha=0) +
  facet_grid(~G, scales="free_x") +
  theme_bw() + labs(title="Cosine between true and estimated signatures")


cos_expos = stats %>%
  ggplot() +
  geom_jitter(aes(x=factor(N), y=cosine_expos, color=inf_type), size=.1) +
  geom_boxplot(aes(x=factor(N), y=cosine_expos, color=inf_type), alpha=0) +
  facet_grid(~G, scales="free_x") +
  theme_bw() + labs(title="Cosine between true and estimated exposures")


pdf(file = "~/GitHub/simbasilica/nobuild/simulations/run_newmodel_0806.pdf",
    height = 8, width = 8)
mse_counts %>% print()
mse_expos %>% print()
cos_sigs %>% print()
cos_expos %>% print()
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
