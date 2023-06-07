devtools::load_all("~/GitHub/basilica/")
devtools::load_all("~/GitHub/simbasilica/")
library(ggplot2)

fits_path = "~/GitHub/simbasilica/nobuild/simulations/fits_new_model_3105/"
data_path = "~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/"

fits = list.files(path=fits_path,
                   pattern="fit.")

stats = lapply(fits, function(fitname)
  compare_single_fit(fitname, fits_path, data_path)
  ) %>% do.call(what=rbind, args=.)


stats %>%
  ggplot() +
  geom_jitter(aes(x=is_hierarchical, y=mse_counts), size=.1) +
  geom_boxplot(aes(x=is_hierarchical, y=mse_counts), alpha=0) +
  facet_grid(N~G) +
  theme_bw() + labs(title="MSE between true and reconstructed counts")


stats %>%
  ggplot() +
  geom_jitter(aes(x=is_hierarchical, y=mse_expos), size=.1) +
  geom_boxplot(aes(x=is_hierarchical, y=mse_expos), alpha=0) +
  facet_grid(N~G) +
  theme_bw() + labs(title="MSE between true and estimated exposures")


stats %>%
  ggplot() +
  geom_jitter(aes(x=is_hierarchical, y=cosine_sigs), size=.1) +
  geom_boxplot(aes(x=is_hierarchical, y=cosine_sigs), alpha=0) +
  facet_grid(N~G) +
  theme_bw() + labs(title="Cosine between true and estimated signatures")

stats %>%
  ggplot() +
  geom_jitter(aes(x=factor(G), y=cosine_sigs), size=.1) +
  geom_boxplot(aes(x=factor(G), y=cosine_sigs), alpha=0) +
  facet_grid(N~is_hierarchical) +
  theme_bw() + labs(title="Cosine between true and estimated signatures")


stats %>%
  ggplot() +
  geom_jitter(aes(x=is_hierarchical, y=cosine_expos), size=.1) +
  geom_boxplot(aes(x=is_hierarchical, y=cosine_expos), alpha=0) +
  facet_grid(N~G) +
  theme_bw() + labs(title="Cosine between true and estimated exposures")




## Example ####
simul_id = "N50.G2.s4"

x.simul = readRDS(paste0(data_path, "simul.", simul_id, ".Rds")) %>%
  create_basilica_obj_simul()
x.fit = readRDS(paste0(fits_path, "fit.", simul_id, ".Rds"))
x.fit.hier = readRDS(paste0(fits_path, "fit.", simul_id, ".Rds"))

x.fit2 = two_steps_inference(x=x.simul$input$counts,
                             k=0:10, py=py, reg_weight=0.,
                             reference_catalogue=COSMIC_filt_merged, )


plot_fit(x.simul)
plot_mutations(x.simul)
plot_exposures(x.simul)
plot_signatures(x.simul)

plot_signatures(x.fit)
plot_signatures(x.fit.hier)
