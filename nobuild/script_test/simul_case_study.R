## Evaluate example #####

fit.hier = readRDS("~/GitHub/simbasilica/nobuild/simulations/run_new_model_0806/fit.hier.N1000.G6.s2.Rds")
simul = readRDS("~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/simul.N1000.G6.s2.Rds") %>%
  create_basilica_obj_simul()
fit = readRDS("~/GitHub/simbasilica/nobuild/simulations/run_new_model_0806/fit.N1000.G6.s2.Rds") %>%
  add_groups(simul$groups)

rare_common = rare_common_sigs(simul)

assign_fit = get_assigned_missing(fit, simul)
assign_hier = get_assigned_missing(fit.hier, simul)

get_groups_with_sigs(simul, c(rare_common$private_rare, rare_common$private_common)) %>% tidyr::unnest(data)

plot_exposures(fit)



### Datasets with diff rare fractions #####
main_path = "~/GitHub/"
simulation_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"
x.simul_old = readRDS("~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/simul.N1000.G6.s1.Rds")
new_model = TRUE

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

comb = tibble(
  N_vals = c(100), n_groups_vals = c(2), samples_per_group = list(10:100), n_priv_comm = c(2), n_priv_rare = c(1)
)

# rare_common = x.simul_old %>% create_basilica_obj_simul() %>% rare_common_sigs()
# rare_fracs = seq(0.05, .3, length.out=3) %>% round(digits=2)
# # rare_fracs = c(0.01)
# private = list("rare"=rare_common$private_rare, "common"=rare_common$private_common)
#
# lapply(rare_fracs, function(rare_f)
#   generate_and_run(shared = rare_common$shared,
#                    private = unlist(private),
#                    private_fracs = list("rare"=rare_f, "common"=0.3),
#                    catalogue = COSMIC_filt_merged,
#                    reference_catalogue = COSMIC_filt_merged[c("SBS1","SBS5"),],
#                    keep_sigs = c("SBS1","SBS5"),
#                    comb_matrix = comb,
#                    py = py,
#                    fits_path = simulation_path,
#                    data_path = simulation_path,
#                    seeds = c(1),
#                    reg_weight = 1.,
#                    regularizer = "cosine",
#                    do.fits = TRUE,
#                    new_model = TRUE,
#                    cohort = paste0("raref", rare_f))
#   )


# simul1 = readRDS("~/GitHub/simbasilica/nobuild/analysis_simul/simul.N1000.G6.s1.raref0.05.Rds") %>%
#   create_basilica_obj_simul()
# fit1 = readRDS("~/GitHub/simbasilica/nobuild/analysis_simul/fit.N1000.G6.s1.raref0.05.Rds")
# fit1.hier = readRDS("~/GitHub/simbasilica/nobuild/analysis_simul/fit.hier.N1000.G6.s1.raref0.05.Rds")

fits = list.files(path=simulation_path,
                  pattern="fit.")

stats = lapply(fits, function(fitname) {
  print(fitname)
  compare_single_fit(fitname, simulation_path, simulation_path, cutoff=0.8)
  } ) %>% do.call(what=rbind, args=.) %>%
  dplyr::mutate(inf_type=ifelse(is_hierarchical, "Hierarchical", "Non-hierarchical")) %>%
  dplyr::mutate(rare_f=stringr::str_replace_all(idd, "N1000.G6.s1.raref",""))


colors_hier = c("darkorange", "dodgerblue4") %>%
  setNames(c("Hierarchical", "Non-hierarchical"))

stats %>%
  ggplot() +
  geom_jitter(aes(x=n_priv_rare_found, y=n_priv_rare, color=inf_type), height=0.1, width=0.1) +
  geom_abline() +
  ggh4x::facet_nested(~G+as.factor(N)+as.factor(rare_f)) +
  scale_color_manual(values=colors_hier) + ylim(0, NA) + xlim(0, NA) +
  theme_bw() + labs(title="N found (x) vs N (y)")

stats %>%
  ggplot() +
  geom_jitter(aes(x=n_priv_common_found, y=n_priv_common, color=inf_type), height=0.1, width=0.1) +
  geom_abline() +
  ggh4x::facet_nested(~G+as.factor(N)+as.factor(rare_f)) +
  scale_color_manual(values=colors_hier) + ylim(0, NA) + xlim(0, NA) +
  theme_bw() + labs(title="N found (x) vs N (y)")

plot_sigs_found(stats, which="rare")
plot_sigs_found(stats, which="common")
plot_sigs_found(stats, which="all", ratio=F)



## evaluate diff rates ####

rare_f = 0.05
simul1 = readRDS(paste0("~/GitHub/simbasilica/nobuild/analysis_simul/simul.N100.G2.s1.raref", rare_f, ".Rds")) %>%
  create_basilica_obj_simul()
fit1 = readRDS(paste0("~/GitHub/simbasilica/nobuild/analysis_simul/fit.N100.G2.s1.raref", rare_f, ".Rds")) %>%
  convert_sigs_names(x.simul=simul1, cutoff=0.8) %>% add_groups(simul1$groups)
hier1 = readRDS(paste0("~/GitHub/simbasilica/nobuild/analysis_simul/fit.hier.N100.G2.s1.raref", rare_f, ".Rds")) %>%
  convert_sigs_names(x.simul=simul1, cutoff=0.8)

assigned_missing_fit1 = get_assigned_missing(x.fit=fit1, x.simul=simul1, cutoff=0.8)
assigned_missing_hier1 = get_assigned_missing(x.fit=hier1, x.simul=simul1, cutoff=0.8)
all_sigs = c(get_signames(simul1), get_signames(fit1), get_signames(hier1)) %>% unique()

rare_common = rare_common_sigs(simul1)

colpal = gen_palette(length(all_sigs)) %>% setNames(all_sigs)

pl_simul = plot_fit(simul1, cls=colpal) +
  patchwork::plot_annotation(title=paste0("SIMUL: Rare sigs ", rare_common$private_rare,
                                          " frequency ", rare_f,
                                          " - Common sigs ", rare_common$private_common))

pl_fit = plot_fit(fit1, cls=colpal, reconstructed = T) +
  patchwork::plot_annotation(title=paste0("FIT: Rare sigs ", rare_common$private_rare,
                                          " frequency ", rare_f,
                                          " - Common sigs ", rare_common$private_common),
                             subtitle=paste0("Missing sigs ", assigned_missing_fit1$missing_fn,
                                             " - Added sigs ", assigned_missing_fit1$added_fp))

pl_hier = plot_fit(hier1, cls=colpal, reconstructed = T) +
  patchwork::plot_annotation(title=paste0("HIER FIT: Rare sigs ", rare_common$private_rare,
                                          " frequency ", rare_f,
                                          " - Common sigs ", rare_common$private_common),
                             subtitle=paste0("Missing sigs ", assigned_missing_hier1$missing_fn,
                                             " - Added sigs ", assigned_missing_hier1$added_fp))

pl_simil_fit1 = plot_similarity_reference(fit1, reference=get_signatures(simul1)) +
  patchwork::plot_annotation(title=paste0("FIT: Rare sigs ", rare_common$private_rare,
                                          " frequency ", rare_f,
                                          " - Common sigs ", rare_common$private_common))
pl_simil_hier1 = plot_similarity_reference(hier1, reference=get_signatures(simul1)) +
  patchwork::plot_annotation(title=paste0("HIER FIT: Rare sigs ", rare_common$private_rare,
                                          " frequency ", rare_f,
                                          " - Common sigs ", rare_common$private_common))

pdf(paste0("GitHub/simbasilica/nobuild/analysis_simul/plots_N100.G2.s1.", rare_f, ".pdf"), height = 10, width = 14)
print(pl_simul)
print(pl_fit)
print(pl_simil_fit1)
print(pl_hier)
print(pl_simil_hier1)
dev.off()

