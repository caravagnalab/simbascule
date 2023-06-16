### Datasets with diff rare fractions #####
main_path = "~/GitHub/"
data_path = "~/GitHub/simbasilica/nobuild/analysis_simul/"
fits_path = "~/GitHub/simbasilica/nobuild/analysis_simul/bug_hier_fixed/"
x.simul_old = readRDS("~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets/simul.N1000.G6.s1.Rds")
new_model = TRUE

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))
library(ggplot2)

comb = tibble(
  N_vals = c(100), n_groups_vals = c(2), samples_per_group = list(10:100), n_priv_comm = c(2), n_priv_rare = c(1)
)

rare_common = x.simul_old %>% create_basilica_obj_simul() %>% rare_common_sigs()
rare_fracs = seq(0.05, .3, length.out=3) %>% round(digits=2)
private = list("rare"=rare_common$private_rare,
               "common"=rare_common$private_common)

fits_path = "~/GitHub/simbasilica/nobuild/analysis_simul/bug_hier_fixed_noreg/"

lapply(rare_fracs, function(rare_f)
  generate_and_run(shared = rare_common$shared,
                   private = unlist(private),
                   private_fracs = list("rare"=rare_f, "common"=0.3),
                   catalogue = COSMIC_filt_merged,
                   reference_catalogue = COSMIC_filt_merged[c("SBS1","SBS5"),],
                   keep_sigs = c("SBS1","SBS5"),
                   comb_matrix = comb,
                   py = py,
                   fits_path = fits_path,
                   data_path = data_path,
                   seeds = c(1),
                   reg_weight = 0.,
                   regularizer = "cosine",
                   do.fits = TRUE,
                   new_model = TRUE,

                   seed_list = c(4,12,7),
                   initializ_seed = FALSE,
                   initializ_pars_fit = TRUE,
                   save_runs_seed = TRUE,

                   cohort = paste0("raref", rare_f))
  )

fits = list.files(path=fits_path,
                  pattern="fit.")

stats = lapply(fits, function(fitname) {
  print(fitname)
  compare_single_fit(fitname, fits_path=fits_path, data_path=data_path, cutoff=0.8)
  } ) %>% do.call(what=rbind, args=.) %>%
  dplyr::mutate(inf_type=ifelse(is_hierarchical, "Hierarchical", "Non-hierarchical")) %>%
  dplyr::mutate(rare_f=stringr::str_replace_all(idd, "N1000.G6.s1.raref",""))


colors_hier = c("darkorange", "dodgerblue4") %>%
  setNames(c("Hierarchical", "Non-hierarchical"))

stats %>%
  ggplot() +
  geom_jitter(aes(x=n_priv_rare_found, y=n_priv_rare, color=inf_type), height=0.1, width=0) +
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



## evaluate diff rates ####

fname = "N100.G2.s1"

make_plots(fname, rare_f=0.05, data_path=data_path, fits_path=fits_path)
make_plots(fname, rare_f=0.17, data_path=data_path, fits_path=fits_path)
make_plots(fname, rare_f=0.3, data_path=data_path, fits_path=fits_path)



## load fits ####
rare_f = 0.17
simul1 = readRDS(paste0("~/GitHub/simbasilica/nobuild/analysis_simul/simul.", fname, ".raref", rare_f, ".Rds")) %>%
  create_basilica_obj_simul()

fit1 = readRDS(paste0("~/GitHub/simbasilica/nobuild/analysis_simul/fit.", fname,
                      ".raref", rare_f, ".Rds")) %>%
  convert_sigs_names(x.simul=simul1, cutoff=0.8) %>%
  add_groups(simul1$groups)

hier1 = readRDS(paste0("~/GitHub/simbasilica/nobuild/analysis_simul/fit.hier.",
                       fname, ".raref", rare_f, ".Rds")) %>%
  convert_sigs_names(x.simul=simul1, cutoff=0.8)



