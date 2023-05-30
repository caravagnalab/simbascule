reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = "~/GitHub/pybasilica/")

devtools::load_all("~/GitHub/basilica")
devtools::load_all("~/GitHub/simbasilica")

source("~/GitHub/simbasilica/nobuild/script_test/helper_fns.R")
library(ggplot2)

comb = tibble(
  N_vals = c(50, 300, 300, 1000, 1000, 5000, 5000),
  n_groups_vals = c(2, 2, 4, 4, 6, 6, 10),
  samples_per_group = list(25:50, 50:300, 50:300, 100:1000, 100:1000, 200:5000, 200:5000),
  n_priv_comm = c(2, 2, 4, 4, 5, 5, 5),
  n_priv_rare = c(1, 1, 1, 1, 2, 2, 3)
)

shared = c("SBS1", "SBS5", "SBS17b")
private = c("SBS10b", "SBS28", "SBS56 SBS10a", "SBS90", "SBS2", "SBS13", "SBS20", "SBS22")

generate_synthetic_datasets(shared = shared,
                            private = private,
                            catalogue = COSMIC_filt_merged,
                            comb_matrix = comb,
                            py = py,
                            CUDA = TRUE,
                            reg_weight = 0,
                            out_path = "~/dati_elenab/signatures/simulations_2905/",
                            seeds = 1,
                            do.fits = TRUE)
