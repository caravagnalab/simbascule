main_path = "/u/cdslab/ebusca00/scratch_shared/basilica_pkgs/"
out_path = paste0(main_path, "simbasilica/nobuild/simulations/run_oldmodel_regul_0806/")
data_path = paste0(main_path, "simbasilica/nobuild/simulations/synthetic_datasets/")
new_model = FALSE

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

comb = tibble(
  N_vals = c(50, 300, 300, 1000, 1000, 5000, 5000),
  n_groups_vals = c(2, 2, 4, 4, 6, 6, 10),
  samples_per_group = list(25:50, 50:300, 50:300, 100:1000, 100:1000, 200:5000, 200:5000),
  n_priv_comm = c(2, 2, 4, 4, 5, 5, 5),
  n_priv_rare = c(1, 1, 1, 1, 2, 2, 3)
)

shared = c("SBS1", "SBS5", "SBS17b")
private = c("SBS10b", "SBS28", "SBS56 SBS10a", "SBS90", "SBS2", "SBS13", "SBS20", "SBS22")

generate_and_run(shared = shared,
                 private = private,
                 catalogue = COSMIC_filt_merged,
                 comb_matrix = comb,
                 py = py,
                 fits_path = out_path,
                 data_path = data_path,
                 seeds = 1:8,
                 mut_range = 10:8000,
                 reference_catalogue = COSMIC_filt_merged,
                 input_catalogue = COSMIC_filt_merged[c("SBS1","SBS5"), ],
                 keep_sigs = c("SBS1", "SBS5"),
                 reg_weight = 1.,
                 regularizer = "cosine",
                 CUDA = TRUE,
                 do.fits = TRUE,
                 verbose = FALSE,
                 new_model = new_model)

