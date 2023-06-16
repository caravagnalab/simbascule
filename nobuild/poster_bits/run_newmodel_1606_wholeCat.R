args = commandArgs(trailingOnly = TRUE)
print(args)

main_path = "/u/cdslab/ebusca00/scratch_shared/basilica_pkgs/"
fits_path = paste0(main_path, "simbasilica/nobuild/poster_bits/run_new_model_wholeCat_1606/")
data_path = paste0(main_path, "simbasilica/nobuild/poster_bits/synthetic_datasets_1606/")

new_model = TRUE

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


if (args[1]=="NULL") {
  comb_i = comb[4:nrow(comb), ] 
} else {
  i = as.integer(args[1])
  comb_i = comb[i+1, ]
  print(i)
  print(comb_i)
}

shared = c("SBS1", "SBS5", "SBS17b")
private = c("SBS10b", "SBS28", "SBS56 SBS10a", "SBS90", "SBS2", "SBS13", "SBS20", "SBS22")


generate_and_run(shared = shared,
                 private = private,
                 catalogue = COSMIC_filt_merged,
                 comb_matrix = comb_i,
                 py = py,
                 private_fracs = list("rare"=0.05, "common"=0.3),
                 fits_path = fits_path,
                 data_path = data_path,
                 seeds = 1:10,
                 mut_range = 10:8000,
                 reference_catalogue = COSMIC_filt_merged,
                 input_catalogue = NULL,
                 keep_sigs = c("SBS1", "SBS5"),
                 hyperparameters = NULL,

                 reg_weight = 1.,
                 regularizer = "cosine",

                 initializ_seed = FALSE,
                 initializ_pars_fit = TRUE,
                 save_runs_seed = TRUE,
                 seed_list = c(10,27,33,92,111),

                 CUDA = TRUE,
                 do.fits = TRUE,
                 verbose = FALSE,
                 new_model = TRUE,
		 cohort = "")

