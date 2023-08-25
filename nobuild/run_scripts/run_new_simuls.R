args = commandArgs(trailingOnly = TRUE)
cat(paste("Arguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
nonparametric = args[2] == "nonparametric"
enforce_sparsity = args[3] == "sparsity"
inference_type = args[4]
run_id = args[5]

new_hier = FALSE
regularizer = "noreg"
reg_weight = 0.

cat(paste("i =", i, "new_hier =", new_hier, "regularizer =", regularizer, "enforce_sparsity=", enforce_sparsity,
          "inference_type =", inference_type, "nonparametric =", nonparametric, "\n"))

main_path = "~/GitHub/"
data_path = "~/signatures/simulations/synthetic_datasets_3107/"
fits_path = paste0("~/signatures/simulations/", "fits_dn.",
                   inference_type, ".", args[2], ".", args[3], ".", regularizer, ".old_hier.", run_id, "/")

cat(paste0("\nSaving in directory: ", fits_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

library(lsa)

cli::cli_process_done()

N = c(150, 500, 1000)
G = c(1, 3, 6)
fracs_rare = 1.

shared_sbs = c("SBS1","SBS5")
private_sbs = c("SBS4","SBS13","SBS10b","SBS7c","SBS7d",
                "SBS11","SBS90","SBS31","SBS10a","SBS22")
catalogue_sbs = COSMIC_filt[c(shared_sbs, private_sbs),]

set.seed(1234)
comb = expand.grid(N_vals=N, n_groups_vals=G, fracs_rare=fracs_rare) %>%
  dplyr::arrange(N_vals)

comb_i = comb[i+1,]
inference_type = c(strsplit(inference_type, "_")[[1]])


# Run model #####
generate_and_run(comb_matrix = comb_i,
                 py = py,

                 fits_path = fits_path,
                 data_path = data_path,
                 seeds = 1:30,

                 catalogue_sbs = catalogue_sbs,
                 alpha_range = c(.15,0.2),
                 alpha_sigma = 0.1,
                 pi_conc = 1.,
                 frac_rare = comb_i$fracs_rare,
                 n_muts_range = 100:5000,
                 shared_sbs = shared_sbs,

                 ## inference
                 reference_catalogue = COSMIC_filt,
                 subset_reference = c("SBS1", "SBS5"),
                 keep_sigs = c("SBS1", "SBS5"),
                 hyperparameters = list("alpha_sigma"=0.1),
                 lr = 0.005,
                 n_steps = 2000,
                 nonparametric = nonparametric,
                 enforce_sparsity = enforce_sparsity,

                 reg_weight = reg_weight,
                 regularizer = regularizer,
                 new_hier = new_hier,
                 regul_denovo = FALSE,
                 regul_fixed = FALSE,

                 initializ_seed = FALSE,
                 do_initial_fit = TRUE,
                 save_runs_seed = TRUE,
                 save_all_fits = TRUE,
                 seed_list = c(4,17,22),

                 CUDA = TRUE,
                 do.fits = TRUE,
                 verbose = FALSE,
                 cohort = i,

                 check_present = TRUE,
		 check_linear_comb = TRUE,
                 inference_type = inference_type
                 )

