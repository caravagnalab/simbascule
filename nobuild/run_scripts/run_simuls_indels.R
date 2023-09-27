args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
inference_type = args[2]
run_id = args[3]

cat(paste("i =", i, "inference_type =", inference_type, "\n"))

main_path = "~/GitHub/"
data_path = "~/signatures/simulations/synthetic_datasets_indels_2609/"
data_path = "~/Dropbox/shared/2022. Basilica/simulations/synthetic_datasets_indels_2609/"
fits_path = paste0("~/signatures/simulations/", "fits_dn.", inference_type, ".", run_id, "/")

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

shared_sbs = c("ID8")
private_sbs = c("ID1",  # correlated with MSI
                "ID2",  # correlated with MSI
                "ID3",  # tobacco smoking
                "ID5",  # no mut process
                "ID6",  # HR deficiency
                "ID7",  # defective DNA mismatch repair
                # "ID8",  # seems to accumulate also in normal cells (prob clock-like)
                "ID13",  # UV light exposure
                "ID17",  # mutations in TOP2A
                "ID18"  # e.coli with pks pathogenicity island
                )
catalogue_sbs = COSMIC_indels[c(shared_sbs, private_sbs),]

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
                 do.fits = FALSE,
                 reference_catalogue = COSMIC_filt,
                 subset_reference = c("SBS1", "SBS5"),
                 keep_sigs = c("SBS1", "SBS5"),
                 hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=5000,
                                      "scale_factor_centroid"=5000, "scale_tau"=0),
                 lr = 0.005,
                 n_steps = 3000,
                 nonparametric = TRUE,
                 enforce_sparsity = TRUE,

                 reg_weight = 0.,
                 regularizer = "cosine",
                 regul_denovo = T,
                 regul_fixed = T,

                 save_all_fits = TRUE,
                 seed_list = c(4,17,22),

                 CUDA = TRUE,
                 cohort = i,

                 check_present = TRUE,
                 check_linear_comb = TRUE,
                 inference_type = inference_type
                 )

