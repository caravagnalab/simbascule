args = commandArgs(trailingOnly = TRUE)
cat(paste("Arguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
# new_hier = args[2] == "new_hier"  # one between "new_hier" and "old_hier"
# enforce_sparsity = args[3] == "sparsity"
nonparametric = args[2] == "nonparametric"
inference_type = args[3]
run_id = args[4]

new_hier = FALSE
enforce_sparsity = TRUE
regularizer = "noreg"
# inference_type = "clust"
reg_weight = 0.

cat(paste("i =", i, "new_hier =", new_hier, "regularizer =", regularizer, "enforce_sparsity=", enforce_sparsity,
          "inference_type =", inference_type, "nonparametric =", nonparametric, "\n"))

main_path = "~/GitHub/"
data_path = paste0(main_path, "simbasilica/nobuild/simulations/synthetic_datasets_2507/")
fits_path = paste0(main_path, "simbasilica/nobuild/simulations/", "fits_dn.",
                   inference_type, ".", args[2], ".", regularizer, ".old_hier.", run_id, "/")

cat(paste0("\nSaving in directory: ", fits_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

cli::cli_process_done()

N = c(150, 500, 1000, 5000)
G = rep(1, times=length(N))
samples_per_group = list(15:150, 50:500, 100:1000, 500:5000) %>% setNames(paste0("N",N))
fracs_rare = seq(0.05, 0.8, length.out=6) %>% round(digits=2)

shared_sbs = c("SBS1","SBS5")
priv_comm = c()
priv_rare = paste0("SBS", c(4, 6, 13))
# priv_rare = "SBS7c" %>% setNames("SBS7d")
private = list("common"=priv_comm, "rare"=priv_rare)


set.seed(1234)
comb = tibble::tibble(N_vals=N, n_groups_vals=G, fracs_rare=list(fracs_rare)) %>%
  tidyr::unnest(fracs_rare) %>%
  tidyr::complete(N_vals, n_groups_vals) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(samples_per_group=list(samples_per_group[[paste0("N", N_vals)]]),
                shared=list(shared_sbs), rare=list(priv_rare), common=list(c()))
                # common=list(sample(priv_comm, size=n_groups_vals, replace=F)))


comb_i = comb[i+1,]
inference_type = c(strsplit(inference_type, "_")[[1]])

# Run model #####
generate_and_run(catalogue = COSMIC_filt,
                 comb_matrix = comb_i,
                 py = py,
                 inference_type = inference_type,
                 private_fracs = list("rare"=comb_i$fracs_rare, "common"=0.3),
                 fits_path = fits_path,
                 data_path = data_path,

                 seeds = 1:10,
                 mut_range = 10:8000,
                 reference_catalogue = COSMIC_filt[c("SBS1","SBS5"), ],
                 input_catalogue = NULL,
                 keep_sigs = c("SBS1", "SBS5"),
                 hyperparameters = list("alpha_sigma"=0.05),
                 lr = 0.005,
                 n_steps = 2000,

                 enforce_sparsity = enforce_sparsity,
                 nonparametric = nonparametric,

                 reg_weight = reg_weight,
                 regularizer = regularizer,
                 new_hier = new_hier,

                 initializ_seed = FALSE,
                 initializ_pars_fit = FALSE,
                 save_runs_seed = TRUE,
                 save_all_fits = TRUE,
                 #   seed_list = c(10),
                 seed_list = c(10, 33, 92),

                 do.fits = TRUE,
                 verbose = FALSE,
                 check_present = FALSE,
                 cohort = i)

