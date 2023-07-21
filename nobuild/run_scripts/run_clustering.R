args = commandArgs(trailingOnly = TRUE)
cat(paste("Arguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
# new_hier = args[2] == "new_hier"  # one between "new_hier" and "old_hier"
# enforce_sparsity = args[3] == "sparsity"
nonparametric = args[2] == "nonparametric"
run_id = args[3]

new_hier = FALSE
enforce_sparsity = TRUE
regularizer = "noreg"
inference_type = "clust"
reg_weight = 0.

cat(paste("i =", i, "new_hier =", new_hier, "regularizer =", regularizer, "enforce_sparsity=", enforce_sparsity, 
          "inference_type =", inference_type, "nonparametric =", args[2], "\n"))

main_path = "/home/elena.buscaroli/GitHub/"
data_path = paste0(main_path, "simbasilica/nobuild/simulations/synthetic_datasets_0507/")
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
G = c(1, 3, 5, 5)
samples_per_group = list(15:150, 50:500, 100:1000, 500:5000) %>% setNames(paste0("N",N))

shared_sbs = c("SBS1","SBS5")
priv_comm = paste0("SBS", c(4, 6, 13, 25))
priv_rare = "SBS7c" %>% setNames("SBS7d")
private = list("common"=priv_comm, "rare"=priv_rare)

set.seed(1234)
comb = tibble::tibble(N_vals=N, n_groups_vals=G) %>%
  tidyr::complete(N_vals, n_groups_vals) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(samples_per_group=list(samples_per_group[[paste0("N", N_vals)]]),
                shared=list(shared_sbs), rare=list(priv_rare),
                common=list(c("SBS7d", sample(priv_comm, size=n_groups_vals-1, replace=F))))


comb_i = comb[i+1,]
inference_type = c(strsplit(inference_type, "_")[[1]])

# Run model #####
generate_and_run(catalogue = COSMIC_filt,
                 comb_matrix = comb_i,
                 py = py,
                 inference_type = inference_type,
                 private_fracs = list("rare"=0.05, "common"=0.3),
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

                 enforce_sparsity = TRUE,
                 nonparametric = nonparametric,

                 reg_weight = 0.,
                 regularizer = "noreg",
                 new_hier = FALSE,

                 initializ_seed = FALSE,
                 initializ_pars_fit = FALSE,
                 save_runs_seed = TRUE,
                 seed_list = c(10),
              #    seed_list = c(10, 33, 92),

                 CUDA = TRUE,
                 do.fits = TRUE,
                 verbose = FALSE,
                 new_model = TRUE,
                 check_present = FALSE,
                 cohort = "")

