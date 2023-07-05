main_path = "~/GitHub/"
data_path = paste0(main_path, "simbasilica/nobuild/simulations/synthetic_datasets_0507/")

cat(paste0("\nSaving in directory: ", data_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

devtools::load_all(paste0(main_path, "simbasilica"))
devtools::load_all(paste0(main_path, "basilica"))

cli::cli_process_done()


# Matrix with combs #####

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


# comb = tibble(
#   N_vals = c(50, 300, 300, 1000, 1000, 5000, 5000),
#   n_groups_vals = c(2, 2, 4, 4, 6, 6, 10),
#   samples_per_group = list(25:50, 50:300, 50:300, 100:1000, 100:1000, 200:5000, 200:5000),
#   n_priv_comm = c(2, 2, 4, 4, 5, 5, 5),
#   n_priv_rare = c(1, 1, 1, 1, 2, 2, 3)
# )

# shared = c("SBS1", "SBS5", "SBS17b")
# private = c("SBS10b", "SBS28", "SBS56 SBS10a", "SBS90", "SBS2", "SBS13", "SBS20", "SBS22")


# Run model #####
generate_and_run(comb_matrix = comb,
                 # shared = shared,
                 # private = private,
                 catalogue = COSMIC_filt,
                 py = NULL,
                 private_fracs = list("rare"=0.05, "common"=0.3),
                 data_path = data_path,

                 seeds = 1:20,
                 mut_range = 10:8000,
                 do.fits = FALSE,
                 cohort = "")

