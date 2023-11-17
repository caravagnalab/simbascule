args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
inference_type = args[2]
run_id = args[3]

cat(paste("i =", i, "inference_type =", inference_type, "\n"))

main_path = "~/GitHub/"
data_path = "~/signatures/simulations/synthetic_datasets_3107/"
fits_path = paste0("~/signatures/simulations/", "fits_dn.", inference_type, ".", run_id, "/")

cat(paste0("\nSaving in directory: ", fits_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

source("~/GitHub/simbasilica/nobuild/run_scripts/fn_run.R")
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
