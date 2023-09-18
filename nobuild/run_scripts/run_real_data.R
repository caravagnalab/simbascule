args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
inference_type = args[2]
run_id = args[3]

cat(paste("i =", i, "inference_type =", inference_type, "\n"))

main_path = "~/GitHub/"
data_path = "~/signatures/real_data/"
fits_path = paste0("~/signatures/real_data/", "fits_dn.", inference_type, ".", run_id, "/")

cat(paste0("\nSaving in directory: ", fits_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))


counts_all = readRDS(paste0(data_path, "counts_all.Rds"))
groups_all = counts_all$organ
input_data = tibble::tibble("counts"=list(counts_all), "groupid"=list(groups_all))

x.fit = fit(x=input.real$counts[[1]],
            k=15:20, # n of denovo signatures
            clusters=35,  # n of clusters
            n_steps=3000, lr=0.005,
            enforce_sparsity=TRUE, # using Beta as prior for alpha centroids
            dirichlet_prior=TRUE,
            reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
            hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=1000,
                                 "scale_factor_centroid"=1000, "scale_tau"=1),  # change default values to hyperparameters
            nonparametric=TRUE, py=py, CUDA=TRUE)

saveRDS(x.fit, paste0(fits_path, "fit.real_data.N",
                      nrow(input.real$counts[[1]]),
                      ".", run_id, ".Rds"))





