args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

inference_type = args[1]
run_id = args[2]

cat(paste("inference_type =", inference_type, "\n"))

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
# input_data = tibble::tibble("counts"=list(counts_all), "groupid"=list(groups_all))

k_list = 15:20
subset_reference = c("SBS1","SBS5")

x.fit = readRDS("/home/elena.buscaroli/signatures/real_data/fit.real_data.N18764.tris.dmm.1809.Rds")

# x.fit = fit(x=counts_all %>% dplyr::select(-organ, -cohort),
#             k=k_list, # n of denovo signatures
#             clusters=35,  # n of clusters
#             n_steps=3000, lr=0.005,
#             enforce_sparsity=TRUE,
#             dirichlet_prior=TRUE,
#             reference_catalogue=COSMIC_filt[subset_reference,],
#             hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=10000,
#                                  "scale_factor_centroid"=5000, "scale_tau"=5),  # change default values to hyperparameters
#             nonparametric=TRUE, py=py, CUDA=TRUE)

lc = filter_signatures_QP(sign1=get_signatures(x.fit),
                          sign2=COSMIC_filt,
                          return_weights=FALSE,
                          filt_pi=0.1)
new_sigs = unique(c(subset_reference, unlist(lc)))
new_min_k = max(0, k_list[1] - length(setdiff(unlist(lc), subset_reference)))
k_list = new_min_k:k_list[length(k_list)]
print(new_sigs)
print(k_list)
x.fit_new = fit(x=counts_all %>% dplyr::select(-organ, -cohort),
            k=k_list, # n of denovo signatures
            clusters=35,  # n of clusters
            n_steps=3000, lr=0.005,
            enforce_sparsity=TRUE,
            dirichlet_prior=TRUE,
            reference_catalogue=COSMIC_filt[new_sigs,],
            hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=10000,
                                 "scale_factor_centroid"=10000, "scale_tau"=0),  # change default values to hyperparameters
            nonparametric=TRUE, py=py, CUDA=TRUE)

x.fit$lc_check = x.fit_new

x.fit$groups_true = groups_all

saveRDS(x.fit, paste0(data_path, "fit.real_data.N",
                      nrow(counts_all),
                      ".", run_id, ".Rds"))





