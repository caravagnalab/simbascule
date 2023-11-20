args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
inference_type = args[2]
run_id = args[3]

cat(paste("i =", i, "inference_type =", inference_type, "\n"))

main_path = "~/GitHub/"
data_path = "~/signatures/real_data/"

# Load packages #####

cli::cli_process_start("Loading packages")

reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

counts_all = readRDS(paste0(data_path, "counts_sbs.Rds"))
groups_all = counts_all$organ

organ_i = unique(groups_all)[i+1]
groups_i = organ_i

k_list = 8:12
g = 10
subset_reference = c("SBS1","SBS5")

set.seed(10)
seed_list = sample(1:100, 10, replace=FALSE)

counts_i = counts_all %>%
  dplyr::filter(organ==organ_i) %>%
  dplyr::select(-organ, -cohort)

if (i == -1) {
  counts_i = counts_all %>%
    dplyr::select(-organ, -cohort)
  groups_i = groups_all
  organs_i = "whole_cohort"
}

x.fit = fit(x=counts_i,
            k=k_list, # n of denovo signatures
            clusters=g,  # n of clusters
            n_steps=3000, lr=0.005,
            enforce_sparsity=TRUE,
            dirichlet_prior=TRUE,
            reference_catalogue=COSMIC_filt[subset_reference,],
            hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=10000,
                                 "scale_factor_centroid"=10000, "scale_tau"=0),  # change default values to hyperparameters
            nonparametric=TRUE,
            py=py,
            CUDA=TRUE,
            seed_list=seed_list)

lc = filter_signatures_QP(sign1=get_signatures(x.fit),
                          sign2=COSMIC_filt,
                          return_weights=FALSE,
                          filt_pi=0.1)
new_sigs = unique(c(subset_reference, unlist(lc)))
new_min_k = max(0, k_list[1] - length(setdiff(unlist(lc), subset_reference)))
k_list = new_min_k:k_list[length(k_list)]

print(new_sigs)
print(k_list)

x.fit_new = fit(x=counts_i,
                k=k_list, # n of denovo signatures
                clusters=g,  # n of clusters
                n_steps=3000, lr=0.005,
                enforce_sparsity=TRUE,
                dirichlet_prior=TRUE,
                reference_catalogue=COSMIC_filt[new_sigs,],
                hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=10000,
                                     "scale_factor_centroid"=10000, "scale_tau"=0),  # change default values to hyperparameters
                nonparametric=TRUE, py=py, CUDA=TRUE, seed_list=seed_list)

x.fit$lc_check = x.fit_new

x.fit$groups_true = groups_i

saveRDS(x.fit, paste0(data_path, "fit.real_data.sbs.", organ_i, ".N",
                      nrow(counts_i),
                      ".", run_id, ".Rds"))





