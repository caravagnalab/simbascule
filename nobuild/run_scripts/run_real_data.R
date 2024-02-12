args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])

devtools::load_all("~/GitHub/simbasilica/")
reticulate::use_condaenv(condaenv="basilica-env")
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica/")

data_path = "~/signatures/"

organ_type = dplyr::case_when(i == 1 ~ "Breast",
                              i == 2 ~ "Colorectal",
                              i == 3 ~ "Lung")

alpha_conc = dplyr::case_when(organ_type == "Breast" ~ list("SBS3"=100),
                              .default=list(10))
print(organ_type)

# Generate organ specific catalogues - include sigs found common in serena's paper
if (organ_type == "Breast") {
  sbs_cat = c(paste0("SBS", c(1, 2, 3, 5, 8, 13, 17, 18, 127)),
              paste0("SBS", c(1, 2, 3, 5, 8, 13, 17, 18, 127)),
              paste0("SBS", c(1, 2, 3, 5, 8, 13, 17, 18, 31)))
  dbs_cat = c(paste0("DBS", c(11, 13, 20)),
              paste0("DBS", c(11, 13, 20)),
              paste0("DBS", c(11, 13, 20)))
} else if (organ_type == "Colorectal") {
  sbs_cat = c(paste0("SBS", c(1, 2, 3, 8, 17, 18, 88, 93, 121)),
              paste0("SBS", c(1, 5, 17, 18, 93, 121)),
              paste0("SBS", c(1, 2, 13, "17a", "17b", 18, 35, 88, 157, 93, 121)))
  dbs_cat = c(paste0("DBS", c(2, 20)),
              paste0("DBS", c(13, 20)),
              paste0("DBS", c(2, 5, 20)))
} else if (organ_type == "Lung") {
  sbs_cat = c(paste0("SBS", c(2, "4a", "4b", 1, 5, 3, 8, 13, 17, 92)),
              paste0("SBS", c(2, "4a", "4b", 1, 5, 3, 8, 13, 17, 18, 92)),
              paste0("SBS", c(2, "4a", "4b", 1, 5, 3, 8, 13, 17, 31, 92)))
  dbs_cat = c(paste0("DBS", c(2, 13, 20)),
              paste0("DBS", c(2, 13)),
              paste0("DBS", c(13, 20, 2, 5)))
}

catalogue_sbs = intersect(rownames(COSMIC_sbs_filt), sbs_cat)
catalogue_dbs = intersect(rownames(COSMIC_dbs), dbs_cat)

# Generate organ specific counts matrices
input_sbs = readRDS(paste0(data_path, "real_data/counts_sbs.Rds")) %>%
  dplyr::filter(organ==organ_type) %>%
  dplyr::select(-organ, -cohort)

input_dbs = readRDS(paste0(data_path, "real_data/counts_dbs.Rds")) %>%
  dplyr::filter(organ==organ_type) %>%
  dplyr::select(-organ, -cohort)

input_sbs = input_sbs[rowSums(input_sbs) > 0, ]
input_dbs = input_dbs[rowSums(input_dbs) > 0, ]
sample_ids = intersect(rownames(input_sbs), rownames(input_dbs))

counts = list("SBS"=input_sbs[sample_ids, colnames(COSMIC_sbs_filt)],
              "DBS"=input_dbs[sample_ids, colnames(COSMIC_dbs)])

reference_cat = list("SBS"=COSMIC_sbs_filt[catalogue_sbs,],
                     "DBS"=COSMIC_dbs[catalogue_dbs,])

x.real.0 = fit(counts=counts,
               k_list=0:25,
               cluster=15,
               reference_cat=reference_cat,
               n_steps=3000, lr=0.005,
               seed_list=c(19,255,18321,331),
               hyperparameters=list("penalty_scale"=0,
                                    "alpha_conc"=alpha_conc),
               store_fits=TRUE,
               py=py, CUDA=TRUE, autoguide=TRUE)

saveRDS(x.real.0, file=paste0("~/signatures/real_data/matched.2011/fits_1102/fit_wcat_penalty0.", organ_type, ".Rds"))

# x.real.lN = fit(counts=counts,
#                 k_list=0:20,
#                 # cluster=10,
#                 reference_cat=reference_cat,
#                 n_steps=3000, lr=0.005,
#                 seed_list=c(19,29,222),
#                 hyperparameters=list("penalty_scale"=nrow(input_sbs)),
#                 store_fits=TRUE,
#                 py=py, CUDA=TRUE)
#
# saveRDS(x.real.lN, file=paste0("~/signatures/real_data/matched.2011/fits_1412/fit_wcat_penaltyLN.", organ_type, ".Rds"))




