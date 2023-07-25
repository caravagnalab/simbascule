args = commandArgs(trailingOnly = TRUE)
i = as.integer(args[1])

print(i)

base_path = "~/GitHub/"
devtools::load_all(paste0(base_path, "simbasilica"))
load_deps(base_path=base_path)

out_path = paste0(base_path, "simbasilica/nobuild/poster_bits/fits_gel/")

tissues = c("Colorectal", "Lung")
N = 100
counts = readRDS(paste0(base_path, "simbasilica/nobuild/poster_bits/fits_gel/input_data.rds")) %>%
  dplyr::filter(cohort == "GEL", organ %in% tissues) %>% dplyr::select(-cohort, -organ)
groups = readRDS(paste0(base_path, "simbasilica/nobuild/poster_bits/fits_gel/input_data.rds")) %>%
  dplyr::filter(cohort == "GEL", organ %in% tissues) %>% dplyr::pull(organ) %>%
  map_groups()

set.seed(13)
idxs = sample(1:nrow(counts), N, replace=F)

counts_n = counts[idxs,]
groups_n = groups[idxs]

if (i == 0) {
  fit.nh = two_steps_inference(x=counts_n, k=0:15,
                               groups=NULL,
                               py=py,
                               seed_list=c(12,22),
                               CUDA=F,
                               reference_catalogue = COSMIC_filt_merged,
                               filtered_catalogue = TRUE,
                               reg_weight = 1.,
                               reg_bic = TRUE,
                               regularizer = "cosine")
  saveRDS(fit.nh$tot, paste0(out_path, "fit_gel.CRC_L_", N, ".nh.Rds"))
} else if (i == 1) {
  fit.h = two_steps_inference(x=counts_n, k=0:15,
                              groups=groups_n, py=py,
                              seed_list=c(12,22,43),
                              CUDA=F,
                              reference_catalogue = COSMIC_filt_merged,
                              filtered_catalogue = TRUE,
                              reg_weight = 1.,
                              reg_bic = TRUE,
                              regularizer = "cosine")$tot
  saveRDS(fit.h, paste0(out_path, "fit_gel.CRC_L_", N, ".h.Rds"))
}



