args = commandArgs(trailingOnly = TRUE)
i = as.integer(args[1])

print(i)

base_path = "/home/elena.buscaroli/GitHub/"
devtools::load_all(paste0(base_path, "simbasilica"))
load_deps(base_path=base_path)

out_path = paste0(base_path, "simbasilica/nobuild/poster_bits/fits_gel/")

tissues = c("Breast", "Colorectal", "Lung")
counts = readRDS(paste0(base_path, "simbasilica/nobuild/poster_bits/fits_gel/input_data.rds")) %>%
  dplyr::filter(cohort == "GEL", organ %in% tissues) %>% dplyr::select(-cohort, -organ)

groups = readRDS(paste0(base_path, "simbasilica/nobuild/poster_bits/fits_gel/input_data.rds")) %>%
  dplyr::filter(cohort == "GEL", organ %in% tissues) %>% dplyr::pull(organ) %>%
  map_groups()


if (i == 0) {
  fit.nh = two_steps_inference(x=counts, k=0:15, groups=NULL, seed_list=c(12,22,43,55), CUDA=TRUE)
  saveRDS(fit.nh, paste0(out_path, "fit_gel.BR_CRC_L.nh.Rds"))
} else if (i == 1) {
  fit.h = two_steps_inference(x=counts, k=0:15, groups=groups, seed_list=c(12,22,43,55), CUDA=TRUE)
  saveRDS(fit.h, paste0(out_path, "fit_gel.BR_CRC_L.h.Rds"))
}



