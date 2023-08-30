devtools::load_all()
load_deps()

main_path = "~/Dropbox/shared/2022. Basilica/real_data/"
save_path = paste0(main_path, "results/")
data_path = paste0(main_path, "processed_data/")

tissues = c("Colorectal")
N = 500
counts_all = readRDS(paste0(data_path, "counts_all.Rds")) %>% dplyr::filter(organ %in% tissues)
groups_all = counts_all$organ

set.seed(13)
idxs = sample(1:nrow(counts_all), N, replace=F)

counts_n = counts_all[idxs,] %>% dplyr::select(-organ, -cohort)
groups_n = groups_all[idxs]


fit_dn = fit(x=counts_n, k=0:15, clusters=10, nonparametric=TRUE, keep_sigs=c("SBS1","SBS5"),
             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),], reg_weight=0.)
fit_cat = fit(x=counts_n, k=0:15, clusters=10, nonparametric=TRUE, keep_sigs=c("SBS1","SBS5"),
              reference_catalogue=COSMIC_filt, reg_weight=0.)




## Generate whole count matrix ####
catalogues = list.files(path=paste0(data_path, "SBS_v2.03/catalogues/"), recursive=TRUE, full.names=TRUE)
counts_all = lapply(catalogues, function(c_i) {
  splitted = strsplit(c_i, "/")[[1]]
  cohort = splitted[length(splitted)-1]
  organ = strsplit(splitted[length(splitted)], "_")[[1]][2]

  return(read.csv(c_i, sep="\t") %>% t() %>% as.data.frame() %>%
           dplyr::mutate(organ=organ, cohort=cohort))
  }) %>% do.call(rbind, .)
saveRDS(counts_all, file=paste0(data_path, "counts_all.Rds"))








