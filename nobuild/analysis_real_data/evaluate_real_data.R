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

input_data = tibble::tibble("counts"=list(counts_n), "groupid"=list(groups_n))

saveRDS(input_data, "~/Dropbox/shared/2022. Basilica/datasets/input.N500_CRC.Rds")


## Dirichlet
fit_dn = fit(x=counts_n, k=15, clusters=7, nonparametric=TRUE, keep_sigs=c("SBS1","SBS5"),
             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),], reg_weight=0., verbose=T,
             enforce_sparsity=TRUE, py=py) #, hyperparameters=list("alpha_sigma"=0.05))

fit_cat = fit(x=counts_n, k=10:15, clusters=10, nonparametric=TRUE, keep_sigs=c("SBS1","SBS5"),
              reference_catalogue=COSMIC_filt, reg_weight=0., verbose=T,
              enforce_sparsity=TRUE, py=py)


fit_dn %>%
  # filter_exposures(0.05) %>%
  # recompute_centroids() %>% merge_clusters() %>%
  plot_exposures() %>%
  patchwork::wrap_plots(fit_dn %>%
                          # filter_exposures(0.05) %>%
                          # recompute_centroids() %>% merge_clusters() %>%
                          plot_exposures(centroids=T),
                        ncol=2, widths=c(9,1), guides="collect") & theme(legend.position="bottom")

# fix_assignments(fit_dn) %>% plot_exposures()
idd = "norm_spars"
saveRDS(fit_dn, paste0(save_path, "fit_CRC_dn.", idd, ".Rds"))
saveRDS(fit_cat, paste0(save_path, "fit_CRC_cat.", idd, ".Rds"))


fit_dn = readRDS(paste0(save_path, "fit_CRC_dn.Rds"))
fit_cat = readRDS(paste0(save_path, "fit_CRC_cat.Rds"))

fit_dn = fit_CRC_dn.dirich_spars

samples = get_group(fit_dn, groupIDs = c("2","3","5","6"), return_idx = T)
fit_dn %>% # convert_sigs_names(reference_cat=COSMIC_filt) %>%
  # merge_clusters() %>%
  # filter_exposures() %>%
  plot_exposures(sampleIDs = samples) %>%
  patchwork::wrap_plots(fit_dn %>% # convert_sigs_names(reference_cat=COSMIC_filt) %>%
                          # merge_clusters() %>%
                          plot_exposures(centroids = T))
fit_dn %>% recompute_centroids() %>% merge_clusters() %>% plot_exposures()
fix_assignments()



fit_cat_mod = fit_cat %>% convert_sigs_names(reference_cat = COSMIC_filt_merged) %>%
  recompute_centroids() %>% merge_clusters()

filter_signatures_QP(get_signatures(fit_cat_mod), COSMIC_filt, filt_pi=0.1, return_weights=T)


## Generate whole count matrix ####
catalogues = list.files(path=paste0(data_path, "SBS_v2.03/catalogues/"), recursive=TRUE, full.names=TRUE)
counts_all = lapply(catalogues, function(c_i) {
  splitted = strsplit(c_i, "/")[[1]]
  cohort = splitted[length(splitted)-1]
  organ = strsplit(splitted[length(splitted)], "_")[[1]][2]

  return(read.csv(c_i, sep="\t", check.names=FALSE) %>% t() %>% as.data.frame() %>%
           dplyr::mutate(organ=organ, cohort=cohort))
  }) %>% do.call(rbind, .)
saveRDS(counts_all, file=paste0(data_path, "counts_all.Rds"))


## Read reference signatures ####
signatures = read.csv(paste0(data_path, "SBS_v2.03/RefSig_SBS_v2.03.tsv"), sep="\t") %>%
  t() %>% as.data.frame()
saveRDS(signatures, file=paste0(data_path, "signatures_all.Rds"))


## Generate whole exposures matrix ####
exposures = list.files(path=paste0(data_path, "SBS_v2.03/organSpecificExposures/"),
                       recursive=TRUE, full.names=TRUE, pattern=".tsv$")
conversion_table = read.csv(paste0(data_path, "SBS_v2.03/RefSig_SBS_conversionMatrix_v2.03.tsv"), sep="\t") %>%
  apply(1, function(sbs) names(sbs)[sbs>0])
referece_expos = readxl::read_xlsx(paste0(data_path, "science.abl9283_tables_s1_to_s33.v2/SupplementaryTables.xlsx"), sheet="Table S23")

expos_all = lapply(exposures, function(e_i) {
  cat(paste0(e_i, "\n"))
  splitted = strsplit(e_i, "/")[[1]]
  cohortname = splitted[length(splitted)-1]
  organname = strsplit(splitted[length(splitted)], "-")[[1]][2] %>%
    stringr::str_replace_all("_SBS_exposures_finalT.tsv","")

  expos_file = read.csv(e_i, sep="\t", check.names=F) %>%
    tibble::rownames_to_column(var="sample")
  ref_sbs = intersect(colnames(expos_file), names(conversion_table))
  ref_expos_i = referece_expos %>% dplyr::filter(cohort==cohortname, organ==organname) %>%
    dplyr::select(sample, unlist(conversion_table[ref_sbs]) %>% setNames(NULL))

  conv_expos = expos_file %>% dplyr::select(-dplyr::all_of(ref_sbs)) %>%
    dplyr::full_join(ref_expos_i, by="sample") %>% tibble::column_to_rownames(var="sample") %>%
    dplyr::mutate(organ=organname, cohort=cohortname)

  assertthat::assert_that(sum(is.na(conv_expos))==0)

  return(conv_expos)
}) %>% dplyr::bind_rows() %>% replace(is.na(.), 0)

saveRDS(expos_all, file=paste0(data_path, "expos_all.Rds"))








