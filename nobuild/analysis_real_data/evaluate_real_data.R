devtools::load_all("~/GitHub/simbasilica/")
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

# saveRDS(input_data, "~/Dropbox/shared/2022. Basilica/datasets/input.N500_CRC.Rds")


## Dirichlet
fit_dn.dir = fit(x=input.simul$counts[[1]], k=4, clusters=6, nonparametric=TRUE,
                 reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                 keep_sigs=c("SBS1","SBS5"), enforce_sparsity=TRUE,
                 dirichlet_prior=TRUE, verbose=T, py=py)

fit_dn.norm = fit(x=input.simul$counts[[1]], k=4, clusters=6, nonparametric=TRUE,
                  reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                  keep_sigs=c("SBS1","SBS5"), enforce_sparsity=TRUE,
                  dirichlet_prior=FALSE, verbose=T, py=py)



fit_dn.norm %>%
  # filter_exposures(0.05) %>%
  # recompute_centroids() %>% merge_clusters() %>%
  plot_exposures() %>%
  patchwork::wrap_plots(fit_dn.norm %>%
                          # filter_exposures(0.05) %>%
                          # recompute_centroids() %>% merge_clusters() %>%
                          plot_exposures(centroids=T),
                        ncol=2, widths=c(9,1), guides="collect") &
  theme(legend.position="bottom", legend.direction="horizontal")

# fix_assignments(fit_dn) %>% plot_exposures()
idd = "N500.simul"
saveRDS(fit_dn.dir, paste0(save_path, "fit_dn.dir.", idd, ".Rds"))
saveRDS(fit_dn.norm, paste0(save_path, "fit_dn.norm.", idd, ".Rds"))

# fit_dn = readRDS(paste0(save_path, "fit_CRC_dn.Rds"))
# fit_cat = readRDS(paste0(save_path, "fit_CRC_cat.Rds"))



## Plots ####
fit1 = fit_dn.dir %>% convert_sigs_names(reference_cat=COSMIC_filt, cutoff=.6)
fit2 = fit_dn.norm %>% convert_sigs_names(reference_cat=COSMIC_filt, cutoff=.6)
ref_cls = COSMIC_color_palette()[c(get_signames(fit1),get_signames(fit2))] %>% purrr::discard(is.na)
dn_cls = c(get_color_palette(fit1),
           get_color_palette(fit2))[grep("D",c(get_signames(fit1),
                                               get_signames(fit2)) %>% unique())]
cls = c(ref_cls, dn_cls)
pp = make_plots_compare(fit1=fit1, fit2=fit2,
                        name1="Dirichlet prior", name2="Normal prior",
                        min_exposure=.05, cls=cls)

centr1 = lapply(unique(fit1$groups), function(gid) {
  idxs = get_group(fit1, groupIDs=gid, return_idx=TRUE)
  if (length(idxs) == 0) next
  plot_exposures(fit1, sampleIDs=idxs, cls=cls) + labs(title="")
})
centr1[["centr"]] = plot_exposures(fit1, centroids=T, cls=cls) + labs(title="")

centr2 = lapply(unique(fit2$groups), function(gid) {
  idxs = get_group(fit2, groupIDs=gid, return_idx=TRUE)
  if (length(idxs) == 0) next
  plot_exposures(fit2, sampleIDs=idxs, cls=cls) + labs(title="")
})
centr2[["centr"]] = plot_exposures(fit2, centroids=T, cls=cls) + labs(title="")


idd = "alpha_logprog.N500_simul"
pdf(paste0(save_path, "plots.", idd, ".pdf"), height=12, width=16)
plot_fit(fit1, fit2, cls=cls, name1="Dirichlet prior", name2="Normal prior") %>% print()
pp$expos_centr %>% print()
patchwork::wrap_plots(centr1, guides="collect") %>% print()
patchwork::wrap_plots(centr2, guides="collect") %>% print()
patchwork::wrap_plots(pp$umap,
                      plot_gradient_norms(fit1),
                      plot_gradient_norms(fit2), ncol=2) %>% print()
plot_posterior_probs(fit1)
plot_posterior_probs(fit2)
dev.off()


tmp %>% plot_signatures()

tmp = fit_dn.dir.alpha_logprog.N500_CRC %>% recompute_centroids() %>% # %>% merge_clusters() %>%
  convert_sigs_names(reference_cat = COSMIC_filt, cutoff = .7)
tmp %>% filter_exposures(min_expos = 0.05) %>% plot_exposures() %>%
  patchwork::wrap_plots(plot_exposures(tmp %>% filter_exposures(min_expos = 0.05), centroids = T), widths = c(9,1),
                        guides="collect")



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










