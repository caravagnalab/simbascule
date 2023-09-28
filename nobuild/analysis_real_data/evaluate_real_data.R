devtools::load_all("~/GitHub/simbasilica/")
load_deps()

## Test on subset ####
main_path = "~/Dropbox/shared/2022. Basilica/real_data/"
save_path = paste0(main_path, "results/")
data_path = paste0(main_path, "processed_data/")

tissues = c("Colorectal", "Lung", "Breast")
N = 2000
counts_all = readRDS(paste0(data_path, "counts_all.Rds")) %>% dplyr::filter(organ %in% tissues)
groups_all = counts_all$organ

set.seed(13)
idxs = sample(1:nrow(counts_all), N, replace=F)

counts_n = counts_all[idxs,] %>% dplyr::select(-organ, -cohort)
groups_n = groups_all[idxs]

input_data = tibble::tibble("counts"=list(counts_n), "groupid"=list(groups_n))
# saveRDS(input_data, "~/Dropbox/shared/2022. Basilica/datasets/input.N1500.CRC_LUNG.Rds")


## Dirichlet
# fit_dn = fit(x=input_data$counts[[1]], k=8, clusters=4, nonparametric=TRUE,
#              reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),], n_steps=2000,
#              keep_sigs=c("SBS1","SBS5"), enforce_sparsity=TRUE,
#              hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=10000,
#                                   "scale_factor_centroid"=10000, "scale_tau"=1),
#              dirichlet_prior=TRUE, py=py)

fit_dn = fit(x=input_data$counts[[1]], k=6:8, clusters=10, nonparametric=TRUE,
             reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),], n_steps=2000,
             keep_sigs=c("SBS1","SBS5"), enforce_sparsity=TRUE,
             hyperparameters=list("alpha_conc"=100, "scale_factor_alpha"=10000,
                                  "scale_factor_centroid"=10000, "scale_tau"=0),
             dirichlet_prior=TRUE, py=py)

saveRDS(fit_dn, "~/Dropbox/shared/2022. Basilica/real_data/results/objects/fit.N2000.BR_CRC_LUNG.sf_test.Rds")

fit_dn %>% convert_sigs_names(reference_cat = COSMIC_filt, cutoff = .7) %>%
  plot_exposures_real(groups_true = input_data$groupid[[1]]) %>%
  patchwork::wrap_plots(fit_dn %>% convert_sigs_names(reference_cat = COSMIC_filt, cutoff = .7) %>%
                          plot_centroids(), widths = c(9,1), guides="collect")
ggsave("~/Dropbox/shared/2022. Basilica/real_data/results/objects/expos.BR_CRC_LUNG.sf_test.png")

lc_sbs = filter_signatures_QP(get_signatures(fit_dn %>% convert_sigs_names(reference_cat=COSMIC_filt)),
                              COSMIC_filt, filt_pi=0.1, return_weights=F) %>% unlist() %>% unique()

fit_dn.lc = fit(x=input_data$counts[[1]], k=4:6, clusters=4, nonparametric=TRUE,
                reference_catalogue=COSMIC_filt[lc_sbs,],
                keep_sigs=c("SBS1","SBS5"), enforce_sparsity=TRUE,
                hyperparameters=list("alpha_conc"=1000, "scale_factor_alpha"=10000,
                                      "scale_factor_centroid"=5000, "scale_tau"=5),
                dirichlet_prior=TRUE, py=py)

fit_dn$lc_check = fit_dn.lc
saveRDS(fit_dn, "~/Dropbox/shared/2022. Basilica/real_data/results/objects/fit.N1500.CRC_LUNG.dmm.Rds")


fit_dn %>% convert_sigs_names(reference_cat=COSMIC_filt) %>% plot_exposures_real(groups_true=input_data$groupid[[1]])

fit_dn.lc %>% convert_sigs_names(reference_cat=COSMIC_filt) %>% plot_exposures_real(groups_true=input_data$groupid[[1]])

filter_signatures_QP(get_signatures(fit_dn %>% convert_sigs_names(reference_cat=COSMIC_filt)),
                     COSMIC_filt, filt_pi=0.1, return_weights=F) %>% unlist() %>% unique()

# grps_true = setNames(input_data$groupid[[1]], rownames(input_data$counts[[1]]))
# grps_fit = setNames(fit2$groups, rownames(fit2$input$counts))
#
# aricode::NMI(grps_true, grps_fit[rownames(input_data$counts[[1]])])
#
# fit2 = fit_dn.dir %>% convert_sigs_names(reference_cat = COSMIC_filt) %>%
#   recompute_centroids() %>% merge_clusters()
#
# fit2 %>% plot_mutations()
# ggsave("~/Dropbox/uni/phd/presentation/data_real.pdf", height=4, width=8)
#
# fit2$groups = ifelse(fit2$groups=="0", "Lung", "Skin")
# fit2 %>% filter_exposures(min_expos = 0.03) %>% plot_exposures()
# ggsave("~/Dropbox/uni/phd/presentation/expos_real.pdf", height=3, width=8)
#
# plot_signatures(fit2, signames = c("SBS4","SBS7a","SBS7b"))
# ggsave("~/Dropbox/uni/phd/presentation/sigs_real.pdf", height=5, width=8)

aricode::NMI(fit2$groups, input_data$groupid[[1]])

fit3 = fit2
fit3$groups = input_data$groupid[[1]]
fit3 %>% plot_exposures()

convert_sigs_names()

# fix_assignments(fit_dn) %>% plot_exposures()
# fit_dn = readRDS(paste0(save_path, "fit_CRC_dn.Rds"))
# fit_cat = readRDS(paste0(save_path, "fit_CRC_cat.Rds"))



## Plots ####
fit1 = fit_dn %>% convert_sigs_names(reference_cat=COSMIC_filt, cutoff=.8)
fit2 = fit_dn$lc_check %>% convert_sigs_names(reference_cat=COSMIC_filt, cutoff=.8)
name1 = "fit1"
name2 = "fit2"
idd = "N1500.CRC_LUNG.dmm"

ref_cls = COSMIC_color_palette()[unique(c(get_signames(fit1),get_signames(fit2)))] %>% purrr::discard(is.na)
dn_cls = c(get_color_palette(fit1),
           get_color_palette(fit2))[grep("D",c(get_signames(fit1),
                                               get_signames(fit2)) %>% unique(), value=TRUE)]
cls = c(ref_cls, dn_cls)
pp = make_plots_compare(fit1=fit1, fit2=fit2,
                        name1=name1, name2=name2,
                        min_exposure=.0, cls=cls)

centr1 = lapply(unique(fit1$groups), function(gid) {
  idxs = get_group(fit1, groupIDs=gid, return_idx=TRUE)
  if (length(idxs) == 0) next
  plot_exposures(fit1, sampleIDs=idxs, cls=cls) + labs(title="")
})
centr1[["centr"]] = plot_centroids(fit1, cls=cls) + labs(title="")

centr2 = lapply(unique(fit2$groups), function(gid) {
  idxs = get_group(fit2, groupIDs=gid, return_idx=TRUE)
  if (length(idxs) == 0) next
  plot_exposures(fit2, sampleIDs=idxs, cls=cls) + labs(title="")
})
centr2[["centr"]] = plot_centroids(fit2, cls=cls) + labs(title="")

pdf(paste0(save_path, "plots.", idd, ".pdf"), height=12, width=16)
plot_fit(fit1, fit2, cls=cls, name1=name1, name2=name2) %>% print()
plot_exposures_real(fit1, groups_true=input_data$groupid[[1]], cls=cls) %>%
  patchwork::wrap_plots(plot_exposures_real(fit2, groups_true=input_data$groupid[[1]], cls=cls), ncol=1)
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



### whole cohort ####
x.real = readRDS("~/Dropbox/shared/2022. Basilica/real_data/results/fit.real_data.N18764.tris.dmm.1809.Rds")

x.real$color_palette = COSMIC_color_palette(get_signatures(x.real))
# x.real$groups = get_obj_initial_params(x.real)$groups

x.real.merged = x.real %>% recompute_centroids() %>% merge_clusters(cutoff = .9)

lapply(x.real.merged$groups %>% unique(), function(gid) {
  idxs = get_group(x.real.merged, groupIDs=gid, return_idx=TRUE)
  x.real.merged %>%
    filter_exposures(0.01) %>%
    convert_sigs_names(reference_cat = COSMIC_filt) %>%
    plot_exposures(sampleIDs=idxs)
}) %>% patchwork::wrap_plots(guides="collect")

all_organs = lapply(x.real.merged$groups_true %>% unique(), function(gid) {
  idxs1 = which(x.real.merged$groups_true == gid)
  idxs2 = rownames(get_exposure(x.real.merged))[idxs1]
  x.real.merged %>%
    filter_exposures(0.01) %>%
    convert_sigs_names(reference_cat = COSMIC_filt) %>%
    plot_exposures_real(sampleIDs=idxs2, groups_true=x.real.merged$groups_true)
})

all_organs2 = lapply(x.real.merged$groups %>% unique(), function(gid) {
  # idxs1 = which(x.real$groups_true == gid)
  # idxs2 = rownames(get_exposure(x.real))[idxs1]
  idxs = get_group(x.real.merged, groupIDs=gid, return_idx=T)
  x.real.merged %>%
    filter_exposures(0.01) %>%
    convert_sigs_names(reference_cat = COSMIC_filt) %>%
    plot_exposures_real(sampleIDs=idxs, groups_true=x.real.merged$groups_true)
})

pdf("~/Dropbox/shared/2022. Basilica/real_data/results/plots.N18764.dmm.1809.pdf", width=16, height=10)
# all_organs %>% patchwork::wrap_plots(guides="collect") %>% print()
# all_organs2 %>% patchwork::wrap_plots(guides="collect") %>% print()
all_organs2 %>% print()
dev.off()



## Whole count matrix ####
catalogues = list.files(path=paste0(data_path, "DBS_v1.01/catalogues/"), recursive=TRUE, full.names=TRUE)
counts_all = lapply(catalogues, function(c_i) {
  splitted = strsplit(c_i, "/")[[1]]
  cohort = splitted[length(splitted)-1]
  organ = strsplit(splitted[length(splitted)], "_")[[1]][2]

  return(read.csv(c_i, sep="\t", check.names=FALSE) %>% t() %>% as.data.frame() %>%
           dplyr::mutate(organ=organ, cohort=cohort))
  }) %>% do.call(rbind, .)
saveRDS(counts_all, file=paste0(data_path, "counts_dbs.Rds"))


## Read reference signatures ####
signatures = read.csv(paste0(data_path, "DBS_v1.01/RefSig_DBS_v1.01.tsv"), sep="\t") %>%
  t() %>% as.data.frame()
saveRDS(signatures, file=paste0(data_path, "signatures_dbs.Rds"))


## Whole exposures matrix ####
exposures = list.files(path=paste0(data_path, "DBS_v1.01/organSpecificExposures/"),
                       recursive=TRUE, full.names=TRUE, pattern=".tsv$")
conversion_table = read.csv(paste0(data_path, "DBS_v1.01/RefSig_DBS_conversionMatrix_v1.01.tsv"), sep="\t") %>%
  apply(1, function(sbs) names(sbs)[sbs>0])
referece_expos = readxl::read_xlsx(paste0(data_path, "science.abl9283_tables_s1_to_s33.v2/SupplementaryTables.xlsx"),
                                   sheet="Table S24")

expos_all = lapply(exposures, function(e_i) {
  cat(paste0(e_i, "\n"))
  splitted = strsplit(e_i, "/")[[1]]
  cohortname = splitted[length(splitted)-1]
  organname = strsplit(splitted[length(splitted)], "-")[[1]][2] %>%
    stringr::str_replace_all("_DBS_exposures_finalT.tsv","")

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

saveRDS(expos_all, file=paste0(data_path, "expos_dbs.Rds"))










