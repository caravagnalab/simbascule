get_true_expos = function(base_path, organ, long=TRUE, mod_names=FALSE) {
  dpath = paste0("processed_data/SBS_v2.03/organSpecificExposures/GEL/GEL-",
                 organ, "_SBS_exposures_finalT.tsv")

  ff = read.csv(paste0(base_path, dpath),
                sep="\t", row.names=1) %>%
    dplyr::select(-unassigned)

  if (!long && !mod_names) return(ff)

  if (mod_names) {
    colnames(ff) = stringr::str_replace_all(colnames(ff),
                                            paste0("GEL.",organ,"_common_"), "")
    colnames(ff) = stringr::str_replace_all(colnames(ff),
                                            paste0("GEL.",organ,"_rare_"), "")
  }

  if (!long) return(ff)

  return(
    ff %>%
      tibble::rownames_to_column(var="samples") %>%
      reshape2::melt(id="samples", variable.name="sbs", value.name="exposure") %>%
      # tidyr::separate("sbs", into=c("rm","type","sbs"), extra="merge", sep="_") %>%
      dplyr::mutate(organ=organ)
  )
}



get_true_sigs = function(base_path, organ, mod_names=F, long=T) {
  ff = read.csv(paste0(base_path,
                  "processed_data/SBS_v2.03/OrganSpecificSigs_GEL_SBS_v2.03.tsv"),
           sep="\t") %>%
    dplyr::select(dplyr::contains(organ))

  if (mod_names) {
    colnames(ff) = stringr::str_replace_all(colnames(ff),
                                            paste0("GEL.",organ,"_common_"), "")
    colnames(ff) = stringr::str_replace_all(colnames(ff),
                                            paste0("GEL.",organ,"_rare_"), "")
  }

  if (long)
    ff = ff %>%
      tibble::rownames_to_column(var="contexts") %>%
      reshape2::melt(id="contexts", variable.name="sbs", value.name="beta") %>%
      dplyr::mutate(organ=organ)

  return(ff)
}



create_basilica_obj_real_data = function(base_path, counts, groupids) {
  sampleids = counts %>% rownames

  true_expos = get_true_expos(base_path, "Lung", mod_names=T) %>%
    rbind(get_true_expos(base_path, "Colorectal", mod_names=T)) %>%
    dplyr::filter(samples %in% sampleids) %>%
    dplyr::group_by(samples, organ) %>%
    dplyr::mutate(exposure=exposure / sum(exposure)) %>% dplyr::ungroup()

  true_signames = true_expos$sbs %>% as.character() %>% unique()
  lung_signames = true_expos %>% dplyr::filter(organ=="Lung") %>% dplyr::pull(sbs) %>%
    as.character() %>% unique()
  crc_signames = true_expos %>% dplyr::filter(organ=="Colorectal") %>% dplyr::pull(sbs) %>%
    as.character() %>% unique()
  common = intersect(lung_signames, crc_signames)

  true_sigs = get_true_sigs(base_path, "Lung", mod_names=T) %>%
    rbind(get_true_sigs(base_path, "Colorectal", mod_names=T)) %>%
    dplyr::filter(sbs %in% true_signames)

  true_sigs_wide = get_true_sigs(base_path, "Lung", mod_names=F, long=F) %>%
    cbind(get_true_sigs(base_path, "Colorectal", mod_names=F, long=F)) %>%
    dplyr::select(dplyr::contains(true_signames))


  COSMIC = renormalize_denovo_thr(COSMIC_catalogue)
  non_ref = setdiff(true_signames, rownames(COSMIC))
  ref = intersect(rownames(COSMIC), true_signames)
  non_ref_common = intersect(non_ref, common)
  non_ref_unique = setdiff(non_ref, common)

  unique_nonref_cat = true_sigs_wide %>%
    dplyr::select(dplyr::contains(non_ref_unique)) %>%
    t() %>% as.data.frame() %>% tibble::rownames_to_column(var="sbs") %>%
    tidyr::separate("sbs", into=c("rm1","rm2","sbs"), sep="_", extra="merge") %>%
    dplyr::select(-dplyr::contains("rm")) %>% tibble::column_to_rownames("sbs")

  common_nonref_cat = true_sigs_wide %>%
    dplyr::select(dplyr::contains(non_ref_common)) %>% t() %>%
    as.data.frame() %>% tibble::rownames_to_column(var="sbs") %>%
    reshape2::melt(id="sbs", variable.name="contexts", value.name="beta") %>%
    tidyr::separate("sbs", into=c("rm1","rm2","sbs"), sep="_", extra="merge") %>%
    dplyr::select(-dplyr::contains("rm")) %>%
    dplyr::group_by(sbs, contexts) %>%
    dplyr::reframe(beta=mean(beta)) %>% dplyr::ungroup() %>%
    tidyr::pivot_wider(id_cols="sbs", values_from="beta", names_from="contexts") %>%
    tibble::column_to_rownames("sbs")

  catalogue = rbind(COSMIC[ref,], unique_nonref_cat, common_nonref_cat) %>%
    renormalize_denovo_thr()
  catalogue["SBS5",] = COSMIC_catalogue["SBS5",]

  catalogue_filt = catalogue[!grepl("\\.|\\_", rownames(catalogue)),]

  true_expos_wide = true_expos %>%
    dplyr::mutate(exposure=replace(exposure, is.na(exposure), 0)) %>%
    dplyr::filter(sbs %in% rownames(catalogue_filt)) %>% dplyr::select(-organ) %>%
    dplyr::group_by(samples) %>%
    dplyr::mutate(exposure=exposure / sum(exposure)) %>%
    tidyr::pivot_wider(id_cols="samples", names_from="sbs",
                       values_from="exposure", values_fill=0) %>%
    tibble::column_to_rownames("samples")

  real_data = list(); class(real_data) = "basilica_obj"
  real_data$input$counts = counts
  real_data$groups = groupids
  real_data$fit$x = real_data$input$counts
  real_data$fit$input_catalogue = catalogue_filt
  real_data$fit$catalogue_signatures = catalogue_filt
  real_data$fit$exposure = true_expos_wide

  return(real_data)
}






create_basilica_obj_sigprofiler = function(fitname, counts=NULL) {
  sigs = read.csv(paste0(fitname, "/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Signatures/COSMIC_SBS96_Signatures.txt"), sep="\t", row.names = 1) %>%
    t() %>% as.data.frame()

  expos = read.csv(paste0(fitname, "/SBS96/Suggested_Solution/COSMIC_SBS96_Decomposed_Solution/Activities/COSMIC_SBS96_Activities.txt"), sep="\t", row.names = 1) %>%
    as.data.frame()
  expos = expos / rowSums(expos)


  obj = list(); class(obj) = "basilica_obj"

  obj$input$counts = counts
  obj$fit$x = counts
  obj$reference_signatures = sigs
  obj$fit$input_catalogue = NULL
  obj$fit$catalogue_signatures = sigs
  obj$fit$exposure = expos

  return(obj)
}







