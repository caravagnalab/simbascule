get_stats_df = function(data_path, fits_path, cutoff=0.8, min_exposure=0.,
                        data_pattern=c("simul."),
                        fits_pattern=c("fit.","fit_hier.","fit_clust.")) {

  fits_pattern = fits_pattern %>% stringr::str_replace_all("\\.","\\\\.")

  return(
    lapply(fits_path, function(path_i) {

      lapply(fits_pattern, function(pattern_i) {
        fits_i = list.files(path=path_i, pattern=pattern_i)

        lapply(fits_i, function(fitname) {
          print(paste(fitname, pattern_i))
          compare_single_fit(fitname=fitname, fits_path=path_i, data_path=data_path,
                             data_pattern=data_pattern, fits_pattern=pattern_i,
                             cutoff=cutoff, min_exposure=min_exposure) %>%
            dplyr::mutate(inf_type=pattern_i,
                          fits_path=path_i)
        } ) %>% do.call(what=rbind, args=.)
      } ) %>% do.call(what=rbind, args=.)
    }) %>% do.call(what=rbind, args=.) %>%
      dplyr::mutate(inf_type=stringr::str_replace_all(inf_type, "\\\\.",""))
  )
}


compare_single_fit = function(fitname, fits_path, data_path, cutoff=0.8,
                              fits_pattern, data_pattern="simul.",
                              filtered_catalogue=TRUE, min_exposure=0.) {
  idd = stringr::str_replace_all(fitname, pattern=paste0(fits_pattern,"|.Rds"), replacement="")
  simulname = paste0(data_pattern, idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()
  x.fit = readRDS(paste0(fits_path, fitname))
  if (x.fit$n_denovo > 0 && filtered_catalogue)
    x.fit$fit$denovo_signatures = renormalize_denovo_thr(x.fit$fit$denovo_signatures)

  x.fit = filter_exposures(x.fit, min_exp=min_exposure) %>%
    convert_sigs_names(x.simul, cutoff=cutoff)

  rare_common = rare_common_sigs(x.simul)

  sigs.fit = get_signatures(x.fit)
  sigs.simul = get_signatures(x.simul)

  expos.fit = get_exposure(x.fit)
  expos.simul = get_exposure(x.simul)

  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  assigned = assigned_missing$assigned_tp
  unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)

  assigned_missing_all = get_assigned_missing(x.fit=x.fit, x.simul=x.simul,
                                              cutoff=0.)

  mse_counts = compute.mse(m_true=x.simul$input$counts,
                           m_inf=get_data(x.fit, reconstructed=T))
  mse_expos = compute.mse(m_true=expos.simul,
                          m_inf=expos.fit,
                          assigned_missing=assigned_missing)
  mse_expos_rare = compute.mse(m_true=expos.simul,
                               m_inf=expos.fit,
                               assigned_missing=assigned_missing,
                               subset_cols=rare_common$private_rare)


  cosine_sigs = compute.cosine(sigs.fit, sigs.simul,
                               # assigned_missing=assigned_missing,
                               assigned_missing=assigned_missing_all,
                               what="sigs")

  cosine_expos = compute.cosine(expos.fit, expos.simul,
                                # assigned_missing=assigned_missing,
                                assigned_missing=assigned_missing_all,
                                what="expos")
  cosine_expos_rare = compute.cosine(expos.fit, expos.simul,
                                     # assigned_missing=assigned_missing,
                                     assigned_missing=assigned_missing_all,
                                     what="expos",
                                     subset_cols=rare_common$private_rare)

  if (have_groups(x.fit)) {
    groups_new = get_groups_rare(x.simul, rare_common)

    ari_rare = aricode::ARI(groups_new, x.fit$groups)
    nmi_rare = aricode::NMI(groups_new, x.fit$groups)

    ari = aricode::ARI(x.simul$groups, x.fit$groups)
    nmi = aricode::NMI(x.simul$groups, x.fit$groups)
  } else {
    ari = nmi = ari_rare = nmi_rare = NA
  }

  n_rare_found = sum(assigned %in% rare_common$private_rare)
  n_common_found = sum(assigned %in% rare_common$private_common)
  n_shared_found = sum(assigned %in% rare_common$shared)

  true_n_sigs = nrow(sigs.simul)
  inf_n_sigs = nrow(sigs.fit)

  return(
    tibble::tibble(
      "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][1] %>%
        as.numeric(),
      "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][2] %>%
        as.numeric(),

      "n_shared"=length(rare_common$shared),
      "n_common"=length(rare_common$private_common),
      "n_rare"=length(rare_common$private_rare),

      "n_shared_found"=n_shared_found,
      "n_common_found"=n_common_found,
      "n_rare_found"=n_rare_found,

      "missing_fn"=list(assigned_missing$missing_fn),
      "added_fp"=list(assigned_missing$added_fp),
      "assigned"=list(assigned_missing$assigned_tp),

      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,
      "mse_expos_rare"=mse_expos_rare,

      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "cosine_expos_rare"=cosine_expos_rare,

      "n_sigs"=true_n_sigs,
      "n_sigs_found"=inf_n_sigs,

      "n_groups_found"=length(x.fit$groups %>% unique()),
      "ari_rare"=ari_rare,
      "nmi_rare"=nmi_rare,
      "ari"=ari,
      "nmi"=nmi,

      "shared"=list(rare_common$shared),
      "priv_common"=list(rare_common$private_common),
      "priv_rare"=list(rare_common$private_rare),

      "idd"=idd,
      "is_hierarchical"=grepl("hier", fitname)
    )
  )
}

