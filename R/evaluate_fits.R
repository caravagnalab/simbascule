get_stats_df = function(data_path, fits_path, cutoff=0.8) {
  fits = list.files(path=fits_path, pattern="fit.")

  return(
    lapply(fits, function(fitname) {
      print(fitname)
      compare_single_fit(fitname, fits_path, data_path, cutoff=cutoff)
    } ) %>%
      do.call(what=rbind, args=.) %>%
      dplyr::mutate(inf_type=ifelse(is_hierarchical, "Hierarchical", "Non-hierarchical"))
  )
}


compare_single_fit = function(fitname, fits_path, data_path, cutoff=0.8,
                              filtered_catalogue=TRUE) {
  idd = stringr::str_replace_all(fitname, pattern="fit.|hier.|.Rds", replacement="")
  simulname = paste0("simul.", idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()
  x.fit = readRDS(paste0(fits_path, fitname))
  if (x.fit$n_denovo > 0 && filtered_catalogue)
    x.fit$fit$denovo_signatures = renormalize_denovo_thr(x.fit$fit$denovo_signatures)

  rare_common = rare_common_sigs(x.simul)

  x.fit = x.fit %>% convert_sigs_names(x.simul, cutoff=cutoff)

  sigs.fit = get_signatures(x.fit)
  sigs.simul = get_signatures(x.simul)

  expos.fit = get_exposure(x.fit)
  expos.simul = get_exposure(x.simul)

  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  assigned = assigned_missing$assigned_tp
  unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)

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
                               assigned_missing=assigned_missing,
                               what="sigs")
  cosine_expos = compute.cosine(expos.fit, expos.simul,
                                assigned_missing=assigned_missing,
                                what="expos")
  cosine_expos_rare = compute.cosine(expos.fit, expos.simul,
                                     assigned_missing=assigned_missing,
                                     what="expos",
                                     subset_cols=rare_common$private_rare)

  n_rare_found = assigned[names(assigned) %in% rare_common$private_rare] %>% length()
  n_common_found = assigned[names(assigned) %in% rare_common$private_common] %>% length()
  n_shared_found = assigned[names(assigned) %in% rare_common$shared] %>% length()

  true_n_sigs = nrow(sigs.simul)
  inf_n_sigs = nrow(sigs.fit)

  return(
    tibble::tibble(
      "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][1] %>%
        as.numeric(),
      "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][2] %>%
        as.numeric(),
      "shared"=list(rare_common$shared),
      "priv_common"=list(rare_common$private_common),
      "priv_rare"=list(rare_common$private_rare),

      "n_shared"=length(rare_common$shared),
      "n_priv_common"=length(rare_common$private_common),
      "n_priv_rare"=length(rare_common$private_rare),

      "n_shared_found"=n_shared_found,
      "n_priv_common_found"=n_common_found,
      "n_priv_rare_found"=n_rare_found,

      "missing_fn"=list(assigned_missing$missing_fn),
      "added_fp"=list(assigned_missing$added_fp),
      "assigned"=list(assigned_missing$assigned_tp),

      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,
      "mse_expos_rare"=mse_expos_rare,

      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "cosine_expos_rare"=cosine_expos_rare,

      "true_K"=true_n_sigs,
      "inf_K"=inf_n_sigs,

      "idd"=idd,
      "is_hierarchical"=grepl("hier", fitname)
    )
  )
}


get_assigned_missing = function(x.fit, x.simul, cutoff=0.8) {
  sigs.fit = get_signatures(x.fit)
  sigs.simul = get_signatures(x.simul)

  assigned = compare_sigs_inf_gt(sigs.fit, sigs.simul, cutoff=cutoff)
  missing = setdiff(rownames(sigs.simul), names(assigned))
  added = setdiff(rownames(sigs.fit), assigned)

  return(list("assigned_tp"=assigned, "missing_fn"=missing, "added_fp"=added))
}


compute.cosine = function(m1, m2, assigned_missing, what, subset_cols=NULL) {
  # m1 = fit
  # m2 = simul

  unassigned = unique(c(assigned_missing$missing_fn,
                        assigned_missing$added_fp))

  if (what == "expos") {
    m1 = as.data.frame(t(m1))
    m2 = as.data.frame(t(m2))
  }

  if (!is.null(subset_cols)) {
    m1 = m1[intersect(rownames(m1), subset_cols), ]
    m2 = m2[intersect(rownames(m2), subset_cols), ]

    if (nrow(m1) == 0) {
      return(0)
    }
    consider = intersect(subset_cols, assigned_missing$assigned_tp)
    unassigned = intersect(unassigned, subset_cols)
  } else
    consider = assigned_missing$assigned_tp

  rownames(m1) = paste0("F_", rownames(m1))
  rownames(m2) = paste0("S_", rownames(m2))

  compare = rbind(m1[paste0("F_",consider), ],
                  m2[paste0("S_",consider), ] )

  cosine_sim = lsa::cosine(t(compare))[paste0("F_",consider),
                                       paste0("S_",consider)]

  if (is.numeric(cosine_sim))
    cosines_tmp = cosines_tmp
  else
    cosines_tmp = sapply(consider, function(i)
      cosine_sim[paste0("F_",i), paste0("S_",i)]
    )

  return(
    mean(c(cosines_tmp,
           rep(0, length.out=length(unassigned))))
  )
}


convert_sigs_names = function(x.fit, x.simul, cutoff=0.8) {
  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  assigned_missing$assigned_tp = assigned_missing$assigned_tp %>% sort

  signames_ref = get_fixed_signames(x.fit)
  signames_dn = get_dn_signames(x.fit)
  which_ref = assigned_missing$assigned_tp %in% signames_ref
  which_dn = assigned_missing$assigned_tp %in% signames_dn

  added_ref = intersect(signames_ref, assigned_missing$added_fp) %>%
    setNames(intersect(signames_ref, assigned_missing$added_fp))
  added_dn = intersect(signames_dn, assigned_missing$added_fp) %>%
    setNames(intersect(signames_dn, assigned_missing$added_fp))

  new_order_ref = c(assigned_missing$assigned_tp[which_ref], added_ref)
  new_order_dn = c(assigned_missing$assigned_tp[which_dn], added_dn)

  # reorder
  x.fit$fit$exposure = x.fit %>% get_exposure() %>%
    dplyr::select(dplyr::contains(new_order_ref), dplyr::contains(new_order_dn))

  x.fit$fit$denovo_signatures = get_denovo_signatures(x.fit)[new_order_dn,]
  x.fit$fit$catalogue_signatures = get_catalogue_signatures(x.fit)[new_order_ref,]
  x.fit$color_palette = x.fit$color_palette[c(new_order_ref, new_order_dn)]

  # modify names
  colnames(x.fit$fit$exposure) = c(new_order_ref, new_order_dn) %>% names
  rownames(x.fit$fit$denovo_signatures) = new_order_dn %>% names
  rownames(x.fit$fit$catalogue_signatures) = new_order_ref %>% names
  names(x.fit$color_palette) = c(new_order_ref, new_order_dn) %>% names

  return(x.fit)
}


compare_sigs_inf_gt = function(sigs.fit, sigs.simul, cutoff=0.8) {
  common = intersect(rownames(sigs.fit), rownames(sigs.simul))
  unique_inf = setdiff(rownames(sigs.fit), common)
  unique_gt = setdiff(rownames(sigs.simul), common)

  total_sigs = rbind(sigs.fit[!rownames(sigs.fit) %in% common,], sigs.simul)
  cosine_matr = lsa::cosine(t(total_sigs))[rownames(sigs.simul), rownames(sigs.fit)]

  assign_similar = cosine_matr %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gt") %>%
    reshape2::melt(id="gt", variable.name="inf", value.name="cosine") %>%
    dplyr::filter(cosine >= cutoff) %>%
    dplyr::group_by(gt) %>%
    dplyr::mutate(inf=as.character(inf)) %>%
    dplyr::filter(cosine == max(cosine)) %>% dplyr::arrange(gt)

  # if (nrow(sigs.simul) > nrow(sigs.fit))
  if (any(duplicated(assign_similar$inf)))
    assign_similar = assign_similar %>% dplyr::group_by(inf) %>%
      dplyr::filter(cosine == max(cosine)) %>% ungroup()

  similar = assign_similar$inf %>% setNames(assign_similar$gt)

  return(similar)
}


rare_common_sigs = function(x.simul) {
  if ("private_sigs" %in% names(x.simul) && length(x.simul$private_sigs) > 0)
    return(list("private_common"=x.simul$private_sigs$private_common,
                "private_rare"=x.simul$private_sigs$private_rare,
                "shared"=rownames(get_fixed_signatures(x.simul))))

  n_grps = length(unique(x.simul$groups))
  expos = get_exposure(x.simul, add_groups=TRUE, long=TRUE)

  rare_comm = expos %>%
    dplyr::mutate(is_present=Exposure > 0) %>%
    dplyr::group_by(Signature, groups) %>%
    dplyr::reframe(sigs_frac=sum(is_present) / length(Sample)) %>%

    dplyr::group_by(Signature) %>%
    dplyr::filter(sigs_frac > 0) %>%
    dplyr::reframe(is_common=any(sigs_frac > 0.1),
                   is_private=length(groups) < n_grps) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(type=dplyr::case_when(
      is_common && is_private ~ "private_common",
      !is_common && is_private ~ "private_rare",
      .default = "shared"
    )) %>% ungroup()

  return(list("private_common"=rare_comm %>%
                dplyr::filter(type=="private_common") %>%
                dplyr::pull(Signature),
              "private_rare"=rare_comm %>%
                dplyr::filter(type=="private_rare") %>%
                dplyr::pull(Signature),
              "shared"=rare_comm %>%
                dplyr::filter(type=="shared") %>%
                dplyr::pull(Signature)))
}



