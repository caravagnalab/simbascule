
get_assigned_missing = function(x.fit, x.simul=NULL, reference_cat=NULL, cutoff=0.8) {
  sigs.fit = get_signatures(x.fit)
  if (!is.null(x.simul)) sigs.simul = get_signatures(x.simul)
  else if (!is.null(reference_cat)) sigs.simul = reference_cat

  assigned = compare_sigs_inf_gt(sigs.fit, sigs.simul, cutoff=cutoff)
  missing = setdiff(rownames(sigs.simul), names(assigned))
  added = setdiff(rownames(sigs.fit), assigned)

  return(list("assigned_tp"=assigned, "missing_fn"=missing, "added_fp"=added))
}


rename_expos = function(exposures, old_names) {
  exposures = exposures %>% dplyr::select(dplyr::contains(old_names))
  colnames(exposures) = old_names %>% names
  return(exposures)
}

rename_sigs = function(signatures, old_names) {
  signatures = signatures[old_names, ]
  rownames(signatures) = old_names %>% names
  return(signatures)
}


rename_params = function(params_list, new_order_ref, new_order_dn) {
  for (parname in names(params_list)) {
    if (is.null(params_list[[parname]])) next
    if (grepl("alpha", parname))
      params_list[[parname]] = rename_expos(exposures=params_list[[parname]],
                                            old_names=c(new_order_ref, new_order_dn))
    else if (grepl("beta_d", parname))
      params_list[[parname]] = rename_sigs(signatures=params_list[[parname]],
                                           old_names=new_order_dn)
    else if (grepl("beta_f", parname))
      params_list[[parname]] = rename_sigs(signatures=params_list[[parname]],
                                           old_names=new_order_ref)
  }
  return(params_list)
}


convert_sigs_names = function(x.fit, x.simul=NULL, reference_cat=NULL, cutoff=0.8) {
  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, reference_cat=reference_cat, cutoff=cutoff)
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

  # rename elements
  x.fit$fit$exposure = rename_expos(exposures=get_exposure(x.fit),
                                    old_names=c(new_order_ref, new_order_dn))
  x.fit$fit$denovo_signatures = rename_sigs(signatures=get_denovo_signatures(x.fit),
                                            old_names=new_order_dn)
  x.fit$fit$catalogue_signatures = rename_sigs(signatures=get_catalogue_signatures(x.fit),
                                            old_names=new_order_ref)
  x.fit$color_palette = x.fit$color_palette[c(new_order_ref, new_order_dn)]

  x.fit$fit$params = rename_params(params_list=x.fit$fit$params,
                                   new_order_ref=new_order_ref,
                                   new_order_dn=new_order_dn)
  x.fit$fit$init_params = rename_params(params_list=x.fit$fit$init_params,
                                        new_order_ref=new_order_ref,
                                        new_order_dn=new_order_dn)

  try(expr = {
    names(x.fit$color_palette) = c(new_order_ref, new_order_dn) %>% names
  }, silent=T)

  return(x.fit)
}


compare_sigs_inf_gt = function(sigs.fit, sigs.simul, cutoff=0.8) {
  common = intersect(rownames(sigs.fit), rownames(sigs.simul))
  unique_inf = setdiff(rownames(sigs.fit), common)
  unique_gt = setdiff(rownames(sigs.simul), common)

  if (length(unique_inf) == 0 || length(unique_gt) == 0)
    return(common %>% setNames(common))

  total_sigs = rbind(sigs.fit[!rownames(sigs.fit) %in% common,],
                     sigs.simul[!rownames(sigs.simul) %in% common,])
  cosine_matr = lsa::cosine(t(total_sigs))[unique_gt, unique_inf]

  if (length(unique_inf) == 1 && length(unique_gt) == 1) {
    cosine_matr = as.data.frame(cosine_matr)
    rownames(cosine_matr) = unique_gt
    colnames(cosine_matr) = unique_inf
  }

  assign_similar = cosine_matr %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gt") %>%
    reshape2::melt(id="gt", variable.name="inf", value.name="cosine") %>%
    dplyr::filter(cosine >= cutoff)

  if (nrow(assign_similar) == 0) return(common %>% setNames(common))

  assign_similar = assign_similar %>%
    dplyr::group_by(gt) %>%
    dplyr::mutate(inf=as.character(inf)) %>%
    dplyr::filter(cosine == max(cosine)) %>% dplyr::arrange(gt)

  # if (nrow(sigs.simul) > nrow(sigs.fit))
  if (any(duplicated(assign_similar$inf)))
    assign_similar = assign_similar %>% dplyr::group_by(inf) %>%
      dplyr::filter(cosine == max(cosine)) %>% ungroup()

  # assigned = assign_similar$inf %>% setNames(assign_similar$gt)
  assigned = c(common, assign_similar$inf) %>% setNames(c(common, assign_similar$gt))

  return(assigned)
}


rare_common_sigs = function(x.simul) {
  if ("private_sigs" %in% names(x.simul) && length(x.simul$private_sigs) > 0)
    return(list("private_common"=x.simul$private_sigs$private_common,
                "private_rare"=x.simul$private_sigs$private_rare,
                "shared"=rownames(get_fixed_signatures(x.simul))))

  if ("sigs" %in% names(x.simul) && length(x.simul$sigs) > 0)
    return(list("private_common"=x.simul$sigs$private_shared,
                "private_rare"=x.simul$sigs$private,
                "shared"=x.simul$sigs$shared))

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




get_groups_rare = function(x.simul, rare_common=NULL) {
  if (is.null(rare_common)) rare_common = rare_common_sigs(x.simul)

  groups_new = x.simul$groups

  sigs_per_group = lapply(unique(x.simul$groups),
                          function(gid)
                            get_sigs_group(x.simul, groupID=gid)) %>%
    setNames(unique(x.simul$groups))

  samples_tmp = lapply(names(sigs_per_group), function(gid) {
    signames = sigs_per_group[[gid]]
    if (any(rare_common$private_rare %in% signames))
      lapply(intersect(rare_common$private_rare, signames),
             function(r) get_samples_with_sigs(x.simul, r, return_idx=TRUE) ) %>%
      setNames(intersect(rare_common$private_rare, signames))
  } ) %>% setNames(names(sigs_per_group)) %>% purrr::discard(is.null)

  for (gid in samples_tmp) for (j in gid) groups_new[j] = max(groups_new)+1

  return(groups_new)
}


compute_ari_nmi = function(x.simul, x.fit) {
  groups_fit = x.fit$groups; groups_simul = x.simul$groups
  if (length(unique(groups_fit)) == 1 || length(unique(groups_simul)) == 1) {
    groups_fit = c(groups_fit, "imolabella")
    groups_simul = c(groups_simul, "imolabella")
  }

  ari = aricode::ARI(as.character(groups_fit), as.character(groups_simul))
  nmi = aricode::NMI(as.character(groups_fit), as.character(groups_simul))
  return(c(ari, nmi))
}
