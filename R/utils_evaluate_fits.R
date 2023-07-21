
get_assigned_missing = function(x.fit, x.simul, cutoff=0.8) {
  sigs.fit = get_signatures(x.fit)
  sigs.simul = get_signatures(x.simul)

  assigned = compare_sigs_inf_gt(sigs.fit, sigs.simul, cutoff=cutoff)
  missing = setdiff(rownames(sigs.simul), names(assigned))
  added = setdiff(rownames(sigs.fit), assigned)

  return(list("assigned_tp"=assigned, "missing_fn"=missing, "added_fp"=added))
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
  try(expr = {
    x.fit$fit$params$alpha_prior = x.fit$fit$params$alpha_prior %>%
      dplyr::select(dplyr::contains(new_order_ref), dplyr::contains(new_order_dn))
    colnames(x.fit$fit$params$alpha_prior) = c(new_order_ref, new_order_dn) %>% names
  }, silent=T)
  try(expr = {
    x.fit$fit$params$alpha_noise = x.fit$fit$params$alpha_noise %>%
      dplyr::select(dplyr::contains(new_order_ref), dplyr::contains(new_order_dn))
    colnames(x.fit$fit$params$alpha_noise) = c(new_order_ref, new_order_dn) %>% names
  }, silent=T)
  try(expr = {
    x.fit$fit$params$beta_d = x.fit$fit$params$beta_d[new_order_dn,]
    rownames(x.fit$fit$params$beta_d) = new_order_dn %>% names
  }, silent=T)
  try(expr = {
    x.fit$fit$params$beta_f = x.fit$fit$params$beta_f[new_order_ref,]
    rownames(x.fit$fit$params$beta_f) = new_order_ref %>% names
  }, silent=T)
  try(expr = {
    names(x.fit$color_palette) = c(new_order_ref, new_order_dn) %>% names
  }, silent=T)
  # names(x.fit$color_palette) = c(new_order_ref, new_order_dn) %>% names

  return(x.fit)
}


compare_sigs_inf_gt = function(sigs.fit, sigs.simul, cutoff=0.8) {
  common = intersect(rownames(sigs.fit), rownames(sigs.simul))
  unique_inf = setdiff(rownames(sigs.fit), common)
  unique_gt = setdiff(rownames(sigs.simul), common)

  if (length(unique_inf) == 0 || length(unique_gt) == 0)
    return(common %>% setNames(common))

  total_sigs = rbind(sigs.fit[!rownames(sigs.fit) %in% common,], sigs.simul)
  cosine_matr = lsa::cosine(t(total_sigs))[rownames(sigs.simul), rownames(sigs.fit)]

  assign_similar = cosine_matr %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gt") %>%
    reshape2::melt(id="gt", variable.name="inf", value.name="cosine") %>%
    dplyr::filter(cosine >= cutoff, !gt%in%common, !inf%in%common) %>%
    dplyr::group_by(gt) %>%
    dplyr::mutate(inf=as.character(inf)) %>%
    dplyr::filter(cosine == max(cosine)) %>% dplyr::arrange(gt)

  if (nrow(assign_similar) == 0) return(common %>% setNames(common))

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
  if (is.null(rare_common)) rare_common = rare_common_sigs(simul)

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

