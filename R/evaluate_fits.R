
compare_single_fit = function(fitname, fits_path, data_path, cutoff=0.8,
                              filtered_catalogue=TRUE) {
  idd = stringr::str_replace_all(fitname, pattern="fit.|hier.|.Rds", replacement="")
  simulname = paste0("simul.", idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()
  x.fit = readRDS(paste0(fits_path, fitname))
  if (x.fit$n_denovo > 0 && filtered_catalogue)
    x.fit$fit$denovo_signatures = renormalize_denovo_thr(x.fit$fit$denovo_signatures)

  rare_common = rare_common_sigs(x.simul)

  sigs.fit = get_signatures(x.fit)
  sigs.simul = get_signatures(x.simul)

  expos.fit = get_exposure(x.fit)
  expos.simul = get_exposure(x.simul)

  # assigned = compare_sigs_inf_gt(sigs.fit, sigs.simul, cutoff=cutoff)
  # unassigned = c(setdiff(rownames(sigs.fit), assigned),
  #                setdiff(rownames(sigs.simul), names(assigned)))
  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  assigned = assigned_missing$assigned_tp
  unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)

  mse_counts = compute.mse(m_true=x.simul$input$counts,
                           m_inf=get_data(x.fit, reconstructed=T))
  mse_expos = compute.mse(m_true=expos.simul,
                          m_inf=expos.fit,
                          assigned=assigned)

  cosine_sigs = compute.cosine(sigs.fit, sigs.simul, assigned, unassigned, what="sigs")
  cosine_expos = compute.cosine(expos.fit, expos.simul, assigned, unassigned, what="expos")

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

      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,

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


compute.cosine = function(m1, m2, assigned, unassigned, what) {
  if (what == "expos") {
    m1 = as.data.frame(t(m1))
    m2 = as.data.frame(t(m2))
  }

  rownames(m1) = paste0("F_", rownames(m1))
  rownames(m2) = paste0("S_", rownames(m2))

  compare = rbind(m1[paste0("F_",assigned), ],
                  m2[paste0("S_",names(assigned)), ] )

  cosine_sim = lsa::cosine(t(compare))[paste0("F_",assigned),
                                       paste0("S_",names(assigned))]

  cosines_tmp = sapply(names(assigned), function(i)
    cosine_sim[paste0("F_",assigned[i]), paste0("S_",i)]
  )

  return(
    mean(c(cosines_tmp,
           rep(0, length.out=length(unassigned))))
  )
}


convert_sigs_names = function(x.fit, x.simul, cutoff=0.8) {
  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  assigned_missing$assigned_tp = assigned_missing$assigned_tp %>% sort

  signames = rownames(get_signatures(x.fit)) %>% sort
  signames_dn = rownames(get_denovo_signatures(x.fit)) %>% sort
  new_names = c(names(assigned_missing$assigned_tp), assigned_missing$added_fp)
  new_names_dn = assigned_missing$assigned_tp[which(assigned_missing$assigned_tp %in% signames_dn)] %>% names

  # reorder
  x.fit$fit$exposure = x.fit %>% get_exposure() %>% dplyr::select(dplyr::contains(signames))
  x.fit$fit$denovo_signatures = get_denovo_signatures(x.fit)[signames_dn,]
  x.fit$color_palette = x.fit$color_palette[signames]

  # modify names
  colnames(x.fit$fit$exposure) = new_names
  rownames(x.fit$fit$denovo_signatures) = new_names_dn
  names(x.fit$color_palette) = new_names

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




plot_sigs_found = function(stats_df, which=c("rare","common","shared","all"), ratio=F) {
  columns = dplyr::case_when(
    which == "rare" ~ c("n_priv_rare_found", "n_priv_rare"),
    which == "common" ~ c("n_priv_common_found", "n_priv_common"),
    which == "shared" ~ c("n_shared_found", "n_shared"),
    which == "all" ~ c("inf_K", "true_K")
  )

  colors_hier = c("darkorange", "dodgerblue4") %>%
    setNames(c("Hierarchical", "Non-hierarchical"))

  if (ratio)
    p = stats_df %>%
    dplyr::mutate(ratio = columns[1] / columns[2]) %>%
    ggplot() +
    geom_jitter(aes(x=as.factor(N), y=ratio, color=inf_type), size=.5, height=0) +
    geom_violin(aes(x=as.factor(N), y=ratio, color=inf_type), alpha=0) +
    facet_grid(~G, scales="free_x") +
    scale_color_manual(values=colors_hier) +
    theme_bw() + labs(title="N found / N")

  if (!ratio)
    p = stats_df %>%
    ggplot() +
    geom_jitter(aes_string(x=columns[1], y=columns[2], color="inf_type"), size=.5, height=0) +
    geom_violin(aes_string(x=columns[1], y=columns[2], color="inf_type"), alpha=0) +
    ggh4x::facet_nested(~G+as.factor(N)) +
    scale_color_manual(values=colors_hier) + ylim(0, NA) + xlim(0, NA) +
    theme_bw() + labs(title="N found (x) vs N (y)")

  return(p)
}

