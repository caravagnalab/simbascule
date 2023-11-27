make_boxplot = function(all_stats, colname) {
  all_stats %>% dplyr::select(N, G, dplyr::contains(colname)) %>%
    reshape2::melt(id=c("N","G"), variable.name="type") %>%
    dplyr::mutate(type=stringr::str_replace_all(type, paste0(colname,"_"),"")) %>%
    ggplot() +
    geom_boxplot(aes(x=as.factor(N), y=value)) +
    ggh4x::facet_nested(type ~ G) +
    theme_bw()
}


eval_single_fit_matched = function(x.fit, x.simul, cutoff=0.8) {
  x.fit = x.fit %>% rename_dn_expos()
  lapply(get_types(x.fit), function(tid) {
    sigs.fit = get_signatures(x.fit, matrix=T)[[tid]]; sigs.simul = get_signatures(x.simul, matrix=T)[[tid]]
    sigs_fixed.fit = get_fixed_signatures(x.fit, matrix=T)[[tid]]
    sigs_dn.fit = get_denovo_signatures(x.fit, matrix=T)[[tid]]

    expos.fit = get_exposure(x.fit, matrix=T)[[tid]]; expos.simul = get_exposure(x.simul, matrix=T)[[tid]]

    assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, type=tid, cutoff=cutoff)
    assigned = assigned_missing$assigned_tp
    unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)

    mse_counts = compute.mse(m_true=get_input(x.simul, matrix=T)[[tid]], m_inf=get_input(x.fit, reconstructed=T, matrix=T)[[tid]])
    mse_expos = compute.mse(m_true=expos.simul, m_inf=expos.fit,
                            assigned_missing=assigned_missing)

    cosine_sigs = compute.cosine(sigs.fit, sigs.simul,
                                 assigned_missing=assigned_missing,
                                 what="sigs")
    cosine_expos = compute.cosine(expos.fit, expos.simul,
                                  assigned_missing=assigned_missing,
                                  what="expos")

    cosine_fixed = list()
    if (!is.null(sigs_dn.fit)) {
      cosine_fixed = lapply(rownames(sigs_dn.fit), function(sid) {
        cosine_sim = lsa::cosine(t(sigs.fit))[sid, rownames(sigs_fixed.fit)]
      }) %>% setNames(rownames(sigs_dn.fit))
    }

    ari_nmi = ari_nmi_km = ari_nmi_km_em = list(NA, NA)
    if (have_groups(x.fit)) {
      # kmeans1 = get_initial_object(x.fit, what="clustering")
      # kmeans2 = get_initial_object(x.fit, what="clustering") %>% run_kmeans()

      ari_nmi = compute_ari_nmi(x.simul=x.simul, x.fit=x.fit)
      # ari_nmi_km1 = compute_ari_nmi(x.simul=x.simul, x.fit=kmeans1)
      # ari_nmi_km2 = compute_ari_nmi(x.simul=x.simul, x.fit=kmeans2)
    }

    res = tibble::tibble(
      "assigned_missing"=list(assigned_missing),

      "sigs_found"=length(assigned),
      # "shared_found"=sum(assigned %in% rare_common$shared),
      # "private_found"=sum(assigned %in% c(rare_common$private_rare, rare_common$private_common)),
      "K_true"=length(get_signames(x.simul)[[tid]]),
      "K_found"=length(get_signames(x.fit)[[tid]]),

      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,

      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "cosine_fixed"=list(cosine_fixed),

      "groups_found"=length(get_cluster_labels(x.fit)),
      "ari"=ari_nmi[[1]],
      "nmi"=ari_nmi[[2]]
      # "ari_km1"=ari_nmi_km1[[1]],
      # "nmi_km1"=ari_nmi_km1[[2]],
      # "ari_km2"=ari_nmi_km2[[1]],
      # "nmi_km2"=ari_nmi_km2[[2]]
    )

    colnames(res) = paste(colnames(res), tid, sep="_")
    return(res)
  }) %>% setNames(get_types(x.fit))
}



make_plots_stats = function(stats) {
  cosine_expos = make_boxplot(stats, "cosine_expos") + labs(title="cosine_expos")
  cosine_sigs = make_boxplot(stats, "cosine_sigs") + labs(title="cosine_sigs")
  mse_counts = make_boxplot(stats, "mse_counts") + labs(title="mse_counts")

  k_ratio = stats %>% dplyr::select(N, G, seed, dplyr::contains("K_")) %>%
    reshape2::melt(id=c("N","G","seed"), variable.name="type") %>%
    tidyr::separate(col="type", into=c("else","what","type")) %>% dplyr::mutate("else"=NULL) %>%
    tidyr::pivot_wider(names_from="what", values_from="value") %>%
    dplyr::mutate(value=found/true) %>%
    ggplot() +
    geom_violin(aes(x=as.factor(N), y=value)) +
    ggh4x::facet_nested(type ~ G) +
    labs(title="K_ratio") +
    theme_bw()

  return(
    patchwork::wrap_plots(mse_counts, cosine_expos, cosine_sigs, k_ratio, ncol=2)
  )
}



stats_single_data = function(fname) {
  cat(paste0(fname, "\n"))
  simul_fit = readRDS(fname)
  x.simul = simul_fit$dataset
  x.fit.0 = simul_fit$fit.0
  x.fit.N = simul_fit$fit.N

  idd = strsplit(fname, "/")[[1]]; idd = idd[[length(idd)]]

  stats_fit0 = stats_fitN = data.frame()

  try(expr = {stats_fit0 = eval_single_fit_matched(x.fit.0, x.simul) %>%
    dplyr::bind_cols() %>% dplyr::mutate(penalty="NoPenalty")})

  try(expr = {stats_fitN = eval_single_fit_matched(x.fit.N, x.simul) %>%
    dplyr::bind_cols() %>% dplyr::mutate(penalty="PenaltyN")})

  return(
    tibble::tibble(fname=fname,
                   "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][2] %>%
                     as.numeric(),
                   "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][3] %>%
                     as.numeric(),
                   "seed"=(stringr::str_replace_all(idd, "s", "") %>% strsplit(split="[.]"))[[1]][4] %>%
                     as.numeric(),
                   "idd"=idd) %>%
      dplyr::select(N, G, seed, idd, dplyr::everything()) %>%
      dplyr::bind_cols(rbind(stats_fit0, stats_fitN))
  )
}

