get_stats_df = function(data_path, fits_path, cutoff=0.8, min_exposure=0.,
                        data_pattern=c("simul."),
                        fits_pattern=c("fit.","fit_hier.","fit_clust."),
                        save_plots=FALSE) {

  fits_pattern = fits_pattern %>% stringr::str_replace_all("\\.","\\\\.")

  return(
    lapply(fits_path, function(path_i) {

      lapply(fits_pattern, function(pattern_i) {
        fits_i = list.files(path=path_i, pattern=pattern_i)

        lapply(fits_i, function(fitname) {
          cat(paste(path_i, fitname, pattern_i, "\n"))
          compare_single_fit(fitname=fitname, fits_path=path_i, data_path=data_path,
                             data_pattern=data_pattern, fits_pattern=pattern_i,
                             cutoff=cutoff, min_exposure=min_exposure,
                             save_plots=save_plots) %>%
            dplyr::mutate(fits_pattern=pattern_i,
                          fits_path=path_i,
                          data_path=data_path)
        } ) %>% do.call(what=rbind, args=.)
      } ) %>% do.call(what=rbind, args=.)
    }) %>% do.call(what=rbind, args=.) %>%
      dplyr::mutate(fits_pattern=stringr::str_replace_all(fits_pattern, "\\\\.",""))
  )
}


compare_single_fit = function(fitname, fits_path, data_path, cutoff=0.8,
                              fits_pattern, data_pattern="simul.",
                              filtered_catalogue=TRUE, min_exposure=0.,
                              save_plots=FALSE) {
  idd = stringr::str_replace_all(fitname, pattern=paste0(fits_pattern,"|.Rds"), replacement="")
  simulname = paste0(data_pattern, idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()
  x.fit = readRDS(paste0(fits_path, fitname)) %>%
    get_new_best(score_name="bic")

  if (x.fit$n_denovo > 0 && filtered_catalogue)
    x.fit$fit$denovo_signatures = renormalize_denovo_thr(x.fit$fit$denovo_signatures)

  x.fit = x.fit %>% convert_sigs_names(x.simul, cutoff=cutoff) # %>% merge_clusters()

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
                               subset_cols=c(rare_common$private_rare,
                                             rare_common$private_common))


  cosine_sigs = compute.cosine(sigs.fit, sigs.simul,
                               # assigned_missing=assigned_missing,
                               assigned_missing=assigned_missing,
                               what="sigs")

  cosine_expos = compute.cosine(expos.fit, expos.simul,
                                # assigned_missing=assigned_missing,
                                assigned_missing=assigned_missing,
                                what="expos")
  cosine_expos_rare = compute.cosine(expos.fit, expos.simul,
                                     # assigned_missing=assigned_missing,
                                     assigned_missing=assigned_missing,
                                     what="expos",
                                     subset_cols=c(rare_common$private_rare,
                                                   rare_common$private_common))

  # rare_freqs = length(get_samples_with_sigs(x.simul, sigs=rare_common$private_rare[1])) / x.simul$n_samples

  if (have_groups(x.fit)) {
    # groups_new = get_groups_rare(x.simul, rare_common)

    groups_fit = x.fit$groups; groups_simul = x.simul$groups
    if (length(unique(groups_fit)) == 1 || length(unique(groups_simul))==1) {
      groups_fit = c(groups_fit, "imolabella")
      groups_simul = c(groups_simul, "imolabella")
    }

    ari = aricode::ARI(x.simul$groups, x.fit$groups)
    nmi = aricode::NMI(groups_fit, groups_simul)

    # ari_rare = aricode::ARI(groups_new, x.fit$groups)
    # nmi_rare = aricode::NMI(groups_new, x.fit$groups)
  } else {
    ari = nmi = ari_rare = nmi_rare = NA
  }

  # n_rare_found = sum(assigned %in% rare_common$private_rare)
  # n_common_found = sum(assigned %in% rare_common$private_common)
  n_private_found = sum(assigned %in% c(rare_common$private_rare, rare_common$private_common))
  n_shared_found = sum(assigned %in% rare_common$shared)

  true_n_sigs = nrow(sigs.simul)
  inf_n_sigs = nrow(sigs.fit)

  similarity_fit = lsa::cosine(t(get_signatures(x.fit))) %>% as.data.frame() %>%
    tibble::rownames_to_column(var="sigs1") %>%
    reshape2::melt(variable.name="sigs2") %>% dplyr::mutate(sigs2=as.character(sigs2)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(combb=paste(sort(c(sigs1,sigs2)), collapse=",")) %>%
    dplyr::filter(sigs1 != sigs2, value > cutoff) %>%
    dplyr::select(-sigs1,-sigs2) %>% unique()
  n_sigs_similar = nrow(similarity_fit)

  if (length(unique(x.fit$groups)) > 1) {
    coss = lsa::cosine(t( x.fit$fit$params$alpha_prior[unique(x.fit$groups)+1,] ))
    mean_centr_simil = coss[upper.tri(coss)] %>% mean
  } else { mean_centr_simil = 0 }

  if (save_plots) {
    all_sigs = unique(c(get_signames(x.simul), get_signames(x.fit)))
    cls = gen_palette(n=length(all_sigs)) %>% setNames(all_sigs)

    plot_counts = plot_mutations(x.fit, reconstructed=T) %>%
      patchwork::wrap_plots(plot_mutations(x.simul, reconstructed=F), ncol=1) &
      patchwork::plot_annotation(title="Counts fit (top) and simulated (bottom)")

    plot_expos = plot_exposures(x.fit %>% filter_exposures(min_expos=min_exposure),
                                add_centroid=TRUE, cls=cls) %>%
      patchwork::wrap_plots(plot_exposures(x.simul, add_centroid=TRUE, cls=cls), ncol=1) &
      patchwork::plot_annotation(title="Exposures fit (top) and simulated (bottom)")

    plot_centroids = plot_exposures(x.fit %>% filter_exposures(min_expos=min_exposure),
                                    centroids=TRUE, cls=cls) %>%
      patchwork::wrap_plots(plot_exposures(x.simul, centroids=T, cls=cls), ncol=1) &
      patchwork::plot_annotation(title="Centroids fit (top) and simulated (bottom)")

    plot_sigs = plot_signatures(x.fit, catalogue=get_signatures(x.simul), cls=cls)

    pdf(paste0(fits_path, "plots.", idd, ".pdf"), height=8, width=14)
    patchwork::wrap_plots(plot_expos, plot_centroids, widths=c(3,1), guides="collect") %>% print()
    plot_counts %>% print()
    plot_sigs %>% print()
    dev.off()
  }

  return(
    tibble::tibble(
      "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][1] %>%
        as.numeric(),
      "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][2] %>%
        as.numeric(),

      "n_shared"=length(rare_common$shared),
      "n_common"=length(rare_common$private_common),
      "n_rare"=length(rare_common$private_rare),

      # "rare_freq"=rare_freqs,

      "n_shared_found"=n_shared_found,
      "n_private_found"=n_private_found,

      "assigned_missing"=list(assigned_missing),

      # "missing_fn"=list(assigned_missing$missing_fn),
      # "added_fp"=list(assigned_missing$added_fp),
      # "assigned"=list(assigned_missing$assigned_tp),

      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,
      "mse_expos_rare"=mse_expos_rare,

      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "cosine_expos_rare"=cosine_expos_rare,

      "n_sigs"=true_n_sigs,
      "n_sigs_found"=inf_n_sigs,

      "n_sigs_similar"=n_sigs_similar,
      "mean_centr_simil"=mean_centr_simil,

      "n_groups_found"=length(x.fit$groups %>% unique()),
      "ari"=ari,
      "nmi"=nmi,

      "shared"=list(rare_common$shared),
      "priv_common"=list(rare_common$private_common),
      "priv_rare"=list(rare_common$private_rare),

      "idd"=idd,
      # "plot_counts"=list(plot_counts),
      # "plot_expos"=list(plot_expos),
      # "plot_centroids"=list(plot_centroids),
      # "plot_sigs"=list(plot_sigs)
    )
  )
}





get_simul_fit = function(stats_df, condition, return_fit=T) {
  filtered = stats_df %>% dplyr::filter(eval(parse(text=condition))) %>%
    dplyr::arrange(dplyr::desc(nmi)) %>% dplyr::sample_n(1)

  if (!return_fit)
    return(filtered)

  fpath = filtered %>% dplyr::pull(fits_path)
  spath = filtered %>% dplyr::pull(data_path)
  x.fit = readRDS(paste0(fpath, filtered$fits_pattern, ".", filtered$idd, ".Rds"))
  x.simul = readRDS(paste0(data_path, "simul.", filtered$idd, ".Rds")) %>% create_basilica_obj_simul()

  return(
    list("filtered"=filtered,
         "x.fit"=x.fit,
         "x.simul"=x.simul,
         "fpath"=fpath,
         "idd"=filtered$idd)
  )
}




get_new_best = function(x.fit, score_name="bic") {
  best = recompute_bic(x.fit, score_name=score_name) %>%
    dplyr::filter(.data[[score_name]]==min(.data[[score_name]]))

  new_best = get_fit_by_id(x=x.fit, idd=paste(best$K, best$groups, sep="."))
  return(new_best)
}


recompute_bic = function(x.fit, score_name="bic") {
  new_scores = get_K_scores(x.fit) %>%
    dplyr::filter(score_id %in% c(score_name, "reg_llik")) %>%
    tidyr::pivot_wider(names_from="score_id", values_from="score") %>%

    dplyr::group_by(K, groups) %>%
    dplyr::slice(which.min(.data[[score_name]])) %>%

    dplyr::rowwise() %>%
    dplyr::mutate(
      new_score=compute_score(x.fit=get_fit_by_id(x.fit,
                                                  idd=paste(K, groups, sep=".")),
                              llik=reg_llik,
                              score_name=score_name)) %>%
    dplyr::ungroup()

  return(new_scores)
}


compute_score = function(x.fit, llik, score_name="bic") {
  n_pars = compute_n_pars(x.fit)

  # k * torch.log(torch.tensor(n, dtype=torch.float64)) - (2 * _log_like)
  if (score_name=="bic") return(n_pars * log(x.fit$n_samples) - 2*llik)
}


compute_n_pars = function(x.fit) {
  return(
    prod(dim(get_denovo_signatures(x.fit))) +
      prod(dim(get_exposure(x.fit))) +
      length(unique(get_groups(x.fit))) +
      length(unique(get_groups(x.fit))) * length(get_signames(x.fit))
  )
}



