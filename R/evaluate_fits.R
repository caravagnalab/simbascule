get_stats_df = function(data_path, fits_path, cutoff=0.8, min_exposure=0.,
                        data_pattern=c("simul."), run_id=c(""),
                        fits_pattern=c("fit.","fit_hier.","fit_clust."),
                        save_plots=FALSE, check_plots=TRUE) {
  if (is.null(names(run_id))) names(run_id) = fits_path
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
                             save_plots=save_plots, check_plots=check_plots) %>%

            dplyr::mutate(fits_pattern=stringr::str_replace_all(pattern_i, "\\\\.",""),
                          unique_id=paste0(fits_pattern, ".", unique_id),
                          fits_path=path_i,
                          data_path=data_path)

        } ) %>% do.call(what=rbind, args=.)
      } ) %>% do.call(what=rbind, args=.)
    }) %>% do.call(what=rbind, args=.) %>%
      dplyr::mutate(run_id=run_id[fits_path],
                    unique_id=paste(fits_pattern, run_id, sep="."))
  )
}


compare_single_fit = function(fitname, fits_path, data_path, fits_pattern,
                              data_pattern="simul.",
                              cutoff=0.8, filtered_catalogue=TRUE,
                              min_exposure=0.,
                              save_plots=FALSE, check_plots=TRUE) {
  idd = stringr::str_replace_all(fitname, pattern=paste0(fits_pattern,"|.Rds"), replacement="")
  simulname = paste0(data_pattern, idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()

  fit.nolc = readRDS(paste0(fits_path, fitname)) %>%
    convert_sigs_names(x.simul, cutoff=cutoff) %>%
    fix_assignments()

  fit.init = fit.nolc %>% get_adjusted_fit_lc() %>%
    convert_sigs_names(x.simul, cutoff=cutoff)

  x.fit = fit.init %>% fix_assignments()

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
                               assigned_missing=assigned_missing,
                               what="sigs")

  cosine_expos = compute.cosine(expos.fit, expos.simul,
                                assigned_missing=assigned_missing,
                                what="expos")
  cosine_expos_rare = compute.cosine(expos.fit, expos.simul,
                                     assigned_missing=assigned_missing,
                                     what="expos",
                                     subset_cols=c(rare_common$private_rare,
                                                   rare_common$private_common))

  if (have_groups(x.fit)) {
    ari_nmi = compute_ari_nmi(x.simul=x.simul, x.fit=x.fit)
    ari_nmi_km = compute_ari_nmi(x.simul=x.simul, x.fit=get_obj_initial_params(x.fit))
    ari_nmi_km_em = compute_ari_nmi(x.simul=x.simul, x.fit=fix_assignments(get_obj_initial_params(x.fit)))

  } else {
    ari_nmi = ari_nmi_km = ari_nmi_km_em = list(NA, NA)
  }

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
    coss = lsa::cosine(t( get_centroids(x.fit)[as.character(unique(x.fit$groups)),] ))
    mean_centr_simil = coss[upper.tri(coss)] %>% mean
  } else { mean_centr_simil = 0 }

  if (save_plots &&
      !(check_plots && paste0("plots.", idd, ".pdf") %in% list.files(fits_path))) {

    all_sigs = c(get_signames(x.fit), get_signames(fit.nolc), get_signames(x.simul)) %>% unique()
    colors_ref = COSMIC_color_palette(catalogue=COSMIC_filt)[all_sigs] %>%
      purrr::discard(is.na)
    colors_dn = gen_palette(n=length(setdiff(all_sigs, names(colors_ref)))) %>%
      setNames(setdiff(all_sigs, names(colors_ref)))
    cls = c(colors_ref, colors_dn)

    pp = make_plots_compare(fit1=x.fit, fit2=x.simul, name1="fit linear comb", name2="simul",
                            min_exposure=min_exposure, cls=cls)
    pdf(paste0(fits_path, "plots.", idd, ".pdf"), height=8, width=14)

    plot_fit(x.fit, x.simul, cls=cls, title="Fit with lcomb vs Simulated") %>% print()
    plot_fit(fit.nolc, x.simul, cls=cls, title="Fit no lcomb vs Simulated") %>% print()
    plot_fit(x.fit, fit.init, cls=cls, title="Fit with lcomb vs Fit initial") %>% print()
    patchwork::wrap_plots(pp$umap,
                          plot_scores(x.fit),
                          plot_gradient_norms(x.fit), ncol=2) %>% print()
    dev.off()
  }

  unique_id = strsplit(fits_path, "/")[[1]] %>% tail(1)

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
      "n_private_found"=n_private_found,

      "assigned_missing"=list(assigned_missing),

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
      "ari"=ari_nmi[[1]],
      "nmi"=ari_nmi[[2]],
      "ari_km"=ari_nmi_km[[1]],
      "nmi_km"=ari_nmi_km[[2]],
      "ari_km_em"=ari_nmi_km_em[[1]],
      "nmi_km_em"=ari_nmi_km_em[[2]],

      "shared"=list(rare_common$shared),
      "priv_common"=list(rare_common$private_common),
      "priv_rare"=list(rare_common$private_rare),

      "idd"=idd,
      "unique_id"=unique_id
    )
  )
}



make_plots_compare = function(fit1, fit2, name1="fit1", name2="fit2",
                              min_exposure=0., cls=NULL) {
  if (is.null(cls)) {
    all_sigs = unique(c(get_signames(fit1), get_signames(fit2)))
    cls = gen_palette(n=length(all_sigs)) %>% setNames(all_sigs)
  }

  ttitle = paste0(" ", name1, " (top) and ", name2, " (bottom)")

  plot_umap = plot_umap_output(umap::umap(get_exposure(fit1)),
                               groups=get_groups(fit1)) %>%
    patchwork::wrap_plots(
      plot_umap_output(umap::umap(get_exposure(fit2)),
                       groups=get_groups(fit2))
    ) & theme(legend.position="bottom")

  plot_counts = plot_mutations(fit1, reconstructed=T) %>%
    patchwork::wrap_plots(plot_mutations(fit2, reconstructed=F), ncol=1) &
    patchwork::plot_annotation(title=paste0("Counts", ttitle))

  plot_expos = plot_exposures(fit1 %>% filter_exposures(min_expos=min_exposure),
                              add_centroid=TRUE, cls=cls) %>%
    patchwork::wrap_plots(plot_exposures(fit2, add_centroid=TRUE, cls=cls),
                          ncol=1, guides="collect") &
    patchwork::plot_annotation(title=paste0("Exposures", ttitle))

  plot_centroids = plot_exposures(fit1 %>% filter_exposures(min_expos=min_exposure),
                                  centroids=TRUE, cls=cls) %>%
    patchwork::wrap_plots(plot_exposures(fit2, centroids=T, cls=cls),
                          ncol=1, guides="collect") &
    patchwork::plot_annotation(title=paste0("Centroids", ttitle))

  plot_sigs = plot_signatures(fit1, catalogue=get_signatures(fit2), cls=cls)

  plot_expos_centr = patchwork::wrap_plots(plot_expos + theme(legend.position="none"),
                                           plot_centroids,
                                           widths=c(3,1), guides="collect") &
    theme(legend.position="bottom") &
    patchwork::plot_annotation(title=paste0("Exposures and centroids", ttitle))

  return(list("exposures"=plot_expos,
              "centroids"=plot_centroids,
              "expos_centr"=plot_expos_centr,
              "signatures"=plot_sigs,
              "counts"=plot_counts,
              "umap"=plot_umap))
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



