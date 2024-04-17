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


stats_fit_quality = function(x.fit, x.simul, suffix_name="") {
  rare_common = rare_common_sigs(x.simul)

  sigs.fit = get_signatures(x.fit); sigs.simul = get_signatures(x.simul)

  expos.fit = get_exposure(x.fit); expos.simul = get_exposure(x.simul)

  assigned_missing = get_assigned_missing(x.fit=x.fit, x.simul=x.simul, cutoff=cutoff)
  assigned = assigned_missing$assigned_tp
  unassigned = c(assigned_missing$missing_fn, assigned_missing$added_fp)

  assigned_missing_all = get_assigned_missing(x.fit=x.fit, x.simul=x.simul,
                                              cutoff=0.)

  mse_counts = compute.mse(m_true=x.simul$input$counts, m_inf=get_data(x.fit, reconstructed=T))
  mse_expos = compute.mse(m_true=expos.simul, m_inf=expos.fit,
                          assigned_missing=assigned_missing)

  cosine_sigs = compute.cosine(sigs.fit, sigs.simul,
                               assigned_missing=assigned_missing,
                               what="sigs")
  cosine_expos = compute.cosine(expos.fit, expos.simul,
                                assigned_missing=assigned_missing,
                                what="expos")

  ari_nmi = ari_nmi_km = ari_nmi_km_em = list(NA, NA)
  if (have_groups(x.fit)) {
    kmeans1 = get_obj_initial_params(x.fit)
    kmeans2 = get_obj_initial_params(x.fit) %>% run_kmeans()

    ari_nmi = compute_ari_nmi(x.simul=x.simul, x.fit=x.fit)
    ari_nmi_km1 = compute_ari_nmi(x.simul=x.simul, x.fit=kmeans1)
    ari_nmi_km2 = compute_ari_nmi(x.simul=x.simul, x.fit=kmeans2)
  }

  res = tibble::tibble(
    "assigned_missing"=list(assigned_missing),

    "sigs_found"=length(assigned),
    "shared_found"=sum(assigned %in% rare_common$shared),
    "private_found"=sum(assigned %in% c(rare_common$private_rare, rare_common$private_common)),
    "K_found"=length(get_signames(x.fit)),

    "mse_counts"=mse_counts,
    "mse_expos"=mse_expos,

    "cosine_sigs"=cosine_sigs,
    "cosine_expos"=cosine_expos,

    "groups_found"=length(x.fit$groups %>% unique()),
    "ari"=ari_nmi[[1]],
    "nmi"=ari_nmi[[2]],
    "ari_km1"=ari_nmi_km1[[1]],
    "nmi_km1"=ari_nmi_km1[[2]],
    "ari_km2"=ari_nmi_km2[[1]],
    "nmi_km2"=ari_nmi_km2[[2]]
  )

  colnames(res) = paste(colnames(res), suffix_name, sep="_")
  return(res)
}


compare_single_fit = function(fitname, fits_path, data_path, fits_pattern,
                              data_pattern="simul.",
                              cutoff=0.8, filtered_catalogue=TRUE,
                              min_exposure=0.,
                              save_plots=FALSE, check_plots=TRUE) {
  idd = stringr::str_replace_all(fitname, pattern=paste0(fits_pattern,"|.Rds"), replacement="")
  simulname = paste0(data_pattern, idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()

  x.fit.nolc = readRDS(paste0(fits_path, fitname)) %>%
    convert_sigs_names(x.simul, cutoff=cutoff) %>%
    # fix_assignments()
    recompute_centroids() %>% merge_clusters()

  x.fit.noadj = readRDS(paste0(fits_path, fitname)) %>%
    get_adjusted_fit_lc() %>%
    convert_sigs_names(x.simul, cutoff=cutoff)

  x.fit = x.fit.noadj %>%
    # fix_assignments()
    recompute_centroids() %>% merge_clusters()

  rare_common = rare_common_sigs(x.simul)

  if (save_plots &&
      !(check_plots && paste0("plots.", idd, ".pdf") %in% list.files(fits_path))) {

    all_sigs = c(get_signames(x.fit), get_signames(x.fit.nolc), get_signames(x.simul)) %>% unique()
    colors_ref = COSMIC_color_palette(catalogue=COSMIC_filt)[all_sigs] %>% purrr::discard(is.na)
    colors_dn = gen_palette(n=length(setdiff(all_sigs, names(colors_ref)))) %>%
      setNames(setdiff(all_sigs, names(colors_ref)))
    cls = c(colors_ref, colors_dn)

    pp = make_plots_compare(fit1=x.fit, fit2=x.simul,
                            name1="fit linear comb", name2="simul",
                            min_exposure=min_exposure, cls=cls)

    pdf(paste0(fits_path, "plots.", idd, ".pdf"), height=8, width=14)
    plot_fit(x.fit, x.simul, cls=cls, name1="Fit LC", name2="Simulated") %>% print()
    plot_fit(x.fit.nolc, x.simul, cls=cls, name1="Fit no LC", name2="Simulated") %>% print()
    plot_fit(x.fit, x.fit.noadj, cls=cls, name1="Fit LC", name2="Fit LC no EM") %>% print()
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
      "idd"=idd,
      "unique_id"=unique_id,

      "shared_true"=length(rare_common$shared),
      "common_true"=length(rare_common$private_common),
      "rare_true"=length(rare_common$private_rare),
      "K_true"=length(get_signames(x.simul)),

      "shared"=list(rare_common$shared),
      "priv_common"=list(rare_common$private_common),
      "priv_rare"=list(rare_common$private_rare)
    ) %>%
      tibble::add_column(stats_fit_quality(x.fit.nolc, x.simul, suffix_name="noLC")) %>%
      tibble::add_column(stats_fit_quality(x.fit, x.simul, suffix_name="LC")) %>%
      tibble::add_column(stats_fit_quality(x.fit.noadj, x.simul, suffix_name="LCnomerge"))
  )
}


make_plot_single_fit = function(fit1, titlee="fit1", min_exposure=0.,
                                cls=NULL, variant_type="SBS",
                                ref_catalogue=COSMIC_filt) {
  if (is.null(cls)) cls = get_color_palette(fit1)

  plot_umap = plot_umap_output(umap::umap(get_exposure(fit1)),
                               groups=get_groups(fit1))

  plot_counts = plot_mutations(fit1, reconstructed=T) + labs(title=paste0("Counts", titlee))

  plot_expos_centr = plot_exposures(fit1 %>% filter_exposures(min_expos=min_exposure),
                                    add_centroid=T, cls=cls) +
    labs(title=paste0("Exposures and centroids", titlee))

  plot_expos_real = plot_exposures_real(fit1, groups_true=fit2$groups,
                                        titlee=paste0(titlee, " exposures"))

  plot_sigs = plot_signatures(fit1, cls=cls, what=variant_type)

  plot_ref = plot_similarity_reference(fit1, reference=ref_catalogue, context=FALSE)

  return(list("expos_centr"=plot_expos_centr,
              "expos_real"=plot_expos_real,
              "signatures"=plot_sigs,
              "counts"=plot_counts,
              "umap"=plot_umap))
}




make_plots_compare = function(fit1, fit2, name1="fit1", name2="fit2",
                              min_exposure=0., cls=NULL,
                              variant_type="SBS", ref_catalogue=COSMIC_filt) {
  if (is.null(cls)) {
    cls = merge_colors_palette(fit1, fit2, ref_catalogue)
  }

  ttitle = paste0(" ", name1, " (top) and ", name2, " (bottom)")

  plot_umap = plot_umap_output(umap::umap(get_exposure(fit1)),
                               groups=get_groups(fit1)) %>%
    patchwork::wrap_plots(
      plot_umap_output(umap::umap(get_exposure(fit2)),
                       groups=get_groups(fit2))
    ) & theme(legend.position="bottom")

  plot_counts = plot_mutations(fit1, reconstructed=T, what=variant_type) %>%
    patchwork::wrap_plots(plot_mutations(fit2, reconstructed=F, what=variant_type), ncol=1) &
    patchwork::plot_annotation(title=paste0("Counts", ttitle))

  plot_expos_centr = plot_exposures(fit1 %>% filter_exposures(min_expos=min_exposure),
                              add_centroid=T, cls=cls) %>%
    patchwork::wrap_plots(plot_exposures(fit2, add_centroid=T, cls=cls),
                          ncol=1, guides="collect") &
    patchwork::plot_annotation(title=paste0("Exposures and centroids", ttitle))

  plot_expos_real = plot_exposures_real(fit1, groups_true=fit2$groups,
                      titlee=paste0(name1, " exposures"), cls=cls) %>%
    patchwork::wrap_plots(plot_exposures_real(fit2, groups_true=fit1$groups,
                                              titlee=paste0(name2, " exposures"), cls=cls),
                          guides="collect", ncol=1)

  plot_sigs = plot_signatures(fit1, cls=cls, what=variant_type) %>%
    patchwork::wrap_plots(plot_signatures(fit2, cls=cls, what=variant_type)) &
    patchwork::plot_annotation(title=paste0("Signatures", ttitle))

  return(list("expos_centr"=plot_expos_centr,
              "expos_real"=plot_expos_real,
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




# get_new_best = function(x.fit, score_name="bic") {
#   best = recompute_bic(x.fit, score_name=score_name) %>%
#     dplyr::filter(.data[[score_name]]==min(.data[[score_name]]))
#
#   new_best = get_fit_by_id(x=x.fit, idd=paste(best$K, best$groups, sep="."))
#   return(new_best)
# }


# recompute_bic = function(x.fit, score_name="bic") {
#   new_scores = get_K_scores(x.fit) %>%
#     dplyr::filter(score_id %in% c(score_name, "reg_llik")) %>%
#     tidyr::pivot_wider(names_from="score_id", values_from="score") %>%
#
#     dplyr::group_by(K, groups) %>%
#     dplyr::slice(which.min(.data[[score_name]])) %>%
#
#     dplyr::rowwise() %>%
#     dplyr::mutate(
#       new_score=compute_score(x.fit=get_fit_by_id(x.fit,
#                                                   idd=paste(K, groups, sep=".")),
#                               llik=reg_llik,
#                               score_name=score_name)) %>%
#     dplyr::ungroup()
#
#   return(new_scores)
# }


# compute_score = function(x.fit, llik, score_name="bic") {
#   n_pars = compute_n_pars(x.fit)
#
#   # k * torch.log(torch.tensor(n, dtype=torch.float64)) - (2 * _log_like)
#   if (score_name=="bic") return(n_pars * log(x.fit$n_samples) - 2*llik)
# }


compute_n_pars = function(x.fit) {
  return(
    prod(dim(get_denovo_signatures(x.fit))) +
      prod(dim(get_exposure(x.fit))) +
      length(unique(get_groups(x.fit))) +
      length(unique(get_groups(x.fit))) * length(get_signames(x.fit))
  )
}



