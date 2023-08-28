plot_sigs_clusters_found = function(stats_df, what=c("private", "shared", "all", "similar", "clusters"),
                                    ratio=F, fill="inf_type", facet=FALSE, scales="fixed", ylim=NULL,
                                    cols=NULL, suffix_name="") {
  if (suffix_name != "" && !grepl("^_", suffix_name)) suffix_name = paste0("_", suffix_name)
  stats_df = stats_df %>% dplyr::select(-dplyr::contains("plot")) %>%
    dplyr::mutate(n_private=n_rare+n_common)

  columns = dplyr::case_when(
    what == "private" ~ c(paste0("private_found", suffix_name), "n_private"),
    what == "shared" ~ c(paste0("shared_found", suffix_name), "n_shared"),
    what == "all" ~ c(paste0("sigs_found", suffix_name), "n_sigs"),
    what == "clusters" ~ c(paste0("groups_found", suffix_name), "G")
  )


  if (is.null(cols))
    cols = gen_palette(length(stats_df[[fill]] %>% unique())) %>%
      setNames(stats_df[[fill]] %>% unique())

  if (ratio)
    p = stats_df %>%
      dplyr::mutate(ratio = .data[[columns[1]]] / .data[[columns[2]]]) %>%
      ggplot() +
      geom_hline(yintercept=1, linetype="dashed", linewidth=0.5, color="grey") +
      geom_violin(aes_string(x="as.factor(N)", y="ratio", fill=fill, color=fill),
                  alpha=0.5, position=position_dodge(width=0.7)) +
      geom_jitter(aes_string(x="as.factor(N)", y="ratio", color=fill),
                  position=position_jitterdodge(jitter.width=0.1,
                                                jitter.height=0.,
                                                dodge.width=0.7),
                  size=1, alpha=0.5) +
      scale_color_manual(values=cols, name="") +
      scale_fill_manual(values=cols, name="") +
      theme_bw() + xlab("# samples") + ylab("Ratio") +
      labs(title=paste0("Ratio between ", stringr::str_replace_all(columns[1], "_", " "),
                        " and ", stringr::str_replace_all(columns[2], "_", " ")))

  valss = unique(round( c(stats_df[[columns[1]]], stats_df[[columns[2]]]) ))
  breakss = seq(0, max(valss), by=2)
  if (!ratio)
    p = stats_df %>%
      ggplot() +
      geom_abline(linetype="dashed", linewidth=0.5, color="grey") +
      geom_jitter(aes_string(x=columns[1], y=columns[2], color=fill),
                  size=1, height=0, width=0.1, alpha=0.5) +
      scale_color_manual(values=cols, name="") +
      theme_bw() +
      labs(title=paste0(stringr::str_replace_all(columns[1], "_", " "), " (x) and ",
                        stringr::str_replace_all(columns[2], "_", " "), " (y)")) +
      xlab(paste0(axis_lab, " found")) + ylab(paste0(axis_lab, " true")) +
      scale_x_continuous(breaks=breakss, limits=c(0, max(valss)+1)) +
      scale_y_continuous(breaks=breakss, limits=c(0, max(valss)+1))

  if (facet && ratio)
    p = p + ggh4x::facet_nested(G ~ clust_type, scales=scales)
  else if (facet && !ratio)
    p = p + ggh4x::facet_nested(G ~ clust_type, scales=scales)

  if (!is.null(ylim)) p = p + ylim(ylim)

  return(p)
}



plot_mse_cosine = function(stats_df, colname, facet=T, fill="inf_type",
                           cols=NULL, scales="fixed", ylim=NULL) {
  stats_df = stats_df %>% dplyr::select(-dplyr::contains("plot"))

  metric = strsplit(colname, "_")[[1]][1]
  quantity = strsplit(colname, "_")[[1]][2]
  if (length(strsplit(colname, "_")[[1]]) == 3)
    type = strsplit(colname, "_")[[1]][3] else type = ""

  metric = dplyr::case_when(
    metric == "mse" ~ "MSE",
    metric == "cosine" ~ sub('^(\\w?)', '\\U\\1', metric, perl=T),
    .default = metric
  )

  quantity = dplyr::case_when(
    quantity == "expos" ~ "exposures",
    quantity == "sigs" ~ "signatures",
    is.na(quantity) ~ "",
    .default = quantity
  )

  if (is.null(cols))
    cols = gen_palette(length(stats_df[[fill]] %>% unique())) %>%
      setNames(stats_df[[fill]] %>% unique())

  p = stats_df %>%
    ggplot() +
    geom_boxplot(aes_string(x="as.factor(N)", y=colname, color=fill),
                fill="white",
                alpha=0.5, position=position_dodge(width=0.7),
                outlier.size=0.5, width=0.5) +

    scale_color_manual(values=cols, name="") +
    theme_bw() + xlab("# samples") + ylab(metric) +
    labs(title=paste(metric, type, quantity))

  if (facet)
    p = p + ggh4x::facet_nested(G ~ clust_type, scales=scales)

  if (!is.null(ylim)) p = p + ylim(ylim)

  return(p)
}


report_stats = function(stats_df, fname, save_path,
                        fill="fits_pattern", suffix_name="LC") {
  pdf(paste0(save_path, "stats_report.sim", fname, ".pdf"), height=6, width=14)

  (plot_sigs_clusters_found(stats_df, what="all", facet=T, ratio=T,
                           scales="free_y", ylim=c(0,NA), fill=fill,
                           suffix_name=suffix_name) %>%
    patchwork::wrap_plots(
      plot_sigs_clusters_found(stats_df, what="private", facet=T, ratio=T,
                               scales="free_y", ylim=c(0,NA), fill=fill, suffix_name=suffix_name),
      guides="collect") & theme(legend.position="bottom")) %>% print()

  (plot_mse_cosine(stats_df, paste0("mse_counts_", suffix_name), scales="free_y",
                   ylim=c(0,1), fill=fill) %>%
    patchwork::wrap_plots(
      plot_mse_cosine(stats_df, paste0("mse_expos_", suffix_name), scales="free_y",
                      ylim=c(0,1), fill=fill),
      guides="collect") & theme(legend.position="bottom")) %>% print()

  (plot_mse_cosine(stats_df, paste0("cosine_sigs_", suffix_name), scales="free_y",
                   ylim=c(0,1), fill=fill) %>%
    patchwork::wrap_plots(
      plot_mse_cosine(stats_df, paste0("cosine_expos_", suffix_name), scales="free_y",
                      ylim=c(0,1), fill=fill),
      guides="collect") & theme(legend.position="bottom")) %>% print()

  (plot_sigs_clusters_found(stats_df %>% dplyr::filter(clust_type!="flat"),
                           what="clusters", facet=T, ratio=T, scales="free_y",
                           ylim=c(0,NA), fill=fill, suffix_name=suffix_name) %>%
    patchwork::wrap_plots(
      plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("nmi_", suffix_name),
                      scales="free_y", ylim=c(0,1), fill=fill),
      plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("ari_", suffix_name),
                      scales="free_y", ylim=c(0,1), fill=fill),
      guides="collect") & theme(legend.position="bottom")) %>% print()

  (plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("nmi_", suffix_name),
                   scales="free_y", ylim=c(0,1), fill=fill) %>%
      patchwork::wrap_plots(
        plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("nmi_km_", suffix_name),
                        scales="free_y", ylim=c(0,1), fill=fill),
        plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("nmi_km_em_", suffix_name),
                        scales="free_y", ylim=c(0,1), fill=fill),
        guides="collect", nrow=1) & theme(legend.position="bottom")) %>% print()

  (plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("ari_", suffix_name),
                   scales="free_y", ylim=c(0,1), fill=fill) %>%
      patchwork::wrap_plots(
        plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("ari_km_", suffix_name),
                        scales="free_y", ylim=c(0,1), fill=fill),
        plot_mse_cosine(stats_df %>% dplyr::filter(clust_type!="flat"), paste0("ari_km_em_", suffix_name),
                        scales="free_y", ylim=c(0,1), fill=fill),
        guides="collect", nrow=1) & theme(legend.position="bottom")) %>% print()

  dev.off()
}




# make_plots = function(fname, data_path, fits_path, save_path=NULL,
#                       title_id="", idd_name="") {
#   simul1 = readRDS(paste0(data_path, "simul.", fname, ".Rds")) %>%
#     create_basilica_obj_simul()
#
#   fit1 = readRDS(paste0(fits_path, "fit.", fname, idd_name, ".Rds")) %>%
#     convert_sigs_names(x.simul=simul1, cutoff=0.8) %>%
#     add_groups(simul1$groups)
#
#   hier1 = readRDS(paste0(fits_path, "fit.hier.", fname, idd_name, ".Rds")) %>%
#     convert_sigs_names(x.simul=simul1, cutoff=0.8)
#
#   assigned_missing_fit1 = get_assigned_missing(x.fit=fit1, x.simul=simul1, cutoff=0.8)
#   assigned_missing_hier1 = get_assigned_missing(x.fit=hier1, x.simul=simul1, cutoff=0.8)
#
#   all_sigs = c(get_signames(simul1), get_signames(fit1), get_signames(hier1)) %>% unique()
#
#   rare_common = rare_common_sigs(simul1)
#
#   colpal = gen_palette(length(all_sigs)) %>% setNames(all_sigs)
#
#   title_simul = paste0("SIMUL: Rare sigs ", toString(rare_common$private_rare),
#                        " - Common sigs ", toString(rare_common$private_common),
#                        " ", title_id)
#
#   title_fit = paste0("FIT: Rare sigs ", toString(rare_common$private_rare),
#                      " - Common sigs ", toString(rare_common$private_common),
#                      " ", title_id)
#   subtitle_fit = paste0("Missing sigs ", toString(assigned_missing_fit1$missing_fn),
#                         " - Added sigs ", toString(assigned_missing_fit1$added_fp))
#
#   title_hier = paste0("HIER FIT: Rare sigs ", toString(rare_common$private_rare),
#                       " - Common sigs ", toString(rare_common$private_common),
#                       " ", title_id)
#   subtitle_hier = paste0("Missing sigs ", toString(assigned_missing_hier1$missing_fn),
#                          " - Added sigs ", toString(assigned_missing_hier1$added_fp))
#
#   pl_simul = plot_fit(simul1, cls=colpal) +
#     patchwork::plot_annotation(title=title_simul)
#
#   pl_fit = plot_fit(fit1, simul, cls=colpal, reconstructed = T) +
#     patchwork::plot_annotation(title=title_fit,
#                                subtitle=subtitle_fit)
#
#   pl_hier = plot_fit(hier1, simul, cls=colpal, reconstructed = T) +
#     patchwork::plot_annotation(title=title_hier, subtitle=subtitle_hier)
#
#   pl_simil_fit1 = plot_similarity_reference(fit1, reference=get_signatures(simul1)) +
#     patchwork::plot_annotation(title=title_fit, subtitle=subtitle_fit)
#
#   pl_simil_hier1 = plot_similarity_reference(hier1, reference=get_signatures(simul1)) +
#     patchwork::plot_annotation(title=title_hier, subtitle=subtitle_hier)
#
#   if (!is.null(save_path)) {
#     pdf(paste0(save_path, "plots_", fname, idd_name, ".pdf"), height = 10, width = 14)
#     print(pl_simul)
#     print(pl_fit)
#     print(pl_simil_fit1)
#     print(pl_hier)
#     print(pl_simil_hier1)
#     dev.off()
#   }
#
#   return(
#     list("simul"=list(pl_simul),
#          "fit"=list("fit"=pl_fit, "simil"=pl_simil_fit1),
#          "hier"=list("fit"=pl_hier, "simil"=pl_simil_hier1))
#   )
#
# }
#
# plot_sigs_stats = function(stats_df, wrap=T, ratio=F, facet_groups=T, compare_runs=T) {
#   p1 = plot_sigs_found(stats_df, "rare", ratio, facet_groups, compare_runs)
#   p2 = plot_sigs_found(stats_df, "common", ratio, facet_groups, compare_runs)
#   p3 = plot_sigs_found(stats_df, "shared", ratio, facet_groups, compare_runs)
#   p4 = plot_sigs_found(stats_df, "all", ratio, facet_groups, compare_runs)
#
#   if (wrap)
#     return(patchwork::wrap_plots(p1, p2, p3, p4, ncol=2, guides="collect") &
#              theme(legend.position = "bottom"))
#
#   return(list(p1, p2, p3, p4))
# }
#
# plot_metrics_stats = function(stats_df, wrap=T, facet_groups=T) {
#   p1 = plot_mse_cosine(stats_df, "mse_counts", facet_groups)
#   p2 = plot_mse_cosine(stats_df, "mse_expos", facet_groups)
#   p3 = plot_mse_cosine(stats_df, "cosine_sigs", facet_groups)
#   p4 = plot_mse_cosine(stats_df, "cosine_expos", facet_groups)
#   p5 = plot_mse_cosine(stats_df, colname="mse_expos_rare", facet_groups)
#
#   if (wrap)
#     return(patchwork::wrap_plots(p1, p2, p3, p4, p5, ncol=2, guides="collect") &
#              theme(legend.position = "bottom"))
#
#   return(list(p1, p2, p3, p4))
# }
#
# stats_plots = function(stats_df, save_path, out_name) {
#   stats_plots_utils(stats_df, save_path, out_name, T, T)
#   stats_plots_utils(stats_df, save_path, paste0("nh.", out_name), T, F)
#   stats_plots_utils(stats_df, save_path, paste0("nh.ng.", out_name), F, F)
# }
#
# stats_plots_utils = function(stats_df, save_path, out_name,
#                              facet_groups=T, compare_runs=T) {
#   pdf(paste0(save_path, "plots.", out_name, ".pdf"), height=8, width=14)
#   plot_sigs_stats(stats_df, wrap=T, facet_groups=facet_groups,
#                   compare_runs=compare_runs) %>% print()
#   plot_sigs_stats(stats_df, wrap=T, ratio=T, facet_groups=facet_groups,
#                   compare_runs=compare_runs) %>% print()
#   plot_metrics_stats(stats_df, wrap=T, facet_groups=facet_groups,
#                      compare_runs=compare_runs) %>% print()
#   dev.off()
# }



plot_umap_output = function(umap_obj, groups) {
  umap_obj$layout %>% as.data.frame() %>%
    tibble::rownames_to_column() %>%
    dplyr::mutate(groupid=groups) %>%
    ggplot() +
    geom_point(aes(x=V1, y=V2, color=as.factor(groupid))) +
    theme_bw()
}

