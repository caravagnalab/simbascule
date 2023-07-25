images_path = paste0("~/Google Drive/My Drive/work/signatures/bits/images/")

simul = readRDS("~/GitHub/simbasilica/nobuild/poster_bits/example_simul/simul.N100.G2.s1.raref0.05.Rds") %>%
  create_basilica_obj_simul()
fit = readRDS("~/GitHub/simbasilica/nobuild/poster_bits/example_simul/fit.N100.G2.s1.raref0.05.Rds")
  # convert_sigs_names(simul) %>% add_groups(simul$groups)
fit.hier = readRDS("~/GitHub/simbasilica/nobuild/poster_bits/example_simul/fit.hier.N100.G2.s1.raref0.05.Rds") %>%
  convert_sigs_names(simul)

sigs = unique(c(get_signames(fit %>% convert_sigs_names(simul)), get_signames(simul)))
col_palette = gen_palette(length(sigs)) %>% setNames(sigs)
assigned_dn = grepl("D", get_assigned_missing(fit, simul)$assigned_tp)
cols_dn = col_palette[names(get_assigned_missing(fit, simul)$assigned_tp)[assigned_dn]] %>%
  setNames(get_assigned_missing(fit, simul)$assigned_tp[assigned_dn])
col_palette = c(col_palette, cols_dn)

simul$groups = rep("Simulated", simul$n_samples)
fit$groups = rep("Predicted", simul$n_samples)
muts_simul = plot_mutations(simul, by_sig=T, cls=col_palette)
muts_fit = plot_mutations(fit, by_sig=T, cls=col_palette, epsilon=FALSE)
muts = patchwork::wrap_plots(muts_simul / muts_fit, guides="collect") &
  ylab("# mutations")
ggsave(paste0(images_path, "simul_muts.pdf"), height=6, width=12)


sigs_all = plot_signatures(fit, catalogue=get_signatures(simul)["SBS22",], cls=col_palette)
expos_simul = plot_exposures(simul, cls=col_palette)
expos_fit = plot_exposures(fit, cls=col_palette)
expos = patchwork::wrap_plots(expos_simul / expos_fit, guides="collect")
similarity = plot_similarity_reference(fit, reference=get_signatures(simul))

prevalence = get_exposure(simul, long=T) %>%
  dplyr::group_by(Signature) %>%
  dplyr::arrange(Signature, Exposure) %>%
  dplyr::mutate(nnn=1, cs=cumsum(nnn)/simul$n_samples) %>%
  ggplot(aes(x=Exposure, y=1-cs, color=Signature)) +
  geom_point(size=.5) +
  geom_line() +
  scale_color_manual(values=col_palette) +
  theme_bw() + ylab("Fraction of samples") +
  xlim(0,1) + ylim(0,1) + theme(legend.position = "bottom") + coord_fixed()

ggsave(paste0(images_path, "simul_expos.pdf"), plot=expos, height=6, width=12)
ggsave(paste0(images_path, "simul_similarity.pdf"), plot=similarity, height=7, width=12)
ggsave(paste0(images_path, "simul_sigs.pdf"), plot=sigs_all, height=7, width=12)
ggsave(paste0(images_path, "simul_preval.pdf"), plot=prevalence, height=6, width=4)



