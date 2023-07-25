base_path = "~/GitHub/simbasilica/nobuild/"
devtools::load_all()
load_deps()

## Fit on CRC and Lung ####
fname = "CRC_L_100"
fit.cl.h = readRDS(paste0(base_path, "poster_bits/fits_gel/fit_gel.", fname, ".h.Rds"))
fit.cl.nh = readRDS(paste0(base_path, "poster_bits/fits_gel/fit_gel.", fname, ".nh.Rds")) %>%
  add_groups(fit.cl.h$groups)

real_d = create_basilica_obj_real_data(base_path,
                                       counts=fit.cl.nh$input$counts,
                                       groupids=fit.cl.nh$groups)

real_d_f = filter_exposures(real_d, min_exp = 0.1)
# assigned_missing = get_assigned_missing(x.fit=fit.cl.nh, x.simul=real_d_f, cutoff=0.6)
fit.cl.nh = fit.cl.nh %>% convert_sigs_names(real_d, cutoff=0.6)
plot_fit(x=fit.cl.nh, x.true=real_d)

expos1 = plot_exposures(fit.cl.nh) + theme(legend.position = "bottom")
expos2 = plot_exposures(real_d) + theme(legend.position = "bottom")

patchwork::wrap_plots(expos1 / expos2)

plot_mutations(fit.cl.nh, by_sig = T)

a = plot_similarity_reference(fit.cl.nh, reference=get_signatures(real_d_f))
ggsave(paste0(base_path, "poster_bits/fits_gel/simil.", fname, ".nh.pdf"),
       height = 20, width = 10)




plot_similarity_reference(fit.cl.h, reference=get_signatures(fit.cl.h))
fit.cl.h %>% plot_signatures()



fit.cl.nh %>% filter_exposures(min_exp=.15) %>%
  plot_fit()





