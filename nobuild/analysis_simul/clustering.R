devtools::load_all()
load_deps()
data_path ="~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/"

simul = readRDS(list.files(data_path, full.names=T, pattern="N300")[1]) %>% create_basilica_obj_simul()
fit.h.old = readRDS("~/GitHub/simbasilica/nobuild/simulations/run_newmodel_2606_dn/fit.hier.N300.G2.s1.Rds")

fit.cl = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 1,
                             k=0:10, clusters=1:5, seed_list=c(10),
                             reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5"),],
                             py=py, new_hier=TRUE)

fit.h = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 1,
                            k=0:5, groups=simul$groups-1, # seed_list=c(10),
                            reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5"),],
                            py=py, new_hier=TRUE)

fit = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 1,
                          k=0:5, groups=NULL, # seed_list=c(10),
                          reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5"),],
                          py=py, new_hier=FALSE, regul_denovo=F)

fit.cl2 = two_steps_inference(x=get_data(simul),
                             k=0:10, clusters=c(2), seed_list=c(10),
                             reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5"),],
                             py=py, new_hier=FALSE)

fit.cl %>% convert_sigs_names(simul) %>%
  plot_fit(x.true=simul,
           cls=gen_palette(n=6) %>%
             setNames(get_signames(simul)),
           sample_name=T, sampleIDs=5)

