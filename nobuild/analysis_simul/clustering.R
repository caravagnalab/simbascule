devtools::load_all()
load_deps()
data_path ="~/GitHub/simbasilica/nobuild/simulations/synthetic_datasets_1606/"

simul = readRDS(list.files(data_path, full.names=T, pattern="N50")[1]) %>% create_basilica_obj_simul()

fit.cl = two_steps_inference(x=get_data(simul), regularizer = "cosine", reg_weight = 1,
                             k=0:3, clusters=2:3, seed_list=c(10),
                             reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5"),],
                             py=py, new_hier=TRUE)

fit.cl2 = two_steps_inference(x=get_data(simul),
                             k=0:10, clusters=c(2), seed_list=c(10),
                             reference_catalogue=COSMIC_filt_merged[c("SBS1","SBS5"),],
                             py=py, new_hier=FALSE)

fit.cl %>% convert_sigs_names(simul) %>%
  plot_fit(x.true=simul,
           cls=gen_palette(n=6) %>%
             setNames(get_signames(simul)),
           sample_name=T, sampleIDs=5)

