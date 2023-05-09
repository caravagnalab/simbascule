simul_path = "~/Desktop/accento_fs/signatures/simul_data_obj"

get_synthetic_datas = function(simul_path) {

  simul_files = list.files(path=simul_path, pattern="simul.", full.names=T)

  simul_df = data.frame() %>%
    dplyr::mutate(N=as.character(NA), G=as.character(NA), s=as.character(NA),
                  sigs_true=NA, exp_true=NA, counts_true=NA,
                  dn_inf=NA, fixed_inf=NA, exp_inf=NA)

  for (sf in simul_files) {
    x = readRDS(sf)
    splt = strsplit(sf, "/")[[1]]
    splt = strsplit(splt[length(splt)], "[.]")[[1]]

    N = splt[2]; G = splt[3]; s = splt[4]

    x.fit = readRDS(paste0(simul_path, "/fit.", N, ".", G, ".", s, ".Rds"))
    x.fit.hier = readRDS(paste0(simul_path, "/fit.hier.", N, ".", G, ".", s, ".Rds"))
  }

}
