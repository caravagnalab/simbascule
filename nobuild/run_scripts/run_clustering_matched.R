args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
run_id = args[2]

cat(paste("i =", i, "\n"))

main_path = "~/GitHub/"
fits_path = paste0("~/signatures/simulations/", "fits_dn.", run_id, "/")
save_path = paste0(fits_path, "clustering_test/")

cat(paste0("\nSaving in directory: ", save_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

source("~/GitHub/simbasilica/nobuild/run_scripts/fn_run.R")
reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

library(lsa)

cli::cli_process_done()



fits_names = list.files(fits_path, pattern=paste0("s", i, ".", run_id, ".Rds"), full.names=FALSE)
dir.create(save_path)


for (fname in fits_names) {
  G = strsplit(fname, "[.]")[[1]][3] %>% stringr::str_replace_all("G","") %>% as.integer()

  simul_fit = readRDS(paste0(fits_path, fname))
  x.nmf = simul_fit$fit.0

  x.cl.autoguide = fit_clustering(x.nmf, cluster=G*2,
                                  n_steps=3000, lr=0.005,
                                  nonparametric=TRUE,
                                  autoguide=TRUE,
                                  store_parameters=FALSE,
                                  store_fits=TRUE,
                                  seed_list=c(11,33,4392),
                                  py=py)

  x.cl.manguide = fit_clustering(x.nmf, cluster=G*2,
                                 n_steps=3000, lr=0.005,
                                 nonparametric=TRUE,
                                 autoguide=FALSE,
                                 store_parameters=FALSE,
                                 store_fits=TRUE,
                                 seed_list=c(11,33,4392),
                                 py=py)

  simul_fit$x.fit0.auto = x.cl.autoguide
  simul_fit$x.fit0.man = x.cl.manguide

  saveRDS(object=simul_fit, file=paste0(save_path, fname))
}




