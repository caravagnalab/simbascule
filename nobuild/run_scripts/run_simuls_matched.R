args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
run_id = args[2]

cat(paste("i =", i, "\n"))

main_path = "~/GitHub/"
fits_path = paste0("~/signatures/simulations/", "fits_dn.", run_id, "/")

cat(paste0("\nSaving in directory: ", fits_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

source("~/GitHub/simbasilica/nobuild/run_scripts/fn_run.R")
reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module = "pybasilica", path = paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

library(lsa)

cli::cli_process_done()

N = c(150, 500, 1000)
G = c(1, 3, 6)
seed_list = 1:30

shared = list("SBS"=c("SBS1","SBS5"),
              "DBS"=c("DBS3","DBS5"))
private = list("SBS"=paste("SBS", c(2,4,6,"7a","7c","10a",13,"17b",18,20,22,31,35,36), sep=""),
               "DBS"=setdiff(rownames(COSMIC_dbs), shared$DBS))
catalogue_sbs = list("SBS"=COSMIC_filt[c(shared[["SBS"]], private[["SBS"]]),],
                     "DBS"=COSMIC_dbs[c(shared[["DBS"]], private[["DBS"]]),])

set.seed(1234)
comb = expand.grid(N_vals=N, G_vals=G) %>%
  dplyr::arrange(N_vals)

comb_i = comb[i,]

# a = gen_run_aux(N=comb_i[1,"N_vals"],
#                 G=comb_i[1,"G_vals"],
#                 seed=1,
#                 private=private,
#                 shared=shared,
#                 n_steps=10,
#                 path=NULL)
# a$dataset %>% plot_fit()
# a$fit.0 %>% plot_fit()
# a$fit.N %>% plot_fit()


# Run model #####
lapply(seed_list, function(s) {
  gen_run_aux(N=comb_i[1,"N_vals"],
              G=comb_i[1,"G_vals"],
              seed=s,
              private=private,
              shared=shared,
              run_fits=TRUE,
              run_name=run_id,
              path=fits_path)
})


