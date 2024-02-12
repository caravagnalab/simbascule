args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])
run_id = args[2]

cat(paste("i =", i, "\n"))

main_path = "~/GitHub/"
fits_path = "~/fast/signatures/runtimes/"
data_path = "/orfeo/scratch/area/evillegas/mutational_signatures/simulations/fits_dn_matched_2011/"

cat(paste0("\nSaving in directory: ", fits_path, "\n\n"))

# Load packages #####

cli::cli_process_start("Loading packages")

source("~/GitHub/simbasilica/nobuild/run_scripts/fn_run.R")
reticulate::use_virtualenv("~/fast/virtualenv/basilica-env/", required=TRUE)
py = reticulate::import_from_path(module="pybasilica", path=paste0(main_path,"pybasilica/"))

devtools::load_all(paste0(main_path, "basilica"))
devtools::load_all(paste0(main_path, "simbasilica"))

library(lsa)

cli::cli_process_done()

N = c(150, 500, 1000)
G = c(1, 3, 6)
seed_list = 1:10

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


# Run model #####
tryCatch(
  {
    lapply(seed_list, function(s) {
      gen_run_aux(N=comb_i[1,"N_vals"],
                  G=comb_i[1,"G_vals"],
                  seed=s,
                  n_steps=3000,
                  private=private,
                  shared=shared,
                  reference_cat=catalogue_sbs,
                  run_fits=TRUE,
                  run_name=run_id,
                  path=fits_path,
                  data_path=data_path,
                  CUDA=TRUE)
    })
  }, error=function(e) print(reticulate::py_last_error())
)
