args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i_tmp = as.integer(args[1])  # i from 1 to 540 - n of datasets to fit
to_add = as.integer(args[2])

i = i_tmp + to_add

# run_id = args[2]
run_id = "SigFitTest"

cat(paste("i =", i, "\n"))

main_path = "~/GitHub/"
dest_path = "~/scratch/signatures/SigFitTest/"
data_path = paste0(dest_path, "/data.", run_id, "/")
fits_path = paste0(dest_path, "/fits_dn.", run_id, "/")
cosmic_path = list(
  "WGS"= paste0(dest_path, "/input/COSMIC_v3.3.1_SBS_GRCh38.txt"),
  "WES" = paste0(dest_path, "/input/COSMIC_v3_SBS_GRCh38-WES.txt")
)

cli::cli_text("\nSaving in directory: {fits_path}")

# Load packages #####

cli::cli_process_start("Loading packages")

source("~/GitHub/simbascule/nobuild/run_scripts/fn_run_SigFitTest.R")
reticulate::use_condaenv("basilica-env")
py = reticulate::import_from_path(module="pybascule", path=paste0(main_path,"pybascule/"))

devtools::load_all(paste0(main_path, "bascule"))
devtools::load_all(paste0(main_path, "simbascule"))

library(lsa)

cli::cli_process_done()

# N = c(150, 500, 1000)
# G = c(1, 3, 6)
# seed_list = 1:30

cli::cli_process_start("Loading file and reference catalogs")

file_i = list.files(data_path, full.names=TRUE)[i]

shared = list("SBS"=c("SBS1","SBS5"))

if (grepl("WGS", file_i)) {
  catalogue_sbs = read.csv(cosmic_path[["WGS"]], sep="\t")
} else if (grepl("WES", file_i)) {
  catalogue_sbs = read.csv(cosmic_path[["WES"]], sep="\t")
}

reference_cat = list(
  "SBS"=t(tibble::column_to_rownames(catalogue_sbs, var="Type"))[shared[["SBS"]],] %>% as.data.frame()
)

cli::cli_process_done()


cli::cli_process_start("BASCULE fit of {file_i}\n")

# Run model #####
tryCatch(
  {
    gen_run_aux(fname=file_i,
                reference_cat=reference_cat,
                keep_sigs=shared[["SBS"]],
                n_steps=2000,
                run_fits=TRUE,
                save_path=fits_path,
                py=py)
  }, error=function(e) {
    print(e)
    print(reticulate::py_last_error())
    cli::cli_process_failed(msg_failed=e)
    }
  )

cli::cli_process_done()


