# Sys.setenv(LD_LIBRARY_PATH=paste0("/u/elenab/.conda/envs/bascule-env/lib:", Sys.getenv("LD_LIBRARY_PATH")))
# export LD_LIBRARY_PATH=/u/elenab/.conda/envs/bascule-env/lib:$LD_LIBRARY_PATH

# library(signature.tools.lib)
devtools::load_all("~/GitHub/signature.tools.lib/")

#to read data
library(dplyr)
library(tidyr)
library(tibble)
library(Rcpp)
# library(readr)

args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i_tmp = as.integer(args[1])  # i from 1 to 540 - n of datasets to fit
to_add = as.integer(args[2])
i = i_tmp + to_add
# run_id = "SigFitTest"
run_id = args[3]  # SigFitTest or generative_model


cat(paste("i =", i, "\n"))

# Read input data
# dest_path = "~/share/signatures/"  # Demetra
# input_path = paste0(dest_path, "/data.", run_id, "/")

dest_path = paste0("~/scratch/signatures/fits_", run_id)  # Orfeo
input_path = paste0(dest_path, "/data.", run_id, "/")
output_path = paste0(dest_path, "/fits_sigtoolslib_E.", run_id, "/")
timing_file = paste0(output_path, "timing.txt")

if (!file.exists(output_path)) dir.create(output_path)

cosmic_path = list(
  "WGS"= paste0(dest_path, "/input/COSMIC_v3.3.1_SBS_GRCh38.txt"),
  "WES" = paste0(dest_path, "/input/COSMIC_v3_SBS_GRCh38-WES.txt")
)

cat(paste("\nSaving in directory:", output_path, "\n"))

cli::cli_process_start("Loading file and reference catalogs")

files = list.files(input_path, pattern="^simul_fit")
file_name = files[i]

if (grepl("WGS", file_name)) {
  catalogue_sbs = read.csv(cosmic_path[["WGS"]], sep="\t")
} else if (grepl("WES", file_name)) {
  catalogue_sbs = read.csv(cosmic_path[["WES"]], sep="\t")
} else {
  catalogue_sbs = read.csv(cosmic_path[["WGS"]], sep="\t")
}

mut_order = c('A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T')


# file name
file_path = paste0(input_path, file_name)

cat(paste0("\nFit for file: ", file_path, "\n"))

# SNV_catalogues = get_input(x$dataset, matrix=TRUE)$SBS %>% t() %>% as.data.frame()
# sigsToUse = getOrganSignatures(organ="Breast", version="latest")
# signatures = getCOSMICSignatures(version="latest")

x = readRDS(file_path)

data = x$dataset$input$SBS$counts %>%
  pivot_wider(names_from="features", values_from="value") %>%
  tibble::column_to_rownames(var="samples")

# rotate dataframe
catalogue_raw = as.data.frame(t(data))

# reorder sbs types
catalogue = catalogue_raw[match(mut_order, rownames(catalogue_raw)), ]

# use custom signatures
signatures = catalogue_sbs %>% tibble::column_to_rownames(var="Type")
signatures = signatures[mut_order, ]

if (file.exists(paste0(output_path, file_name))) {
  cli::cli_text("File already present. Not performing inference.")
} else {
  cli::cli_process_start("\nSignatureToolsLib fit\n")

  output_rds = paste0(output_path, file_name)
  output_dir = output_rds %>% stringr::str_replace_all(".Rds", "/")

  unlink(output_dir, recursive=TRUE)

  cat(paste0("\nSave extraction in: ", output_dir, "\n"))
  cat(paste0("\nSave RDS in: ", output_rds, "\n"))
  cat(paste0("\nSave time in: ", timing_file, "\n"))

  true_K = length(x$dataset$nmf$SBS$beta_denovo$sigs %>% unique())
  min_K = max(true_K - 2, 2)
  max_K = max(true_K + 2, min_K+1)

  nboots = 10 # CHANGE
  # n_cores = parallel::detectCores()

  start.time = Sys.time()
  # fit cosmic signatures
  subs_fit_ext = SignatureExtraction(cat = catalogue,
                                      matrix_of_fixed_signatures = signatures[, c("SBS1","SBS5")],
                                      nsig = c(min_K:max_K),
                                      nrepeats = 20, # CHANGE  # they suggest at least 100 repeats
                                      nboots = nboots,
				      parallel = TRUE,
                                      nparallel = 4,
                                      outFilePath = output_dir,
				      plotResultsFromAllClusteringMethods = FALSE
  )
  after.extraction = Sys.time()

  # find optimal number of signatures based on low "norm.Error orig. cat." and high "ave.SilWid.MC"

  normalize = function(x) {
    return((x - min(x)) / (max(x) - min(x)))
  }

  metrics_path = list.files(output_dir, pattern=glob2rx("^Sigs_OverallMetrics_extraction_nboots*.tsv$"))
  metrics_tmp = read.csv(file.path(output_dir, metrics_path), sep="\t") %>% as_tibble()

  metrics = metrics_tmp %>%
    mutate(NORM.ave.SilWid.MC=normalize(ave.SilWid.MC),
           NORM.error.orig=normalize(ave.KLD.orig/max(ave.KLD, ave.KLD.orig))) %>%
    select(nsig, starts_with("NORM"))

  best_K = metrics %>%
    mutate(ssum=NORM.ave.SilWid.MC + NORM.error.orig) %>%
    filter(ssum==max(ssum)) %>%
    pull(nsig)

  estimated_signatures = readTable(paste0(output_dir, "round_1/sig_", best_K,
                                          "/Sigs_plot_extraction_ns", best_K, "_nboots",
                                          nboots, ".tsv"))

  subs_fit_res = Fit(catalogues = catalogue,
                      signatures = estimated_signatures,
                      useBootstrap = TRUE,
                      nboot = nboots,
                      nparallel = 4,
                      randomSeed = 42,
                      exposureFilterType = "giniScaledThreshold",
                      verbose = TRUE
  )
  end.time = Sys.time()

  cli::cli_process_done()

  time_difference = end.time - start.time
  print(time_difference)

  # save results
  subs_fit_res$time = time_difference
  subs_fit_res$extraction_time = after.extraction - start.time
  saveRDS(subs_fit_res, output_rds)

  # write execution time
  line_to_write = paste0(file_name, ", ", time_difference)
  write(line_to_write, file = timing_file, append = TRUE)

}

