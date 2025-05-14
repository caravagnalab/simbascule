library(signature.tools.lib)

#to read data
library(dplyr)
library(tidyr)
library(tibble)
library(readr)


args = commandArgs(trailingOnly = TRUE)
cat(paste("\nArguments:", paste(args, collapse=", "), "\n"))

i = as.integer(args[1])  # from 1 to 3 if run_id is GEL, from 1 to 7 if run_id is ICGC
run_id = args[2]  # GEL or TCGA


cat(paste("i =", i, "\n"))
cat(paste("Dataset ID =", run_id, "\n"))


# Read input data
main_path = "/orfeo/scratch/cdslab/ebusca00/signatures/fits_real_data/"
input_path = paste0(main_path, "real_data_", run_id, "/")
output_path = paste0(main_path, "fits_sigtoolslib.real_data_", run_id, "/")
# timing_file = paste0(output_path, "timing.txt")

files = list.files(path=input_path)

mut_order = c('A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T')

file_name = files[i]
file_path = file.path(input_path, file_name)

cat("\nFit for file: ", file_path, "\n")
cat("\nSave RDS in: ", output_path, "\n")
# cat("\nSave time in: ", timing_file, "\n")

# read data
if (grepl(".csv$", file_path)) {
  data = read.csv(file_path, header=TRUE, check.names=FALSE) %>%
    as_tibble() %>%
    column_to_rownames(var="samples")
} else if (grepl(".Rds$", file_path)) {
  data = readRDS(file_path)$input$SBS$counts %>%
    pivot_wider(names_from="features", values_from="value") %>%
    tibble::column_to_rownames(var="samples") %>% as.matrix()
}

# rotate dataframe
catalogue_raw = as.data.frame(t(data))
# reorder sbs types
catalogue = catalogue_raw[match(mut_order, rownames(catalogue_raw)), ]

organ_name = file_name %>% stringr::str_remove_all("surv.|sbs_|.csv|.Rds") %>% rapportools::tocamel(upper=TRUE)

catalogue_sbs = getOrganSignatures(organ=organ_name, version="latest", typemut="subs", verbose=TRUE)

# fit cosmic signatures
start.time = Sys.time()
subs_fit_res <- Fit(catalogues = catalogue,
                    signatures = catalogue_sbs,
                    useBootstrap = TRUE,
                    nboot = 20,
                    nparallel = 8,
                    randomSeed = 42,
                    exposureFilterType = "giniScaledThreshold",
                    verbose = TRUE
                    )

end.time = Sys.time()
time_difference = end.time - start.time
subs_fit_res$time = time_difference

print(time_difference)

#save results
output_rds = paste0(output_path, file_name) %>% stringr::str_replace_all(".csv",".Rds")
saveRDS(subs_fit_res, output_rds)

# # write execution time
# line_to_write = paste0(file_name, ", ", time_difference)
# write(line_to_write, file = timing_file, append = TRUE)
