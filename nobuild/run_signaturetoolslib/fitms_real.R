library(signature.tools.lib)

#to read data
library(dplyr)
library(tidyr)
library(tibble)
library(readr)

#Read input data
input_path = "/orfeo/scratch/area/evillegas/mutational_signatures/fitms/real_data/input_data/"
output_path = "/orfeo/scratch/area/evillegas/mutational_signatures/fitms/real_data/fitms_fits/"
timing_file = paste0(output_path, "timing.txt")

files = list.files(path=input_path)

mut_order = c('A[C>A]A', 'A[C>A]C', 'A[C>A]G', 'A[C>A]T', 'C[C>A]A', 'C[C>A]C', 'C[C>A]G', 'C[C>A]T', 'G[C>A]A', 'G[C>A]C', 'G[C>A]G', 'G[C>A]T', 'T[C>A]A', 'T[C>A]C', 'T[C>A]G', 'T[C>A]T', 'A[C>G]A', 'A[C>G]C', 'A[C>G]G', 'A[C>G]T', 'C[C>G]A', 'C[C>G]C', 'C[C>G]G', 'C[C>G]T', 'G[C>G]A', 'G[C>G]C', 'G[C>G]G', 'G[C>G]T', 'T[C>G]A', 'T[C>G]C', 'T[C>G]G', 'T[C>G]T', 'A[C>T]A', 'A[C>T]C', 'A[C>T]G', 'A[C>T]T', 'C[C>T]A', 'C[C>T]C', 'C[C>T]G', 'C[C>T]T', 'G[C>T]A', 'G[C>T]C', 'G[C>T]G', 'G[C>T]T', 'T[C>T]A', 'T[C>T]C', 'T[C>T]G', 'T[C>T]T', 'A[T>A]A', 'A[T>A]C', 'A[T>A]G', 'A[T>A]T', 'C[T>A]A', 'C[T>A]C', 'C[T>A]G', 'C[T>A]T', 'G[T>A]A', 'G[T>A]C', 'G[T>A]G', 'G[T>A]T', 'T[T>A]A', 'T[T>A]C', 'T[T>A]G', 'T[T>A]T', 'A[T>C]A', 'A[T>C]C', 'A[T>C]G', 'A[T>C]T', 'C[T>C]A', 'C[T>C]C', 'C[T>C]G', 'C[T>C]T', 'G[T>C]A', 'G[T>C]C', 'G[T>C]G', 'G[T>C]T', 'T[T>C]A', 'T[T>C]C', 'T[T>C]G', 'T[T>C]T', 'A[T>G]A', 'A[T>G]C', 'A[T>G]G', 'A[T>G]T', 'C[T>G]A', 'C[T>G]C', 'C[T>G]G', 'C[T>G]T', 'G[T>G]A', 'G[T>G]C', 'G[T>G]G', 'G[T>G]T', 'T[T>G]A', 'T[T>G]C', 'T[T>G]G', 'T[T>G]T')

#Loop over all files
organs = c("Breast", "Colorectal", "Lung")
for (organ in organs){
    file_path = paste0(input_path, "sbs_", organ, ".csv")
    file_name = paste0(organ, ".Rds")
    print(file_path)

    #read data
    data <- read.csv(file_path, header = TRUE, check.names = FALSE)
    #set samples as row names
    row.names(data) <- data$samples
    data$samples <- NULL
    #rotate dataframe
    catalogue_raw = as.data.frame(t(data))
    #reorder sbs types
    catalogue <- catalogue_raw[match(mut_order, rownames(catalogue_raw)), ]

    
    start.time = Sys.time()
    #fit cosmic signatures
    subs_fit_res <- FitMS(catalogues = catalogue,
                      exposureFilterType = "giniScaledThreshold",
                      useBootstrap = FALSE,
                      organ = organ,
                      nboot=20,
                     )

    end.time = Sys.time()
    time_difference = end.time - start.time
    print(time_difference)

    #save results
    output_rds = paste0(output_path, file_name)
    saveRDS(subs_fit_res, output_rds)
    # plotFit(subs_fit_res, outdir = output_dir)
    
    #write execution time
    line_to_write = paste0(file_name, ", ", time_difference)
    write(line_to_write, file = timing_file, append = TRUE)
    
    }
