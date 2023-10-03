devtools::load_all()
load_deps()

main_path = "~/Dropbox/shared/2022. Basilica/real_data/"
mapping_table = readxl::read_xlsx(paste0(main_path, "processed_data/DBS_v1.01/DBS-doublet-base-substitution-classification.xlsx"),
                                  skip=2) %>% dplyr::select(-`Broad category`, -Count) %>%
  dplyr::rename(Reverse=`Reverse-complemented mutation`)

## map counts to cosmic colnames ####
counts = readRDS(paste0(main_path, "processed_data/catalogues/counts_dbs.Rds"))
signatures = COSMIC_dbs

curr_names = colnames(counts %>% dplyr::select(-organ, -cohort))
mapped_names = sapply(curr_names, function(i)
  if (i %in% colnames(signatures)) {
    return(i)
  } else if (i %in% mapping_table$Mutation) {
    return(mapping_table %>%
      dplyr::filter(Mutation==i) %>%
      dplyr::pull(Reverse))
  } else if (i %in% mapping_table$Reverse) {
    return(mapping_table %>%
      dplyr::filter(Reverse==i) %>%
      dplyr::pull(Mutation))
  }
  )

renaming = names(mapped_names) %>% setNames(mapped_names)
counts_renamed = counts %>% dplyr::rename(dplyr::all_of(renaming))

intersect(colnames(counts_renamed), colnames(signatures)) %>% unique() %>% length()

saveRDS(counts_renamed, paste0(main_path, "processed_data/catalogues/counts_dbs.Rds"))

