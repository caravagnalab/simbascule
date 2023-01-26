
#----------------------------------------------------------------------QC:PASSED
# split reference catalogue to 2 sub catalogue:
# reference catalogue (SBS1 included) + denovo catalogue
split.reference <- function(reference, ratio, seed=NULL) {

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  # read csv file as data.frame
  #reference <- read.table(reference_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
  num_ref <- round(ratio * nrow(reference))

  SBS1 <- reference['SBS1', ]   # save SBS1 (data.frame)
  reference <- reference[!(rownames(reference) %in% c("SBS1")), ] # excludes SBS1

  # shuffle the reference catalogue
  shuffled_reference = reference[sample(1:nrow(reference)), ]

  ref <- shuffled_reference[1:(num_ref-1), ]
  ref <- ref[order(rownames(ref)), ]
  ref <- rbind(SBS1, ref) # includes SBS1

  denovo <- shuffled_reference[num_ref:nrow(shuffled_reference), ]
  denovo <- denovo[order(rownames(denovo)), ]

  obj <- list(reference=ref, denovo=denovo)
  return(obj)
}
