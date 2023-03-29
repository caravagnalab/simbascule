#----------------------------------------------------------------------QC:PASSED
# generate signatures which includes:
#   * fixed signatures (SBS1 included)
#   * denovo signatures
# similarity between catalogue signatures are less than threshold
# similarity between denovo signatures are less than threshold
# but similarity between catalogue and denovo signatures are not taken to the account (future work)
generate.signatures <- function(
    reference_catalogue,
    denovo_catalogue,
    reference_cosine, # cosine similarity matrix of reference signatures (SBS1 excluded)
    denovo_cosine,    # cosine similarity matrix of denovo signatures
    complexity,       # c(3,2) | "low" | "medium" | "high" |
    similarity_limit,
    seed=NULL
) {

  if (complexity == -1)
    return(list(fixed = reference_catalogue, denovo = denovo_catalogue))

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  if (is.numeric(complexity) & length(complexity)==2) {
    fixed_num <- complexity[1]
    denovo_num <- complexity[2]
  } else if (is.character(complexity)) {
    if (complexity=='low') {
      fixed_num <- sample(3:5, 1)
      denovo_num <- sample(0:2, 1)
    }
    else if (complexity=='medium') {
      fixed_num <- sample(1:2, 1)
      denovo_num <- sample(3:5, 1)
    }
    else if (complexity=='high') {
      fixed_num <- sample(3:5, 1)
      denovo_num <- sample(3:5, 1)
    }
    else {
      stop("complexity argument should be selected from {'low', 'medium', 'high'}")
    }
  } else {
    stop("wrong complexity argument!")
  }

  SBS1 <- reference_catalogue['SBS1', ] # save SBS1 (data.frame)
  reference <- reference_catalogue[!(rownames(reference_catalogue) %in% c("SBS1")), ] # excludes SBS1

  # catalogue signatures -------------------------------------------------------

  if (fixed_num > 1) {

    while (TRUE) {
      shuffled_reference = reference[sample(1:nrow(reference)), ]
      signatures <- rownames(shuffled_reference[1:(fixed_num-1), ])
      cos_matrix <- reference_cosine[c("SBS1", signatures), c("SBS1", signatures)]
      for (i in 1:nrow(cos_matrix)) {
        cos_matrix[i, i] <- 0
      }
      max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
      if (cos_matrix[max][1] < similarity_limit) {
        fixed_df <- rbind(SBS1, reference[signatures, ])
        break
      }
    }
  }
  else {
    fixed_df <- SBS1
  }

  # denovo signatures ----------------------------------------------------------

  if (denovo_num > 1) {

    while (TRUE) {
      shuffled_denovo = denovo_catalogue[sample(1:nrow(denovo_catalogue)), ]
      signatures <- rownames(shuffled_denovo[1:denovo_num, ])
      cos_matrix <- denovo_cosine[signatures, signatures]
      for (i in 1:nrow(cos_matrix)) {
        cos_matrix[i, i] <- 0
      }
      max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
      if (cos_matrix[max][1] < similarity_limit) {
        denovo_df <- denovo_catalogue[signatures, ]
        rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
        break
      }
    }
  }
  else if (denovo_num==1) {
    shuffled_denovo = denovo_catalogue[sample(1:nrow(denovo_catalogue)), ]
    denovo_df <- shuffled_denovo[1, ]
    rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
  }
  else {
    denovo_df <- NULL
  }

  obj <- list(fixed = fixed_df, denovo = denovo_df)
  return(obj)
}

