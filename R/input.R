#----------------------------------------------------------------------QC:PASSED

generate.input <- function(
    reference_catalogue,
    beta_fixed,
    complexity,
    seed=NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  
  k_fixed <- nrow(beta_fixed)
  
  if (is.null(complexity)) {
    return(NULL)
  } else if (is.numeric(complexity) & length(complexity)==2) {
    n_overlap <- complexity[1]
    n_extra <- complexity[2]
    if (n_overlap > k_fixed) {
      stop("overlap signatures are more than fixed signatures!")
    }
  } else if (is.character(complexity)) {
    if (complexity=='low') {
      n_overlap <- sample(1:k_fixed, 1)
      n_extra <- 0
    }
    else if (complexity=='medium') {
      n_overlap <- sample(1:k_fixed, 1)
      n_extra <- sample(1:k_fixed, 1)
    }
    else if (complexity=='high') {
      n_overlap <- 0
      n_extra <- sample(1:(k_fixed), 1)
    }
    else {
      stop("complexity argument should be selected from {'low', 'medium', 'high'}")
    }
  } else {
    stop("wrong complexity argument!")
  }
  
  extra_ref <- reference_catalogue[setdiff(rownames(reference_catalogue), rownames(beta_fixed)), ]
  
  if (n_overlap > 0) {
    overlap <- sample(rownames(beta_fixed))[1:n_overlap]
  } else {
    overlap <- NULL
  }
  
  if (n_extra > 0) {
    extra <- sample(rownames(extra_ref))[1:n_extra]
  } else {
    extra <- NULL
  }
  
  df <- reference_catalogue[c(overlap, extra), ]
  
  return(df)
}
