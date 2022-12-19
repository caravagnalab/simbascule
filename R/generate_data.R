#----------------------------------------------------------------------QC:PASSED
generate.data <- function(
    reference_catalogue,
    denovo_catalogue,
    reference_cosine,
    denovo_cosine,
    targetX,
    inputX,
    similarity_limit,
    groups,
    mut_range,
    seed=NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  
  # SIGNATURES --------------------------------
  signatures <- generate.signatures(
    reference_catalogue = reference_catalogue,
    denovo_catalogue  = denovo_catalogue,
    reference_cosine = reference_cosine,
    denovo_cosine = denovo_cosine,
    complexity = targetX,
    similarity_limit = similarity_limit,
    seed = seed
  )
  beta <- rbind(signatures$fixed, signatures$denovo)
  
  # INPUT ------------------------------------
  input <- generate.input(
    reference_catalogue = reference_catalogue,
    beta_fixed = signatures$fixed,
    complexity = inputX,
    seed = seed
  )
  
  # EXPOSURE ---------------------------------
  alpha <- generate.exposure(beta=beta, groups=groups, seed=seed) # include group column
  
  # THETA ------------------------------------
  num_samples <- length(groups)
  theta <- generate.theta(mut_range=mut_range, num_samples=num_samples, seed=seed)
  
  # COUNT MATRIX -----------------------------
  #m <- generate.counts(alpha=alpha, beta=beta, theta=theta, seed=seed)
  
  # removing group column
  if (!is.null(alpha$group)) {
    alpha <- subset(alpha, select = -c(group))
  }
  M <- as.data.frame(round(as.matrix(alpha*theta) %*% as.matrix(beta), digits = 0))
  rownames(M) <- rownames(alpha)
  colnames(M) <- colnames(beta)
  
  # MODIFY COMPLEXITY VALUES -----------------
  if (is.numeric(targetX)) {
    targetX <- paste(targetX, collapse = "|")
  }
  if (is.numeric(inputX)) {
    inputX <- paste(inputX, collapse = "|")
  }
  
  # CREATE TIBBLE ----------------------------
  obj <- tibble::tibble(
    x = list(M),
    input_cat = list(input),
    ref_cat = list(reference_catalogue),
    exp_exposure = list(alpha),
    exp_fixed = list(signatures$fixed),
    exp_denovo = list(signatures$denovo),
    targetX = targetX,
    inputX = inputX,
  )
  return(obj)
}
