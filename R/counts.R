#----------------------------------------------------------------------QC:PASSED
# may use theta*alpha*beta or generate.counts (should be discussed to )
generate.counts <- function(alpha, beta, theta, seed=NULL) {
  
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  
  if (!is.null(alpha$group)) {
    alpha <- subset(alpha, select = -c(group))
  }
  
  alpha <- alpha[, order(colnames(alpha))]
  beta <- beta[order(rownames(beta)), ]
  if (!identical(colnames(alpha), rownames(beta))) {
    print(colnames(alpha))
    print(rownames(beta))
    stop("alpha and beta are NOT valid!")
  }
  
  num_samples <- nrow(alpha)
  
  M <- matrix(rep(0, num_samples*96) , ncol = 96)
  
  # iterate over samples
  for (sample in 1:num_samples) {
    
    p <- as.numeric(alpha[sample, ]) # select sample i
    
    # iterate over number of mutations in sample i
    for (j in 1:theta[sample]) {
      
      # sample signature profile index from categorical data
      signature_idx <- extraDistr::rcat(1, p)
      signature <- beta[signature_idx, ]
      
      # sample mutation feature index for corresponding signature from categorical data
      mutation_idx <- extraDistr::rcat(1, as.numeric(signature))
      
      # add +1 to the mutation feature in position j in branch i
      M[sample, mutation_idx] <- M[sample, mutation_idx] + 1
      
    }
  }
  M <- as.data.frame(M)
  colnames(M) <- colnames(beta)
  rownames(M) <- rownames(alpha)
  return(M)
}
