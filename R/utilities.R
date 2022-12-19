
#-------------------------------------------------------------------------------
cosine.vector <- function(a, b) {
  
  if (!identical(colnames(a), colnames(b))) {
    a = a[names(b)]
  }
  
  numerator <- sum(a * b)
  denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  return(numerator / denominator)
}

#-------------------------------------------------------------------------------
cosine.matrix <- function(a, b) {
  # a and b are data.frame
  
  df <- data.frame(matrix(0, nrow(a), nrow(b)))
  rownames(df) <- rownames(a)
  colnames(df) <- rownames(b)
  
  for (i in 1:nrow(a)) {
    denovo <- a[i, ]
    for (j in 1:nrow(b)) {
      ref <- b[j, ]
      
      score <- cosine.vector(denovo, ref)
      df[i,j] <- score
    }
  }
  
  return(df)
}

#----------------------------------------------------------------------QC:PASSED

fixed.accuracy <- function(reference, input, expected_fixed, inferred_fixed) {
  ref_list <- rownames(reference)
  if (is.null(expected_fixed)) {exp_list <- c()} else {exp_list <- rownames(expected_fixed)}
  if (is.null(inferred_fixed)) {inf_list <- c()} else {inf_list <- rownames(inferred_fixed)}
  if (is.null(input)) {input_list <- c()} else {input_list <- rownames(input)}
  
  TP <- length(intersect(inf_list, exp_list))
  FP <- length(setdiff(inf_list, exp_list))
  #TN <- length( setdiff( setdiff(ref_list, exp_list), inf_list) )
  TN <- length( setdiff(setdiff(input_list, exp_list), inf_list)  )
  FN <- length(setdiff(exp_list, inf_list))
  
  accuracy <- list(TP=TP, FP=FP, TN=TN, FN=FN)
  
  #accuracy <- (TP + TN) / (TP + TN + FP + FN)
  return(accuracy)
}

#----------------------------------------------------------------------QC:PASSED

reconstruct.count <- function(m, alpha, beta) {
  # all args are data.frame
  theta <- diag(rowSums(m))               # matrix
  alpha <- theta %*% as.matrix(alpha)     # matrix
  beta <- as.matrix(beta)                 # matrix
  
  mr_matrix <- alpha %*% beta
  mr <- round(as.data.frame(mr_matrix))
  rownames(mr) <- rownames(m)
  return(mr)
}

#----------------------------------------------------------------------QC:PASSED

compute.mae <- function(m , mr) {
  mae <- sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}

#----------------------------------------------------------------------QC:PASSED

compute.mse <- function(m , mr) {
  mse <- sum((m - mr)^2) / (dim(m)[1] * dim(m)[2])
  return(mse)
}

#----------------------------------------------------------------------QC:PASSED

denovo.similarity <- function(expected_denovo, inferred_denovo) {
  
  if (length(expected_denovo)==0 | length(inferred_denovo)==0) {
    return(NULL)
  } else {
    df <- data.frame(matrix(nrow = nrow(inferred_denovo), ncol = nrow(expected_denovo)))
    colnames(df) <- rownames(expected_denovo)
    rownames(df) <- rownames(inferred_denovo)
    
    for (i in 1:nrow(inferred_denovo)) {
      inferred <- inferred_denovo[i,]
      inferred_name <- rownames(inferred)
      for (j in 1:nrow(expected_denovo)) {
        target <- expected_denovo[j, ]
        target_name <- rownames(target)
        score <- cosine.vector(inferred, target)
        df[inferred_name, target_name] <- score
      }
    }
    
    #------------------------------
    #match_list <- list()
    match_df <- data.frame(matrix(nrow = nrow(inferred_denovo), ncol = 2))
    colnames(match_df) <- c("match", "similarity")
    rownames(match_df) <- rownames(inferred_denovo)
    
    similarity <- 0
    iter <- min(nrow(inferred_denovo), nrow(expected_denovo))
    for (i in 1:iter) {
      
      max = which(df == max(df), arr.ind = TRUE)
      similarity <- similarity + df[max]
      
      row <- row.names(df[max[,1],])
      column <- names(df[max[,2]])
      
      #match_list[row] <- column
      match_df[row, 'match'] <- column
      match_df[row, 'similarity'] <- df[max]
      
      df[row, column] <- 0
    }
    
    return( list( similarity_average=(similarity / iter), match_df=match_df ) )
  }
}

#----------------------------------------------------------------------QC:PASSED

denovo.ratio <- function(expected_denovo, inferred_denovo) {
  
  if (is.null(expected_denovo)) {n_exp <- 0} else {n_exp <- nrow(expected_denovo)}
  if (is.null(inferred_denovo)) {n_inf <- 0} else {n_inf <- nrow(inferred_denovo)}
  
  denovo_ratio <- (n_inf + 1) / (n_exp + 1)
  
  return(denovo_ratio)
}
