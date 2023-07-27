cosine.vector <- function(a, b) {

  if (!identical(colnames(a), colnames(b))) {
    a = a[names(b)]
  }

  numerator <- sum(a * b)
  denominator <- sqrt(sum(a^2)) * sqrt(sum(b^2))
  return(numerator / denominator)
}


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


edit.exposure <- function(alpha) {

  n <- nrow(alpha)
  val <- 0.05

  while ( TRUE ) {

    exp <- alpha

    # fixed
    newFixed <- runif(n, 0, val)
    fixed_diff <- exp[, 4] - newFixed
    exp[, 4] <- newFixed
    exp[, 3] <- exp[, 3] + fixed_diff

    # de-novo
    newDenovo <- runif(n, 0, val)
    denovo_diff <- exp[, 7] - newDenovo
    exp[, 7] <- newDenovo
    exp[, 6] <- exp[, 6] + denovo_diff

    if ( sum(rowSums(exp))==n & all(apply(exp, 1, function(x) all(x > 0))) ) {
      return(exp)
    }
  }
}


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


compute.mae <- function(m , mr) {
  mae <- sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}


compute.mse <- function(m_inf, m_true, subset_cols=NULL, assigned_missing=NULL) {
                        # assigned=NULL, subset_cols=NULL) {
  if (!is.null(assigned_missing)) {
    # m_true[, assigned_missing$added_fp] = 0
    # m_inf[, assigned_missing$missing_fn] = 0
    m_true = m_true %>% dplyr::select(assigned_missing$assigned_tp)
    m_inf = m_inf %>% dplyr::select(assigned_missing$assigned_tp)
  }

  if (!is.null(subset_cols)) {
    m_true = m_true %>%
      dplyr::select(intersect(subset_cols, assigned_missing$assigned_tp))
    m_inf = m_inf %>%
      dplyr::select(intersect(subset_cols, assigned_missing$assigned_tp))
  }

  mse = sum((m_inf - m_true)^2) / (dim(m_inf)[1] * dim(m_inf)[2])
  mse = sqrt(mse) / ( max(m_true) - min(m_true) )
  return(mse)
}



compute.cosine = function(m1, m2, assigned_missing, what, subset_cols=NULL) {
  # m1 = fit
  # m2 = simul

  unassigned = unique(c(assigned_missing$missing_fn,
                        assigned_missing$added_fp))

  if (what == "expos") {
    # m1[, assigned_missing$missing_fn] = 1e-10
    # m2[, assigned_missing$added_fp] = 1e-10
    m1 = as.data.frame(t(m1[, assigned_missing$assigned_tp]))
    m2 = as.data.frame(t(m2[, assigned_missing$assigned_tp]))
  }

  if (!is.null(subset_cols)) {
    m1 = m1[intersect(rownames(m1), assigned_missing$assigned_tp[subset_cols]), ]
    m2 = m2[intersect(rownames(m2), subset_cols), ]

    if (nrow(m1) == 0) {
      return(0)
    }
    tmp = intersect(subset_cols, names(assigned_missing$assigned_tp))
    consider = assigned_missing$assigned_tp[tmp] %>%
      setNames(tmp)
    unassigned = intersect(unassigned, subset_cols)
  } else if (what == "expos") {
    # consider = c(assigned_missing$assigned_tp, unassigned) %>%
    #   setNames(c(names(assigned_missing$assigned_tp), unassigned))
    consider = c(assigned_missing$assigned_tp) %>%
      setNames(c(names(assigned_missing$assigned_tp)))
  } else {
    consider = assigned_missing$assigned_tp
  }

  rownames(m1) = paste0("F_", rownames(m1))
  rownames(m2) = paste0("S_", rownames(m2))

  compare = rbind(m1[paste0("F_",consider), ],
                  m2[paste0("S_",names(consider)), ] )

  cosine_sim = lsa::cosine(t(compare))[paste0("F_",consider),
                                       paste0("S_",names(consider))]

  if (is.null(dim(cosine_sim)))
    cosines_tmp = cosine_sim
  else
    cosines_tmp = sapply(names(consider), function(i) {
      cosine_sim[paste0("F_",consider[i]), paste0("S_",i)]
    })

  return(mean(cosines_tmp))

  # return(
  #   mean(c(cosines_tmp,
  #          rep(0, length.out=length(unassigned))))
  # )
}



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

      match_df[row, 'match'] <- column
      match_df[row, 'similarity'] <- df[max]

      df[row, column] <- 0
    }

    return( list( similarity_average=(similarity / iter), match_df=match_df ) )
  }
}


denovo.ratio <- function(expected_denovo, inferred_denovo) {

  if (is.null(expected_denovo)) {n_exp <- 0} else {n_exp <- nrow(expected_denovo)}
  if (is.null(inferred_denovo)) {n_inf <- 0} else {n_inf <- nrow(inferred_denovo)}

  denovo_ratio <- (n_inf + 1) / (n_exp + 1)

  return(denovo_ratio)
}
