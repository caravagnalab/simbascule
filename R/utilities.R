compute.mae = function(m , mr) {
  mae = sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}


compute.mse = function(m_inf, m_true, subset_cols=NULL, assigned_missing=NULL) {
                        # assigned=NULL, subset_cols=NULL) {
  if (!is.null(assigned_missing)) {
    m_true = m_true %>% dplyr::select(names(assigned_missing$assigned_tp))
    m_inf = m_inf %>% dplyr::select(assigned_missing$assigned_tp)
  }

  if (!is.null(subset_cols)) {
    m_true = m_true %>%
      dplyr::select(intersect(subset_cols, assigned_missing$assigned_tp))
    m_inf = m_inf %>%
      dplyr::select(intersect(subset_cols, assigned_missing$assigned_tp))
  }

  if (ncol(m_inf) != ncol(m_true)) {
    common = intersect(colnames(m_inf), colnames(m_true))
    m_inf = m_inf[, common]
    m_true = m_true[, common]
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

  m1 = as.data.frame(m1)
  m2 = as.data.frame(m2)

  if (what == "expos") {
    if (any(rownames(m1) != rownames(m2)))
      rownames(m2) = rownames(m1)

    m1 = as.data.frame(t(m1[, assigned_missing$assigned_tp]))
    m2 = as.data.frame(t(m2[, names(assigned_missing$assigned_tp)]))
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
}


