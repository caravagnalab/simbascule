compute.mae = function(m , mr) {
  mae = sum(abs(m - mr)) / (dim(m)[1] * dim(m)[2])
  return(mae)
}


compute.mse = function(m_inf, m_true, subset_cols=NULL, assigned_missing=NULL,
                       keep_missing=FALSE) {
                        # assigned=NULL, subset_cols=NULL) {
  if (!is.null(assigned_missing) &&
      !keep_missing &&
      length(assigned_missing$assigned_tp) == 0) return(1)
  if (!is.null(assigned_missing)) {
    if (keep_missing) {
      m_true[, assigned_missing$added_fp] = 0
      m_inf[, assigned_missing$missing_fn] = 0
      m_true = m_true %>% dplyr::select(names(assigned_missing$assigned_tp),
                                        assigned_missing$added_fp,
                                        assigned_missing$missing_fn)
      m_inf = m_inf %>% dplyr::select(assigned_missing$assigned_tp %>% setNames(NULL),
                                      assigned_missing$added_fp,
                                      assigned_missing$missing_fn)
    } else {
      m_true = dplyr::select(m_true, names(assigned_missing$assigned_tp))
      m_inf = dplyr::select(m_inf, assigned_missing$assigned_tp)
    }
  } else if (ncol(m_inf) != ncol(m_true)) {
    common = intersect(colnames(m_inf), colnames(m_true))
    m_inf = m_inf[, common]
    m_true = m_true[, common]
  }

  if (!is.null(subset_cols)) {
    m_true = dplyr::select(m_true, intersect(subset_cols, assigned_missing$assigned_tp))
    m_inf = dplyr::select(m_inf, intersect(subset_cols, assigned_missing$assigned_tp))
  }

  mse = sum((m_inf - m_true)^2) / (dim(m_inf)[1] * dim(m_inf)[2])
  mse = sqrt(mse) / ( max(m_true) - min(m_true) )
  return(mse)
}



compute.cosine = function(m_true, m_inf, assigned_missing, what, subset_cols=NULL,
                          keep_missing=FALSE) {
  # m1 = fit
  # m2 = simul

  if (!keep_missing && (length(assigned_missing$assigned_tp) == 0)) return(0)

  unassigned = unique(c(assigned_missing$missing_fn, assigned_missing$added_fp))

  m_true = as.data.frame(m_true); m_inf = as.data.frame(m_inf)

  if (what == "expos") {
    if (any(rownames(m_inf) != rownames(m_true))) rownames(m_true) = rownames(m_inf)

    if (keep_missing) {
      m_true[, assigned_missing$added_fp] = 1e-10
      m_inf[, assigned_missing$missing_fn] = 1e-10
      m_true = m_true %>% dplyr::select(names(assigned_missing$assigned_tp),
                                        assigned_missing$added_fp,
                                        assigned_missing$missing_fn)
      m_inf = m_inf %>% dplyr::select(assigned_missing$assigned_tp %>% setNames(NULL),
                                      assigned_missing$added_fp,
                                      assigned_missing$missing_fn)
      m_true = as.data.frame(t(m_true)); m_inf = as.data.frame(t(m_inf))
    } else {
      m_inf = m_inf %>%
        dplyr::select(dplyr::matches(paste0("^",assigned_missing$assigned_tp,"$"))) %>%
        t() %>% as.data.frame()
      m_true = m_true %>%
        dplyr::select(dplyr::matches(paste0("^",names(assigned_missing$assigned_tp),"$"))) %>%
        t() %>% as.data.frame()

    }
  }

  if (!is.null(subset_cols)) {
    m_inf = m_inf[intersect(rownames(m_inf), assigned_missing$assigned_tp[subset_cols]), ]
    m_true = m_true[intersect(rownames(m_true), subset_cols), ]

    if (nrow(m_inf) == 0) return(0)

    tmp = intersect(subset_cols, names(assigned_missing$assigned_tp))
    consider = assigned_missing$assigned_tp[tmp] %>% setNames(tmp)
    unassigned = intersect(unassigned, subset_cols)
  } else if (what == "expos") {
    if (keep_missing) {
      consider = c(assigned_missing$assigned_tp,assigned_missing$added_fp,assigned_missing$missing_fn) %>%
        setNames(c(names(assigned_missing$assigned_tp),assigned_missing$added_fp,assigned_missing$missing_fn))
    } else {
      consider = c(assigned_missing$assigned_tp) %>% setNames(c(names(assigned_missing$assigned_tp)))
      }
  } else {
     consider = assigned_missing$assigned_tp
  }

  rownames(m_inf) = paste0("F_", rownames(m_inf)); rownames(m_true) = paste0("S_", rownames(m_true))

  compare = rbind(m_inf[paste0("F_",consider), ], m_true[paste0("S_",names(consider)), ] )

  cosine_sim = lsa::cosine(t(compare))[paste0("F_",consider), paste0("S_",names(consider))]
  cosine_sim[is.na(cosine_sim)] = 0

  if (is.null(dim(cosine_sim))) cosines_tmp = cosine_sim
  else cosines_tmp = sapply(names(consider), function(i) {
      cosine_sim[paste0("F_",consider[i]), paste0("S_",i)]
    })

  return(mean(cosines_tmp))
}


