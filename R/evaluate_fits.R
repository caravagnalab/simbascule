# empty_stats_df = function() {
#   num = as.numeric(NA)
#   chr = as.character(NA)
#   lgc = as.logical(NA)
#
#   stats = data.frame() %>%
#     dplyr::mutate(N=num, G=num,
#                   n_shared=num, n_private=num,
#                   mse_counts=num, mse_expos=num,
#                   cosine_sigs=num, cosine_expos=num,
#                   true_K=num, inf_K=num,
#                   idd=chr, is_hierarchical=lgc)
#   return(stats)
# }


compare_single_fit = function(fitname, fits_path, data_path, cutoff=0.8) {
  idd = stringr::str_replace_all(fitname, pattern="fit.|hier.|.Rds", replacement="")
  simulname = paste0("simul.", idd, ".Rds")

  x.simul = readRDS(paste0(data_path, simulname)) %>% create_basilica_obj_simul()
  x.fit = readRDS(paste0(fits_path, fitname))

  sigs.fit = get_signatures(x.fit)
  sigs.simul = get_signatures(x.simul)

  expos.fit = get_exposure(x.fit)
  expos.simul = get_exposure(x.simul)

  assigned = compare_sigs_inf_gt(sigs.fit, sigs.simul, cutoff=cutoff)
  unassigned = c(setdiff(rownames(sigs.fit), assigned),
                 setdiff(rownames(sigs.simul), names(assigned)))

  mse_counts = compute.mse(m_true=x.simul$input$counts, m_inf=get_data(x.fit, reconstructed=T))
  mse_expos = compute.mse(m_true=expos.simul, m_inf=expos.fit, assigned=assigned)
  cosine_sigs = compute.cosine(sigs.fit, sigs.simul, assigned, unassigned, what="sigs")
  cosine_expos = compute.cosine(expos.fit, expos.simul, assigned, unassigned, what="expos")

  true_n_sigs = nrow(sigs.simul)
  inf_n_sigs = nrow(sigs.fit)

  return(
    data.frame(
      "N"=(stringr::str_replace_all(idd, "N", "") %>% strsplit(split="[.]"))[[1]][1] %>%
        as.numeric(),
      "G"=(stringr::str_replace_all(idd, "G", "") %>% strsplit(split="[.]"))[[1]][2] %>%
        as.numeric(),
      "n_shared"=nrow(x.simul$input$reference_catalogue),
      "n_private"=nrow(x.simul$fit$denovo_signatures),
      "mse_counts"=mse_counts,
      "mse_expos"=mse_expos,
      "cosine_sigs"=cosine_sigs,
      "cosine_expos"=cosine_expos,
      "true_K"=true_n_sigs,
      "inf_K"=inf_n_sigs,
      "idd"=idd,
      "is_hierarchical"=grepl("hier", fitname)
    )
  )
}


compute.cosine = function(m1, m2, assigned, unassigned, what) {
  if (what == "expos") {
    m1 = as.data.frame(t(m1))
    m2 = as.data.frame(t(m2))
  }

  rownames(m1) = paste0("F_", rownames(m1))
  rownames(m2) = paste0("S_", rownames(m2))

  compare = rbind(m1[paste0("F_",assigned), ],
                  m2[paste0("S_",names(assigned)), ] )

  cosine_sim = lsa::cosine(t(compare))[paste0("F_",assigned),
                                       paste0("S_",names(assigned))]

  cosines_tmp = sapply(names(assigned), function(i)
    cosine_sim[paste0("F_",assigned[i]), paste0("S_",i)]
  )

  return(
    mean(c(cosines_tmp,
           rep(0, length.out=length(unassigned))))
  )
}


convert_sigs_names = function(x, assigned) {
  signames = rownames(get_signatures(x))
  signames_dn = intersect(rownames(x$fit$denovo_signatures), assigned)

  missing = c(setdiff(signames, assigned),
              setdiff(assigned, signames))
  new_names = c(names(assigned), missing)
  new_names_dn = c(names(assigned[assigned %in% signames_dn]), missing)

  x$fit$exposure = x$fit$exposure[, c(assigned, missing)]
  x$fit$denovo_signatures = x$fit$denovo_signatures[c(), ]
  colnames(x$fit$exposure) = new_names

}


compare_sigs_inf_gt = function(sigs.fit, sigs.simul, cutoff=0.8) {
  common = intersect(rownames(sigs.fit), rownames(sigs.simul))
  unique_inf = setdiff(rownames(sigs.fit), common)
  unique_gt = setdiff(rownames(sigs.simul), common)

  total_sigs = rbind(sigs.fit[!rownames(sigs.fit) %in% common,], sigs.simul)
  cosine_matr = lsa::cosine(t(total_sigs))[rownames(sigs.simul), rownames(sigs.fit)]

  assign_similar = cosine_matr %>% as.data.frame() %>%
    tibble::rownames_to_column(var="gt") %>%
    reshape2::melt(id="gt", variable.name="inf", value.name="cosine") %>%
    dplyr::filter(cosine >= cutoff) %>%
    dplyr::group_by(gt) %>%
    dplyr::mutate(inf=as.character(inf)) %>%
    dplyr::filter(cosine == max(cosine)) %>% dplyr::arrange(gt)

  # if (nrow(sigs.simul) > nrow(sigs.fit))
  if (any(duplicated(assign_similar$inf)))
    assign_similar = assign_similar %>% dplyr::group_by(inf) %>%
      dplyr::filter(cosine == max(cosine)) %>% ungroup()

  similar = assign_similar$inf %>% setNames(assign_similar$gt)

  return(similar)
}


