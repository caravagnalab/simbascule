generate_simulation_dataset = function(G, N,
                                       sbs_catalogue, private_sbs,
                                       private_shared_sbs, py,
                                       shared_sbs=c("SBS1","SBS5"),
                                       n_muts_range=500:5000, frac_rare=1.,
                                       alpha_range=c(.15,0.2), alpha_sigma=0.1,
                                       pi_conc=1, seed=1234) {

  set.seed(seed)

  stopifnot(alpha_range[1] > alpha_sigma)

  n_rare = min(ceiling(frac_rare * N), N)
  min_n = min(10, n_rare)
  repeat {
    pi = rdirichlet(1, alpha=rep(pi_conc, G))
    n_gs = sapply(1:G, function(i) round(N * pi[i]))
    if (all(n_gs >= min_n) && any(n_gs <= n_rare)) break
  }

  private_shared = lapply(1:floor(G/2), function(i) data.frame(group=sample(G, 2), idd=private_shared_sbs[i])) %>%
    do.call(rbind, .)
  shared = expand.grid(group=1:G, idd=shared_sbs)
  private = data.frame(group=1:G, idd=sample(private_sbs, G, replace=F))
  sbs_groups = rbind(private_shared, shared, private)


  counts = alpha = alpha_prior = tibble::tibble()
  beta = data.frame()
  for (gid in 1:G) {
    n_g = n_gs[gid]
    set.seed(seed)
    n_muts_g = reticulate::r_to_py(sample(n_muts_range, size=n_g))

    beta_g = sbs_catalogue[c(shared_sbs,
                             private %>% dplyr::filter(group==gid) %>% dplyr::pull(idd),
                             private_shared %>% dplyr::filter(group==gid) %>% dplyr::pull(idd)), ]
    # private_sbs = setdiff(private_sbs, rownames(beta_g))
    # private_shared_sbs = setdiff(private_shared_sbs, rownames(beta_g))

    repeat {
      alpha_prior_g = sample(1:100, size=nrow(beta_g)) %>% setNames(rownames(beta_g))
      alpha_prior_g = alpha_prior_g / sum(alpha_prior_g)
      w_min = min(alpha_prior_g)

      if (w_min >= alpha_range[1] && w_min <= alpha_range[2]) break
    }

    alpha_sigma_g = min(min(alpha_prior_g) - 0.02, alpha_sigma)

    data_g = py$generate_model(reticulate::r_to_py(alpha_prior_g),
                               beta_g, n_muts_g, N=n_g,
                               seed=as.integer(seed), use_normal=T,
                               alpha_sigma=alpha_sigma_g)

    # data_g$data$group = gid
    counts = dplyr::bind_rows(counts, data_g$data %>%
                                tibble::rownames_to_column(var="sample") %>%
                                dplyr::mutate(sample=paste0("G",gid,"_",sample),
                                              groupid=gid))
    alpha = dplyr::bind_rows(alpha, data_g$alpha %>%
                               tibble::rownames_to_column(var="sample") %>%
                               dplyr::mutate(sample=paste0("G",gid,"_",sample),
                                             groupid=as.character(gid))) %>%
      replace(is.na(.), 0)
    alpha_prior = dplyr::bind_rows(alpha_prior, alpha_prior_g/sum(alpha_prior_g)) %>%
      replace(is.na(.), 0)
    beta = dplyr::bind_rows(beta, beta_g[setdiff(rownames(beta_g), rownames(beta)), ])
  }

  pl_alpha = alpha %>%
    # tibble::rownames_to_column(var="sample") %>%
    reshape2::melt() %>%
    ggplot() + geom_bar(aes(x=sample, y=value, fill=variable), stat="identity") +
    facet_grid(~groupid, scales="free_x", space="free_x")

  return(tibble("counts"=list(counts), "alpha"=list(alpha),
                "alpha_prior"=list(alpha_prior), "beta"=list(beta),
                "sbs_groups"=list(tidyr::nest(sbs_groups, data=idd)),
                "alpha_plot"=list(pl_alpha))
         )
}
