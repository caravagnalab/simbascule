generate_data_aux = function(N, G, catalogue_sbs,
                             alpha_range, alpha_sigma, seed,
                             n_muts_range = 500:5000,
                             pi_conc=1., frac_rare=1.,
                             shared_sbs=c("SBS1","SBS5"),
                             cohort="", out_path=NULL,
                             idd=NULL) {

  if (!is.null(idd))
    out_name = paste0("simul.", idd, ".", cohort, ".Rds") %>% stringr::str_replace_all("\\.\\.", ".")

  if (!is.null(out_path) && paste0(out_name) %in% list.files(paste0(out_path)))
    return(readRDS(paste0(out_path, out_name)))

  set.seed(seed)
  sbs = sample(setdiff(rownames(catalogue_sbs), shared_sbs), G + floor(G/2), replace=F)

  set.seed(seed)
  private_shared_sbs = sample(sbs, floor(G/2))
  private_sbs = setdiff(sbs, private_shared_sbs)

  dataa = generate_simulation_dataset(G=G, N=N, alpha_sigma=alpha_sigma,
                                      catalogue_sbs=catalogue_sbs, private_sbs=private_sbs,
                                      private_shared_sbs=private_shared_sbs,
                                      shared_sbs=shared_sbs, alpha_range=alpha_range,
                                      py=py, seed=seed)

  if (!is.null(out_path)) {
    if (!dir.exists(out_path))
      dir.create(out_path, recursive=T)

    if (cohort == "") saveRDS(dataa, paste0(out_path, "simul.", idd, ".Rds"))
    else saveRDS(dataa, paste0(out_path, "simul.", idd, ".", cohort, ".Rds"))
  }

  return(dataa)
}


generate_simulation_dataset = function(G, N,
                                       catalogue_sbs, private_sbs,
                                       private_shared_sbs, py,
                                       shared_sbs=c("SBS1","SBS5"),
                                       min_samples=2,
                                       n_muts_range=500:5000, frac_rare=1.,
                                       alpha_range=c(.15,.2), alpha_sigma=0.1,
                                       pi_conc=1, seed=1234, cohort="") {

  set.seed(seed)

  stopifnot(alpha_range[1] > alpha_sigma)

  n_rare = min(ceiling(frac_rare * N), N)
  min_n = min(min_samples, n_rare)
  repeat {
    pi = rdirichlet(1, alpha=rep(pi_conc, G))
    n_gs = sapply(1:G, function(i) round(N * pi[i]))
    if (all(n_gs >= min_n) && any(n_gs <= n_rare)) break
  }

  private_shared = data.frame()
  if (length(private_shared_sbs) > 0) {
    private_shared = lapply(1:floor(G/2), function(i)
      data.frame(group=sample(G, 2), idd=private_shared_sbs[i])
    ) %>% do.call(rbind, .)
  }

  shared = expand.grid(group=1:G, idd=shared_sbs)
  private = data.frame(group=1:G, idd=sample(private_sbs, G, replace=F))
  sbs_groups = rbind(private_shared, shared, private) %>%
    dplyr::mutate(idd=as.character(idd))

  counts = alpha = alpha_prior = tibble::tibble()
  beta = data.frame()
  for (gid in 1:G) {
    n_g = n_gs[gid]
    set.seed(seed)
    n_muts_g = reticulate::r_to_py(sample(n_muts_range, size=n_g))

    beta_g = catalogue_sbs[sbs_groups %>%
                             dplyr::filter(group==gid) %>%
                             dplyr::pull(idd), ]

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
                                             groupid=gid)) %>%
      replace(is.na(.), 0)
    alpha_prior = dplyr::bind_rows(alpha_prior, alpha_prior_g/sum(alpha_prior_g)) %>%
      replace(is.na(.), 0)
    beta = dplyr::bind_rows(beta, beta_g[setdiff(rownames(beta_g), rownames(beta)), ])
  }

  pl_alpha = alpha %>%
    reshape2::melt(id=c("sample","groupid")) %>%
    ggplot() + geom_bar(aes(x=sample, y=value, fill=variable), stat="identity") +
    facet_grid(~groupid, scales="free_x", space="free_x")

  return(tibble("counts"=list(counts %>%
                                dplyr::select(-groupid) %>%
                                tibble::column_to_rownames(var="sample")),
                "alpha"=list(alpha %>%
                               dplyr::select(-groupid) %>%
                               tibble::column_to_rownames(var="sample")),
                "groups"=list(counts$groupid),
                "alpha_prior"=list(alpha_prior), "beta"=list(beta),
                "sbs_groups"=list(tidyr::nest(sbs_groups, data=idd)),
                "shared"=list(shared_sbs),
                "private"=list(private_sbs),
                "private_shared"=list(private_shared_sbs),
                "alpha_plot"=list(pl_alpha))
         )
}
