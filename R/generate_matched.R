generate_simulation_dataset_matched = function(N, G, private, py,
                                               shared=list("SBS"=c("SBS1","SBS5"),
                                                           "DBS"=c("DBS2","DBS4")),
                                               reference=list("SBS"=COSMIC_filt,
                                                              "DBS"=COSMIC_dbs),
                                               min_samples=2,
                                               n_muts_range=500:5000,
                                               frac_rare=1.,
                                               alpha_range=c(.15,.2),
                                               alpha_sigma=0.1,
                                               pi_conc=1,
                                               seed=list("SBS"=56,"DBS"=67)) {

  types = names(shared)

  repeat{
    pi = gtools::rdirichlet(1, alpha=rep(1,G)) %>% as.numeric()
    n_g = round(pi*N)
    if (all(n_g > 5)) break
  }

  private_df = select_private_signatures(G=G, types=types, private=private)

  used_sigs = c()
  counts_groups = lapply(1:G, function(gid) {
    lapply(types, function(tid) {
      keep_sigs = setdiff(rownames(reference[[tid]]), used_sigs)
      private_t = private_df %>% dplyr::filter(group==gid, type==tid) %>%
        dplyr::pull(idd)

      obj_tid = generate_simulation_dataset(N=n_g[gid], G=1,
                                            private_sbs=private_t,
                                            private_shared_sbs=c(),
                                            shared_sbs=shared[[tid]],
                                            use_all_sigs=TRUE,
                                            catalogue_sbs=reference[[tid]],
                                            min_samples=min_samples,
                                            n_muts_range=n_muts_range,
                                            frac_rare=frac_rare,
                                            alpha_range=alpha_range,
                                            alpha_sigma=alpha_sigma,
                                            pi_conc=pi_conc,
                                            seed=seed[[tid]]*gid, py=py)
      used_sigs <<- c(used_sigs, obj_tid$beta[[1]] %>% rownames())
      return(obj_tid)
    }) %>% setNames(types)
  })

  counts = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      counts_g = counts_groups[[gid]][[tid]]$counts[[1]]
      rownames(counts_g) = stringr::str_replace_all(rownames(counts_g), "G1", paste0("G",gid))
      return(counts_g)
    }) %>% do.call(rbind, .)
  }) %>% setNames(types)


  betas = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      betas_g = counts_groups[[gid]][[tid]]$beta[[1]] # %>% tibble::rownames_to_column("signame")
    }) %>% do.call(rbind, .) %>% unique()
  }) %>% setNames(types)


  exposures = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      alphas_g = counts_groups[[gid]][[tid]]$alpha[[1]]
      rownames(alphas_g) = stringr::str_replace_all(rownames(alphas_g), "G1", paste0("G",gid))
      return(alphas_g)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(dplyr::across(dplyr::everything(),
                                                              function(i) replace(i, is.na(i), 0)))
  }) %>% setNames(types)

  centroids = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      centroids_g = counts_groups[[gid]][[tid]]$alpha_prior[[1]] %>% data.frame()
      rownames(centroids_g) = paste0("G",gid)
      return(centroids_g)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(dplyr::across(dplyr::everything(),
                                                              function(i) replace(i, is.na(i), 0)))
  }) %>% setNames(types)

  return(tibble::tibble("counts"=counts,
                        "betas"=betas,
                        "exposures"=exposures,
                        "centroids"=centroids,
                        "private_sigs"=list(private_df),
                        "types"=types))
}


select_private_signatures = function(G, types, private) {
  n_privates = lapply(types, function(tid) {
    n_private_tid = max(floor(G/2),1)
    repeat{
      if (length(private[[tid]]) >= n_private_tid*G) break
      n_private_tid = n_private_tid - 1
    }
    return(max(n_private_tid, 1))
  }) %>% setNames(types)

  n_grps_with_private = lapply(types, function(tid) {
    if (n_privates[[tid]] == 1) return(min(G, length(private[[tid]])))
    return(G)
  }) %>% setNames(types)

  private_df = lapply(types, function(tid) {
    private_tid = private[[tid]]
    if (length(private_tid) == 0) return(data.frame())
    grps = sample(1:G, size=n_grps_with_private[[tid]])
    repeat {
      private_tid_df = lapply(grps, function(gid)
        data.frame(group=gid,
                   idd=sample(x=private_tid, size=n_privates[[tid]], replace=F),
                   type=tid)
      ) %>% do.call(rbind, .)
      if (check_private_overlaps(private_tid_df, n_privates=n_privates[[tid]], type=tid)) break
    }
    return(private_tid_df)
  }) %>% do.call(rbind, .)

  return(private_df)
}


check_private_overlaps = function(private_df, n_privates, type) {
  G_list = private_df$group %>% unique() %>% sort()
  for (gid1 in G_list) {
    for (gid2 in G_list) {
      if (gid1 == gid2) next
      prv_g1 = private_df %>% dplyr::filter(group==gid1, type==type) %>% dplyr::pull(idd)
      prv_g2 = private_df %>% dplyr::filter(group==gid2, type==type) %>% dplyr::pull(idd)
      if (n_privates==1 && length(intersect(prv_g1, prv_g2)) == 1) return(FALSE)
      if (length(intersect(prv_g1, prv_g2)) > 1) return(FALSE)
    }
  }
  return(TRUE)
}


create_bascule_obj_simul = function(simul_df) {
  types = simul_df$types
  obj_simul = list(); class(obj_simul) = "bascule_obj"

  obj_simul$prv_sigs = simul_df$private_sigs

  obj_simul$input = lapply(types, function(tid) {
    list("counts"=simul_df$counts[[tid]] %>% wide_to_long(what="counts"))
  }) %>% setNames(types)

  obj_simul$nmf = lapply(types, function(tid) {
    list("exposure"=simul_df$exposures[[tid]] %>% wide_to_long(what="exposures"),
         "beta_denovo"=simul_df$betas[[tid]] %>% wide_to_long(what="beta"))
  }) %>% setNames(types)

  groups = rownames(simul_df$counts[[1]]) %>% strsplit("_") %>%
    sapply(function(i) return(i[[1]][1]))

  obj_simul$clustering = list(
    "clusters"=tibble::tibble("samples"=obj_simul$input[[1]]$counts %>%
                                dplyr::select(samples) %>% unique() %>%
                                dplyr::pull(samples),
                              "clusters"=groups)
  )

  return(obj_simul)
}
