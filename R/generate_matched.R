generate_simulation_dataset_matched = function(N, G, private, private_shared, py,
                                               shared=list("SBS"=c("SBS1","SBS5"),
                                                           "DBS"=c("DBS2","DBS4")),
                                               reference=list("SBS"=COSMIC_filt,
                                                              "DBS"=COSMIC_dbs),
                                               seed=list("SBS"=56,"DBS"=67)) {

  types = names(shared)

  repeat{
    pi = gtools::rdirichlet(1, alpha=rep(1,G)) %>% as.numeric()
    n_g = round(pi*N)
    if (all(n_g > 5)) break
  }

  used_sigs = c()
  counts_groups = lapply(1:G, function(gid) {
    lapply(types, function(tid) {
      keep_sigs = setdiff(rownames(reference[[tid]]), used_sigs)
      private_t = private[[tid]][private[[tid]] %in% keep_sigs]
      private_shared_t = private_shared[[tid]][private_shared[[tid]] %in% keep_sigs]

      obj_tid = generate_simulation_dataset(N=n_g[gid], G=1,
                                            private_sbs=private_t,
                                            private_shared_sbs=private_shared_t,
                                            shared_sbs=shared[[tid]],
                                            catalogue_sbs=reference[[tid]],
                                            seed=seed[[tid]], py=py)
      used_sigs <<- c(used_sigs, obj_tid$beta[[1]] %>% rownames())
    }) %>% setNames(types)

    return(tibble::tibble("SBS"=list(sbs), "DBS"=list(dbs)))
  })

  counts = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      counts_g = counts_groups[[gid]][[tid]][[1]]$counts[[1]]
      rownames(counts_g) = stringr::str_replace_all(rownames(counts_g), "G1", paste0("G",gid))
      return(counts_g)
    }) %>% do.call(rbind, .)
  }) %>% setNames(types)


  betas = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      betas_g = counts_groups[[gid]][[tid]][[1]]$beta[[1]] %>% tibble::rownames_to_column("signame")
      return(betas_g)
    }) %>% do.call(rbind, .) %>% unique() %>% tibble::column_to_rownames(var="signame")
  }) %>% setNames(types)


  exposures = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      alphas_g = counts_groups[[gid]][[tid]][[1]]$alpha[[1]]
      rownames(alphas_g) = stringr::str_replace_all(rownames(alphas_g), "G1", paste0("G",gid))
      return(alphas_g)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(dplyr::across(dplyr::everything(),
                                                              function(i) replace(i, is.na(i), 0)))
  }) %>% setNames(types)


  centroids = lapply(types, function(tid) {
    lapply(1:G, function(gid) {
      centroids_g = counts_groups[[gid]][[tid]][[1]]$alpha_prior[[1]] %>% data.frame()
      rownames(centroids_g) = paste0("G",gid)
      return(centroids_g)
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(dplyr::across(dplyr::everything(),
                                                              function(i) replace(i, is.na(i), 0)))
  }) %>% setNames(types)

  return(tibble::tibble("counts"=counts,
                        "betas"=betas,
                        "exposures"=exposures,
                        "centroids"=centroids,
                        "types"=types))
}



create_basilica_obj_simul = function(simul_df) {
  types = simul_df$types[[1]]
  obj_simul = list()

  obj_simul$input = lapply(types, function(tid) {
    list("counts"=simul_df$counts[[tid]] %>% wide_to_long(what="counts"))
  }) %>% setNames(types)

  obj_simul$nmf = lapply(types, function(tid) {
    list("exposure"=simul_df$exposures[[tid]] %>% wide_to_long(what="exposures"),
         # "beta_fixed"=NULL,
         "beta_denovo"=simul_df$betas[[tid]] %>% wide_to_long(what="beta"))
  }) %>% setNames(types)

  obj_simul$clustering = list(
    "clusters"=tibble::tibble(samples=simul_df$counts[[1]][[1]]$samples,
                              clusters=paste0("G",init_params$init_clusters)),
    "centroids"
  )

  lapply(types, function(tid) {
    centr = centroids[[tid]] %>% tibble::rownames_to_column(var="clusters") %>%
      reshape2::melt(variable.name="sigs")
  })
}
