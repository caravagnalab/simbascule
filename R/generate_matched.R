generate_simulation_dataset_matched = function(N, G, private, py,
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
      set.seed(seed[[tid]])
      private_t = private[[tid]][private[[tid]] %in% keep_sigs] %>%
        sample(size=1, replace=F)

      obj_tid = generate_simulation_dataset(N=n_g[gid], G=1,
                                            private_sbs=private_t,
                                            private_shared_sbs=c(),
                                            shared_sbs=shared[[tid]],
                                            catalogue_sbs=reference[[tid]],
                                            seed=seed[[tid]], py=py)
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
                        "types"=types))
}



create_basilica_obj_simul = function(simul_df) {
  types = simul_df$types
  obj_simul = list()

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

  # lapply(types, function(tid) {
  #   centr = centroids[[tid]] %>% tibble::rownames_to_column(var="clusters") %>%
  #     reshape2::melt(variable.name="sigs")
  # })

  return(obj_simul)
}
