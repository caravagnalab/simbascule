devtools::load_all("~/GitHub/simbasilica/")
load_deps()


G = 2; N = 150
repeat{
  pi = gtools::rdirichlet(1, alpha=rep(1,G)) %>% as.numeric()
  n_g = round(pi*N)
  if (all(n_g > 5)) break
}

private_sbs=c("SBS17b"); private_shared_sbs=c("SBS4"); shared_sbs=c("SBS1","SBS5")
private_dbs=c("DBS7"); private_shared_dbs=c("DBS11","DBS5"); shared_dbs=c("DBS3")

used_sbs = used_dbs = c()
counts_groups = lapply(1:G, function(gid) {
  sbs_keep = setdiff(rownames(COSMIC_filt), used_sbs)
  dbs_keep = setdiff(rownames(COSMIC_dbs), used_dbs)

  print(sbs_keep)
  print(dbs_keep)

  sbs = generate_simulation_dataset(N=n_g[gid], G=1,
                                    private_sbs=private_sbs[private_sbs %in% sbs_keep],
                                    private_shared_sbs=private_shared_sbs[private_shared_sbs %in% sbs_keep],
                                    shared_sbs=shared_sbs,
                                    catalogue_sbs=COSMIC_filt,
                                    seed=54637, py=py)

  dbs = generate_simulation_dataset(N=n_g[gid], G=1,
                                    private_sbs=private_dbs[private_dbs %in% dbs_keep],
                                    private_shared_sbs=private_shared_dbs[private_shared_dbs %in% dbs_keep],
                                    shared_sbs=shared_dbs,
                                    catalogue_sbs=COSMIC_dbs,
                                    seed=4739, py=py)

  used_sbs <<- c(used_sbs, sbs$beta[[1]] %>% rownames())
  used_dbs <<- c(used_dbs, dbs$beta[[1]] %>% rownames())

  return(tibble::tibble("SBS"=list(sbs), "DBS"=list(dbs)))
})


types = c("SBS","DBS")
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

# plot_exposures_aux(exposures$SBS %>% wide_to_long(what="exposures") %>% dplyr::mutate(type="SBS"))
# plot_exposures_aux(exposures$DBS %>% wide_to_long(what="exposures") %>% dplyr::mutate(type="DBS"))


test = fit(counts=counts, k_list=4, cluster=6, n_steps=3,
           reference_cat=list("SBS"=COSMIC_filt[c("SBS1","SBS5"),], "DBS"=COSMIC_dbs["DBS4",]),
           keep_sigs=c("SBS1","SBS5","DBS4"),
           hyperparameters=list("scale_factor_centroid"=5000,
                                "scale_factor_alpha"=5000, "tau"=0),
           seed_list=c(10,33,4), filter_dn=TRUE, store_fits=TRUE)

get_alternative_run(x, params=list("seed"=10))


x = readRDS("~/GitHub/simbasilica/nobuild/analysis_multiple_signals/ex.Rds")

# x = fit(counts=counts, k_list=3:5, cluster=6, n_steps=3000,
#         reference_cat=list("SBS"=COSMIC_filt[c("SBS1","SBS5"),], "DBS"=COSMIC_dbs["DBS4",]),
#         keep_sigs=c("SBS1","SBS5","DBS4"),
#         hyperparameters=list("scale_factor_centroid"=5000,
#                              "scale_factor_alpha"=5000, "tau"=0),
#         seed_list=c(10,33,4), filter_dn=TRUE, store_fits=TRUE)
# saveRDS(x, "~/GitHub/simbasilica/nobuild/analysis_multiple_signals/ex.Rds")




test$clustering$pyro$alternatives$runs_seed$`seed:10`

test$nmf$SBS$pyro$alternatives$all_fits














