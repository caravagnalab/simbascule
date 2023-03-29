# library(lsa)

## Function to select "n_fixed" signatures from a reference catalogue (i.e. COSMIC)
## with a cosine similarity lower than "cosine_limit"
## It returns the subset reference
select_fixed_sbs = function(referece_cat, n_fixed=4, cosine_limit=.5) {
  ref_cosine.all = lsa::cosine(reference_cat %>% t())

  while (TRUE) {
    idxs = sample(1:nrow(ref_cosine.all), n_fixed)
    sbss = rownames(ref_cosine.all)[idxs]
    sbs_cosine = ref_cosine.all[sbss, sbss]
    if (any(sbs_cosine != 1. & sbs_cosine > cosine_limit)) next

    sbs_fixed = sbss
    sbs_ref = reference_cat[sbs_fixed,]
    break
  }
  return(sbs_ref)
}


single_dataset = function(N, n_groups, samples_per_group,
                          reference_cat, denovo_cat,
                          ref_cosine, denovo_cosine,
                          private_sigs, private_fracs,
                          cosine_limit, seed,
                          out_path=NULL) {

  groups = sample(1:n_groups, N, replace=T)
  while (!all(lapply(1:n_groups, function(n) length(groups[groups==n]) %in% samples_per_group) %>% unlist()))
    groups = sample(1:n_groups, N, replace=T)

  idd = paste0("N", N, ".G", n_groups, ".s", seed)

  x = generate.data(
    reference_catalogue=reference_cat,
    denovo_catalogue=denovo_cat,
    reference_cosine=ref_cosine,
    denovo_cosine=denovo_cosine,
    targetX=-1,
    inputX=NULL,
    similarity_limit=cosine_limit,
    groups=groups,
    private_sigs=private_sigs,
    private_fracs=private_fracs,
    mut_range=1:5000,
    seed=seed)

  if (is.null(out_path)) return(x)

  if (!dir.exists(out_path))
    dir.create(out_path, recursive=T)

  saveRDS(x, paste0(out_path, "simul.", idd, ".Rds"))
  return(x)
}




plot_alpha = function(alpha)
  return(
    alpha %>%
      reshape2::melt(id="group", variable.name="sbs", value.name="alpha") %>%
      # dplyr::mutate(type=dplyr::case_when(
      #   sbs %in% private_rare ~ "private_rare",
      #   sbs %in% private_common ~ "private_common",
      #   sbs %in% shared ~ "shared")) %>%
      ggplot() +
      geom_point(aes(x=sbs, y=alpha), size=.5) +
      facet_grid(group~.) + ylim(0,1) +
      theme_bw() + theme(legend.position="bottom")
  )


plot_beta = function(beta)
  return(
    beta %>%
      tibble::rownames_to_column(var="sbs") %>%
      reshape2::melt(id="sbs", variable.name="context", value.name="beta") %>%
      ggplot() +
      geom_bar(aes(x=context, y=beta), stat="identity") +
      facet_grid(sbs~.)
  )
