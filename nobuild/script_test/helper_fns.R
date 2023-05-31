generate_synthetic_datasets = function(shared,
                                       private,
                                       catalogue,
                                       comb_matrix,
                                       py,
                                       out_path = NULL,
                                       data_path = NULL,
                                       seeds = 1:30,
                                       mut_range = 10:8000,
                                       input_catalogue = NULL,
                                       reg_weight = 0.,
                                       regularizer = "cosine",
                                       CUDA = FALSE,
                                       do.fits = FALSE,
                                       verbose = FALSE,
                                       new_model = TRUE) {

  failed = file(paste0(out_path, "failed_runs.txt"), open="w")

  reference_cat = catalogue[shared,]

  for (i in 1:nrow(comb)) {

    private_common = sample(private, comb$n_priv_comm[i])
    tmp = setdiff(private, private_common)
    private_rare = sample(tmp, comb$n_priv_rare[i])
    denovo_cat = catalogue[c(private_common, private_rare),]

    for (j in seeds) {

      x = single_dataset(
        N = comb$N_vals[i][[1]],
        n_groups = comb$n_groups_vals[i][[1]],
        samples_per_group = comb$samples_per_group[i][[1]],
        reference_cat = reference_cat,
        denovo_cat = denovo_cat,
        private_sigs = list("rare" = private_rare, "common" = private_common),
        private_fracs = list("rare" = 0.05, "common" = 0.1),
        mut_range = mut_range,
        seed = j,
        out_path = data_path
      )

      min_k = max(1, nrow(reference_cat) + nrow(denovo_cat) - 5)
      max_k = min_k + 10
      k_list = min_k:max_k

      idd = paste0("N", comb$N_vals[i][[1]], ".G", comb$n_groups_vals[i][[1]], ".s", j)

      cat(paste0(idd, "\n"))

      if (do.fits) {
        fits = run_model(x = x$x[[1]],
                  k = k_list,
                  py = py,
                  reference_catalogue = reference_cat,
                  input_catalogue = input_catalogue,
                  reg_weight = reg_weight,
                  CUDA = CUDA,
                  regularizer = regularizer,
                  filtered_cat = TRUE,
                  verbose = verbose,
                  groups = x$groups[[1]] - 1,
                  new_model = new_model,
                  error_file = failed,
                  idd = idd)

        filename1 = paste0("fit.N", comb$N_vals[i][[1]], ".G",
                           comb$n_groups_vals[i][[1]], ".s", seeds[j], ".Rds")

        filename2 = paste0("fit.hier.N", comb$N_vals[i][[1]], ".G",
                           comb$n_groups_vals[i][[1]], ".s", seeds[j], ".Rds")

        save_fit(fits$fit1, out_path, filename1)
        save_fit(fits$fit.hier, out_path, filename2)
      }
    }
  }

  close(failed)
}


save_fit = function(x.fit, path, filename) {
  if (is.null(path)) return()

  if (!dir.exists(path))
    dir.create(path, recursive=T)

  if (filename %in% list.files(paste0(path)))
    return()

  saveRDS(x.fit, paste0(path, filename))
}


run_model = function(...,
                     input_catalogue=NULL,
                     filtered_cat=TRUE,
                     groups=NULL,
                     new_model=FALSE,
                     error_file=NULL,
                     idd="") {
  msg1 = paste0("fit.", idd, "\n")
  msg2 = paste0("fit.hier.", idd, "\n")

  if (!new_model) {
    x.fit = try_run(error_file,
            expr = fit(..., groups=NULL, input_catalogue=input_catalogue),
            msg = msg1)

    x.fit.hier = try_run(error_file,
            expr = fit(..., groups=groups, input_catalogue=input_catalogue),
            msg = msg2)

  } else {
    x.fit = try_run(error_file,
            expr = two_steps_inference(..., groups=NULL)$tot,
            msg = msg1)

    x.fit.hier = try_run(error_file,
            expr = two_steps_inference(..., groups=groups)$tot,
            msg = msg2)

  }

  return(list("fit1"=x.fit, "fit.hier"=x.fit.hier))
}


try_run = function(error_file, expr, msg) {

  tryCatch(expr = expr,
           error = function(e) {
             writeLines(msg, error_file)
             writeLines(paste(e))
             return(NULL)
           })
}



single_dataset = function(N, n_groups, samples_per_group,
                          reference_cat, denovo_cat,
                          private_sigs, private_fracs,
                          cosine_limit, seed,
                          reference_cosine=NULL, denovo_cosine=NULL,
                          mut_range=10:8000, cohort_name="",
                          out_path=NULL) {

  groups = sample(1:n_groups, N, replace=T)
  while (!all(lapply(1:n_groups, function(n) length(groups[groups==n]) %in% samples_per_group) %>% unlist()))
    groups = sample(1:n_groups, N, replace=T)

  idd = paste0("N", N, ".G", n_groups, ".s", seed)

  if (cohort_name == "") out_name = paste0(out_path, "simul.", idd, ".Rds") else
    out_name = paste0(out_path, "simul.", idd, ".", cohort_name, ".Rds")

  if (!is.null(out_path) && paste0("simul.", idd, ".Rds") %in% list.files(paste0(out_path)))
    return(readRDS(out_name))

  x = generate.data(
    reference_catalogue=reference_cat,
    denovo_catalogue=denovo_cat,
    reference_cosine=reference_cosine,
    denovo_cosine=denovo_cosine,
    targetX=-1,
    inputX=NULL,
    similarity_limit=cosine_limit,
    groups=groups,
    private_sigs=private_sigs,
    private_fracs=private_fracs,
    mut_range=mut_range,
    seed=seed)

  if (is.null(out_path)) return(x)

  if (!dir.exists(out_path))
    dir.create(out_path, recursive=T)

  if (cohort_name == "")
    saveRDS(x, paste0(out_path, "simul.", idd, ".Rds"))
  else
    saveRDS(x, paste0(out_path, "simul.", idd, ".", cohort_name, ".Rds"))

  return(x)
}



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




## Visualization functions ####
my_plot_exposure = function(x, cls=RColorBrewer::brewer.pal(n=9, name="Set1")) {
  alpha = x$exp_exposure[[1]]
  alpha$groups = x$groups[[1]]
  n = ncol(alpha)
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cls = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  return(
    alpha %>% as.data.frame() %>%
      dplyr::mutate(sample=paste0(1:nrow(alpha))) %>%
      reshape2::melt(id=c("sample","groups")) %>%
      dplyr::rename(Signature=variable) %>%
      ggplot(aes(x=sample, y=value, fill=Signature)) +
      geom_bar(stat="identity") +
      scale_fill_manual(values=cls) +
      labs(title="Exposure") +
      theme(axis.ticks.x=element_blank(),axis.text.x=element_blank()) + ylab("") +
      facet_grid(~groups, scales="free_x")
  )
}


my_plot_signatures = function(x, cls=RColorBrewer::brewer.pal(n=9, name="Set1")) {
  beta = rbind(x$exp_fixed[[1]], x$exp_denovo[[1]])
  n = ncol(beta)
  qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
  cls = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

  return(
    beta %>%
      as.data.frame() %>%
      tibble::rownames_to_column(var="sbs") %>%
      reshape2::melt(id="sbs") %>%
      dplyr::rename(Var1=sbs, Var2=variable) %>%
      dplyr::mutate(substitution=paste0(substr(start=3, stop=3, Var2),">",substr(start=5, stop=5, Var2)),
             context=paste0(substr(start=1, stop=1, Var2), "_", substr(start=7, stop=7, Var2))) %>%
      ggplot() +
      geom_bar(aes(value, x=context, fill=Var1), stat="identity") +
      facet_grid(Var1 ~ substitution, scales="free") +
      theme_bw() +
      theme(axis.ticks.x=element_blank(),axis.text.x=element_blank()) +
      scale_fill_manual(values=cls) +
      guides(fill="none")  +
      labs(x="Context", y="", title="Signatures")
  )
}


plot_simulated_data = function(x, cls=RColorBrewer::brewer.pal(n=9, name="Set1")) {
  x$x[[1]] %>%
    dplyr::mutate(group=paste0(x$groups[[1]])) %>%
    reshape2::melt(id="group") %>%
    dplyr::mutate(substitution=substr(start=3, stop=5, variable),
                  context = paste0(substr(start=1, stop=1, variable), "_", substr(start=7, stop=7, variable))) %>%
    dplyr::group_by(group) %>%
    dplyr::mutate(n_muts=value/sum(value)) %>%

    ggplot() +
    geom_bar(aes(y=n_muts, x=context, fill=group), stat="identity") +
    theme_bw()  +
    facet_grid(group~ substitution, scales="free") +
    scale_fill_manual(values=cls) +
    theme(strip.text.y=element_text(angle=0),
          axis.text.x=element_text(angle=90, size=4),
          legend.position="none") +
    labs(x="Context", y="Mutation count") + labs(title="Data")
}



## OLD PLOTS ####

# plot_muts = function(x) {
#   return(
#     x$x[[1]] %>%
#       dplyr::mutate(group=x$groups[[1]]) %>%
#       reshape2::melt(id="group", variable.name="context", value.name="n_muts") %>%
#       dplyr::group_by(group) %>%
#       dplyr::mutate(n_muts_d=n_muts/sum(n_muts)) %>%
#       ggplot() +
#       geom_bar(aes(x=context, y=n_muts_d), stat="identity") +
#       facet_grid(group~.)
#   )
# }

plot_alpha = function(x) {
  alpha = x$exp_exposure[[1]]
  if (!"group" %in% colnames(x)) alpha$group = x$groups[[1]]
  return(
    alpha %>%
      as.data.frame() %>%
      reshape2::melt(id="group", variable.name="sbs", value.name="alpha") %>%
      ggplot() +
      geom_point(aes(x=sbs, y=alpha), size=.5) +
      facet_grid(group~.) + ylim(0,1) +
      theme_bw() + theme(legend.position="bottom")
  )
}

# plot_beta = function(x) {
#   beta = rbind(x$exp_fixed[[1]],
#                x$exp_denovo[[1]])
#   return(
#     beta %>%
#       as.data.frame() %>%
#       tibble::rownames_to_column(var="sbs") %>%
#       reshape2::melt(id="sbs", variable.name="context", value.name="beta") %>%
#       ggplot() +
#       geom_bar(aes(x=context, y=beta), stat="identity") +
#       facet_grid(sbs~.)
#   )
# }
