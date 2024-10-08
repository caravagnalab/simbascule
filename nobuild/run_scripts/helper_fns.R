# generate_synthetic_datasets = function(shared,
#                                        private,
#                                        catalogue,
#                                        comb_matrix,
#                                        py,
#                                        out_path = NULL,
#                                        data_path = NULL,
#                                        seeds = 1:30,
#                                        mut_range = 10:8000,
#                                        input_catalogue = NULL,
#                                        reg_weight = 0.,
#                                        regularizer = "cosine",
#                                        CUDA = FALSE,
#                                        do.fits = FALSE,
#                                        verbose = FALSE,
#                                        new_model = TRUE) {
#   if (!dir.exists(out_path))
#     dir.create(out_path, recursive=T)
#
#   failed = file(paste0(out_path, "failed_runs.txt"), open="w")
#
#   reference_cat = catalogue[shared,]
#
#   for (i in 1:nrow(comb)) {
#
#     private_common = sample(private, comb$n_priv_comm[i])
#     tmp = setdiff(private, private_common)
#     private_rare = sample(tmp, comb$n_priv_rare[i])
#     denovo_cat = catalogue[c(private_common, private_rare),]
#
#     for (j in seeds) {
#
#       x = single_dataset(
#         N = comb$N_vals[i][[1]],
#         n_groups = comb$n_groups_vals[i][[1]],
#         samples_per_group = comb$samples_per_group[i][[1]],
#         reference_cat = reference_cat,
#         denovo_cat = denovo_cat,
#         private_sigs = list("rare" = private_rare, "common" = private_common),
#         private_fracs = list("rare" = 0.05, "common" = 0.1),
#         mut_range = mut_range,
#         seed = j,
#         out_path = data_path
#       )
#
#       min_k = max(1, nrow(reference_cat) + nrow(denovo_cat) - 5)
#       max_k = min_k + 10
#       k_list = min_k:max_k
#
#       idd = paste0("N", comb$N_vals[i][[1]], ".G", comb$n_groups_vals[i][[1]], ".s", j)
#
#       cat(paste0(idd, "\n"))
#
#       if (do.fits) {
#         fits = run_model(x = x$x[[1]],
#                   k = k_list,
#                   py = py,
#                   reference_catalogue = reference_cat,
#                   input_catalogue = input_catalogue,
#                   reg_weight = reg_weight,
#                   CUDA = CUDA,
#                   regularizer = regularizer,
#                   filtered_cat = TRUE,
#                   verbose = verbose,
#                   groups = x$groups[[1]] - 1,
#                   new_model = new_model,
#                   error_file = failed,
#                   idd = idd)
#
#         filename1 = paste0("fit.N", comb$N_vals[i][[1]], ".G",
#                            comb$n_groups_vals[i][[1]], ".s", seeds[j], ".Rds")
#
#         filename2 = paste0("fit.hier.N", comb$N_vals[i][[1]], ".G",
#                            comb$n_groups_vals[i][[1]], ".s", seeds[j], ".Rds")
#
#         save_fit(fits$fit1, out_path, filename1)
#         save_fit(fits$fit.hier, out_path, filename2)
#       }
#     }
#   }
#
#   close(failed)
# }
#
#
# save_fit = function(x.fit, path, filename) {
#   if (is.null(path)) return()
#
#   if (is.null(x.fit)) return()
#
#   if (!dir.exists(path))
#     dir.create(path, recursive=T)
#
#   if (filename %in% list.files(paste0(path)))
#     return()
#
#   saveRDS(x.fit, paste0(path, filename))
# }
#
#
# run_model = function(...,
#                      input_catalogue=NULL,
#                      filtered_cat=TRUE,
#                      groups=NULL,
#                      new_model=FALSE,
#                      error_file=NULL,
#                      idd="") {
#   msg1 = paste0("fit.", idd, "\n")
#   msg2 = paste0("fit.hier.", idd, "\n")
#
#   if (!new_model) {
#     x.fit = try_run(error_file,
#             expr = fit(..., groups=NULL, input_catalogue=input_catalogue),
#             msg = msg1)
#
#     x.fit.hier = try_run(error_file,
#             expr = fit(..., groups=groups, input_catalogue=input_catalogue),
#             msg = msg2)
#
#   } else {
#     x.fit = try_run(error_file,
#             expr = two_steps_inference(..., groups=NULL)$tot,
#             msg = msg1)
#
#     x.fit.hier = try_run(error_file,
#             expr = two_steps_inference(..., groups=groups)$tot,
#             msg = msg2)
#
#   }
#
#   return(list("fit1"=x.fit, "fit.hier"=x.fit.hier))
# }
#
#
# try_run = function(error_file, expr, msg) {
#
#   tryCatch(expr = expr,
#            error = function(e) {
#              writeLines(msg, error_file)
#              writeLines(paste(e))
# 	     writeLines(paste(reticulate::py_last_error()))
#              return(c(paste(e), paste(reticulate::py_last_error()) ) )
#            })
# }
#
#
#
# single_dataset = function(N, n_groups, samples_per_group,
#                           reference_cat, denovo_cat,
#                           private_sigs, private_fracs,
#                           cosine_limit, seed,
#                           reference_cosine=NULL, denovo_cosine=NULL,
#                           mut_range=10:8000, cohort_name="",
#                           out_path=NULL) {
#
#   groups = sample(1:n_groups, N, replace=T)
#   while (!all(lapply(1:n_groups, function(n) length(groups[groups==n]) %in% samples_per_group) %>% unlist()))
#     groups = sample(1:n_groups, N, replace=T)
#
#   idd = paste0("N", N, ".G", n_groups, ".s", seed)
#
#   if (cohort_name == "") out_name = paste0(out_path, "simul.", idd, ".Rds") else
#     out_name = paste0(out_path, "simul.", idd, ".", cohort_name, ".Rds")
#
#   if (!is.null(out_path) && paste0("simul.", idd, ".Rds") %in% list.files(paste0(out_path)))
#     return(readRDS(out_name))
#
#   x = generate.data(
#     reference_catalogue=reference_cat,
#     denovo_catalogue=denovo_cat,
#     reference_cosine=reference_cosine,
#     denovo_cosine=denovo_cosine,
#     targetX=-1,
#     inputX=NULL,
#     similarity_limit=cosine_limit,
#     groups=groups,
#     private_sigs=private_sigs,
#     private_fracs=private_fracs,
#     mut_range=mut_range,
#     seed=seed)
#
#   if (is.null(out_path)) return(x)
#
#   if (!dir.exists(out_path))
#     dir.create(out_path, recursive=T)
#
#   if (cohort_name == "")
#     saveRDS(x, paste0(out_path, "simul.", idd, ".Rds"))
#   else
#     saveRDS(x, paste0(out_path, "simul.", idd, ".", cohort_name, ".Rds"))
#
#   return(x)
# }



## Function to select "n_fixed" signatures from a reference catalogue (i.e. COSMIC)
## with a cosine similarity lower than "cosine_limit"
## It returns the subset reference
# select_fixed_sbs = function(referece_cat, n_fixed=4, cosine_limit=.5) {
#   ref_cosine.all = lsa::cosine(reference_cat %>% t())
#
#   while (TRUE) {
#     idxs = sample(1:nrow(ref_cosine.all), n_fixed)
#     sbss = rownames(ref_cosine.all)[idxs]
#     sbs_cosine = ref_cosine.all[sbss, sbss]
#     if (any(sbs_cosine != 1. & sbs_cosine > cosine_limit)) next
#
#     sbs_fixed = sbss
#     sbs_ref = reference_cat[sbs_fixed,]
#     break
#   }
#   return(sbs_ref)
# }
#
