library(lsa)
library(dplyr)
library(ggplot2)
source("~/GitHub/simbasilica/nobuild/script_test/helper_fns.R")
py = reticulate::import_from_path("pybasilica","~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica/")
devtools::load_all()

path = "~/GitHub/simbasilica/nobuild/script_test/simulations/"
# source("./R/generate_data.R")
# source("./R/signatures.R")
# source("./R/input.R")
# source("./R/exposure.R")
# source("./R/theta.R")

make_plots = function(x, x.true=NULL, reconstructed=T, cls = NULL) {
  mm = plot_mutations(x, reconstructed = reconstructed)
  alp = plot_exposures(x, sort_by="SBS1", cls=cls)
  bet = plot_signatures(x, cls=cls)

  if (is.null(x.true)) return(patchwork::wrap_plots(bet + (mm/alp), guides="collect"))

  mm.true = plot_mutations(x.true, reconstructed = F)
  alpha.true = plot_exposures(x.true, sort_by="SBS1", cls=cls)

  return(patchwork::wrap_plots(bet + (mm/mm.true/alp/alpha.true), guides="collect"))
}

# Generate inputs ####

# cosmic = read.csv("./script_test/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t") %>%
#   tibble::column_to_rownames(var="Type") %>% t()

cosine_limit = .8
n_fixed = 2 # n of fixed signatures

shared = c("SBS1", "SBS40 SBS3 SBS5")
private_common = c("SBS17b", "SBS2", "SBS20")
private_rare = c("SBS8 SBS4")

denovo_catalogue = COSMIC_filt_merged[c(private_common, private_rare),]
# reference catalogue is all cosmic - the sigs selected as private common/rare and shared
# reference_cat = cosmic[!rownames(cosmic) %in% keep,]
# reference_cat = select_fixed_sbs(reference_cat, n_fixed=n_fixed, cosine_limit=cosine_limit)
reference_catalogue = COSMIC_filt_merged[shared,]
denovo_cosine = lsa::cosine(denovo_catalogue %>% t())
reference_cosine = lsa::cosine(reference_catalogue %>% t())


# First datasets - rare ####

# x = single_dataset(350, 2, 150:350, reference_catalogue, denovo_catalogue, reference_cosine, denovo_cosine,
#                    private_sigs=list("rare"=private_rare,"common"=private_common),
#                    private_fracs=list("rare"=0.05,"common"=0.3), cohort_name="rare",
#                    cosine_limit, seed=23, out_path="./nobuild/script_test/simulations/")

x.simul = readRDS(paste0(path, "simul.N350.G2.s23.rare.Rds"))
xx = create_basilica_obj_simul(x.simul)
make_plots(xx, reconstructed = F)

## Old model fit, no groups
x.fit.good = x.fit
x.fit = fit(x.simul$x[[1]], k=0:7, py=py,
            input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
            reference_catalogue=COSMIC_filt_merged,
            reg_bic=TRUE,
            filtered_cat=TRUE,
            regularizer = "cosine")


# Old model fit, with groups
x.fit.g = fit(x.simul$x[[1]], k=0:7, py=py,
              input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
              reference_catalogue=COSMIC_filt_merged,
              reg_bic=TRUE, filtered_cat=TRUE, groups=x.simul$groups[[1]]-1,
              regularizer = "cosine")
x.fit.g$groups = x.simul$groups[[1]]-1


sigs = c(rownames(get_signatures(xx)),
         rownames(get_signatures(x.fit)),
         rownames(get_signatures(x.fit.g))) %>% unique()
cls = gen_palette(length(sigs)) %>% setNames(sigs)



true.plots = make_plots(xx, reconstructed=F, cls=cls) & theme(legend.position="bottom")
x1.plots = make_plots(x.fit, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")
x2.plots = make_plots(x.fit.g, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")


pdf("~/GitHub/hBasilica/results/fit_simulated/report.simul_tests.norare.pdf", height = 8, width = 14)
true.plots & patchwork::plot_annotation(title="True data")
x1.plots & patchwork::plot_annotation(title="Fit 1")
x2.plots & patchwork::plot_annotation(title="Fit 2 - hierarchical")
# patchwork::wrap_plots(plot_signatures(x1, signatures="SBS20", cls=cls) /
#                         plot_signatures(x2,
#                                         signatures=setdiff(rownames(get_signatures(x2)),
#                                                            rownames(get_signatures(x1))),
#                                         cls=cls))
dev.off()






## new model
x.fit.new = two_steps_inference(x.simul$x[[1]],
                                k=0:7,
                                input_catalogue=COSMIC_filt_merged,
                                enforce_sparsity1=TRUE,
                                enforce_sparsity2=FALSE,
                                py=py)
x.new = x.fit.new$tot
make_plots(x.new, xx, reconstructed = T)

# saveRDS(x.new, "~/GitHub/simbasilica/nobuild/test_new1.Rds")
#
#
#
# x.fit.new2 = two_steps_inference(x.simul$x[[1]],
#                                 k=0:7,
#                                 input_catalogue=COSMIC_filt_merged[rownames(COSMIC_filt_merged)!="SBS20",],
#                                 enforce_sparsity1=TRUE,
#                                 enforce_sparsity2=FALSE,
#                                 py=py)
# x.new2 = x.fit.new2$tot
# make_plots(x.new2, xx, reconstructed = T)
#
# patchwork::wrap_plots(plot_signatures(x.new2, signatures = c("SBS14","SBS44")) /
#                         plot_signatures(x.new, signatures = "SBS20"))



# Second datasets - no rare #####

cosine_limit = .8
n_fixed = 2

shared = c("SBS1", "SBS40 SBS3 SBS5")
private_common = c("SBS8 SBS4", "SBS2", "SBS20", "SBS17b")
private_rare = c()

reference_catalogue = COSMIC_filt_merged[shared,]
denovo_catalogue = COSMIC_filt_merged[c(private_common, private_rare),]

denovo_cosine = lsa::cosine(denovo_catalogue %>% t())
reference_cosine = lsa::cosine(reference_catalogue %>% t())

# x.simul2 = single_dataset(350, 2, 150:350, reference_catalogue, denovo_catalogue,
#                          reference_cosine, denovo_cosine,
#                          private_sigs=list("rare"=private_rare,"common"=private_common),
#                          private_fracs=list("rare"=0.05,"common"=0.3), cohort_name="norare",
#                          cosine_limit, seed=23, out_path="./nobuild/script_test/simulations/")
x.simul2 = readRDS("~/GitHub/simbasilica/nobuild/script_test/simulations/simul.N350.G2.s23.norare.Rds")
xx2 = create_basilica_obj_simul(x.simul2)
make_plots(xx2, reconstructed = F)

## Old model fit, no groups
x.fit = fit(x.simul2$x[[1]], k=0:7, py=py,
            input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
            reference_catalogue=COSMIC_filt_merged,
            reg_bic=TRUE,
            filtered_cat=TRUE,
            regularizer = "cosine")


# Old model fit, with groups
x.fit.g = fit(x.simul2$x[[1]], k=0:7, py=py,
              input_catalogue=COSMIC_filt_merged[c("SBS1","SBS40 SBS3 SBS5"),],
              reference_catalogue=COSMIC_filt_merged,
              reg_bic=TRUE, filtered_cat=TRUE, groups=x.simul2$groups[[1]]-1,
              regularizer = "cosine")
x.fit.g$groups = x.simul$groups[[1]]-1


sigs = c(rownames(get_signatures(xx)),
         rownames(get_signatures(x.fit)),
         rownames(get_signatures(x.fit.g))) %>% unique()
cls = gen_palette(length(sigs)) %>% setNames(sigs)



true.plots = make_plots(xx, reconstructed=F, cls=cls) & theme(legend.position="bottom")
x1.plots = make_plots(x.fit, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")
x2.plots = make_plots(x.fit.g, xx, reconstructed=T, cls=cls) & theme(legend.position="bottom")


pdf("~/GitHub/hBasilica/results/fit_simulated/report.simul_tests.rare.pdf", height = 8, width = 14)
true.plots & patchwork::plot_annotation(title="True data")
x1.plots & patchwork::plot_annotation(title="Fit 1")
x2.plots & patchwork::plot_annotation(title="Fit 2 - hierarchical")
# patchwork::wrap_plots(plot_signatures(x1, signatures="SBS20", cls=cls) /
#                         plot_signatures(x2,
#                                         signatures=setdiff(rownames(get_signatures(x2)),
#                                                            rownames(get_signatures(x1))),
#                                         cls=cls))
dev.off()




# x.fit2 = two_steps_inference(x.simul2$x[[1]],
#                              k=0:7,
#                              input_catalogue=COSMIC_filt_merged,
#                              enforce_sparsity1=TRUE,
#                              enforce_sparsity2=FALSE,
#                              py=py)
# x.new2 = x.fit2$tot
# saveRDS(x.new2, "~/GitHub/simbasilica/nobuild/test_new1.Rds")




x.new2 %>% plot_exposures()
x.new2 %>% plot_signatures()
x.new2 %>% plot_similarity_reference()
x.new2 %>% plot_mutations()




## Old #########
x$x[[1]] %>%
  dplyr::mutate(group=x$groups[[1]]) %>%
  reshape2::melt(id="group", variable.name="context", value.name="n_muts") %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(n_muts_d=n_muts/sum(n_muts)) %>%
  ggplot() +
  geom_bar(aes(x=context, y=n_muts_d), stat="identity") +
  facet_grid(group~.)

x$exp_fixed[[1]] %>% as.data.frame() %>%
  plot_beta()

x$exp_denovo[[1]] %>% as.data.frame() %>%
  plot_beta()

x$exp_exposure[[1]] %>% dplyr::mutate(group=x$groups[[1]]) %>%
  plot_alpha()







## TODO Sample 50/100 datasets #####

N_vals = c(100, 500, 1000)
n_groups_vals = c(4, 8, 10)
samples_per_group = 15:100




## Read simulation #####

x = readRDS("./script_test/simulations/simul.N100.G5.s23.Rds")

alpha = x$exp_exposure[[1]]

counts = x$x[[1]] %>%
  dplyr::mutate(groups=x$groups[[1]])

write.csv(counts, "./script_test/simul.N100.G5.s23_counts.csv", row.names=F)
write.csv(alpha, "./script_test/simul.N100.G5.s23_alpha.csv", row.names=F)



