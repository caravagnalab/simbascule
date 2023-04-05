library(lsa)
library(dplyr)
library(ggplot2)
source("./script_test/helper_fns.R")


cosmic = read.csv("./script_test/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t") %>%
  tibble::column_to_rownames(var="Type") %>% t()

# cosine_limit = .8
# n_fixed = 4 # n of fixed signatures

## De Novo signatures ####

shared = c("SBS1","SBS5","SBS17b")

private_common = c("SBS4", "SBS6", "SBS10a")

private_rare = c()

denovo_cat = cosmic[c(private_common, private_rare),]

## Refrence signatures ####
reference_cat = cosmic[shared,]



## Cosine similarity ####
denovo_cosine = lsa::cosine(denovo_cat %>% t())
ref_cosine = lsa::cosine(reference_cat %>% t())


## Generate dataset ####

x = single_dataset(100, 2, 50:100, reference_cat, denovo_cat, reference_cosine, denovo_cosine,
                   private_sigs=list("rare"=private_rare,"common"=private_common),
                   private_fracs=list("rare"=0.05,"common"=0.3), cosine_limit=cosine_limit,
                   seed=23, out_path="./script_test/simulations/")

plot_beta(x)
plot_alpha(x)
plot_muts(x)

if (Sys.getenv("GITHUB_PATH") == "") path=paste0("~/dati_elenab/signatures/") else path=Sys.getenv("GITHUB_PATH")
py_path = paste0(path, "pybasilica")
py = reticulate::import_from_path(module="pybasilica", path=py_path)

x.fit.noreg = basilica::fit(x=x$x[[1]], k=0:7, py=py, # groups = x$groups[[1]]-1,
                            # reference_catalogue=basilica::COSMIC_catalogue[c(shared, private_common),],
                            reference_catalogue=basilica::COSMIC_catalogue[shared,],
                            input_catalogue=basilica::COSMIC_catalogue[c("SBS1","SBS5"),])



saveRDS(x.fit.noreg, "./script_test/simulations/fit_noreg_nogroups_no_private.Rds")

x.fit.noreg %>% plot_signatures()
x.fit.noreg %>% plot_exposure()
x.fit.noreg %>% plot_similarity_reference(reference = cosmic[c(shared,private_common, private_rare),])

x.fit.reg = basilica::fit(x=x$x[[1]], k=1:7, py=py, groups = x$groups[[1]]-1,
                          reference_catalogue=basilica::COSMIC_catalogue[c(shared, private_common),],
                          input_catalogue=basilica::COSMIC_catalogue[c("SBS1","SBS5"),],
                          reg_weight=1.)

x.fit.reg %>% basilica::plot_signatures()
x.fit.reg %>% basilica::plot_exposure()
x.fit.reg %>% basilica::plot_similarity_reference()

saveRDS(x.fit.reg, "./script_test/simulations/fit_reg_groups_subcosmic.N100.G5.s23.Rds")


