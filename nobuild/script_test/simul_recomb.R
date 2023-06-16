library(lsa)
library(dplyr)
library(ggplot2)
source("~/GitHub/simbasilica/script_test/helper_fns.R")

devtools::load_all("~/GitHub/basilica/")

# cosmic = read.csv("./script_test/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t") %>%
#   tibble::column_to_rownames(var="Type") %>% t()

# cosmic = COSMIC_filtered
#
# # cosine_limit = .8
# # n_fixed = 4 # n of fixed signatures
#
# ## De Novo signatures ####
#
# shared = c("SBS1","SBS5","SBS17b")
#
# private_common = c("SBS4", "SBS6", "SBS10a")
#
# private_rare = c()
#
# denovo_cat = cosmic[c(private_common, private_rare),]
#
# ## Refrence signatures ####
# reference_cat = cosmic[shared,]


shared = c("SBS1","SBS5", "SBS17b")
sbs.names = c("SBS1","SBS5",
            "SBS2","SBS6","SBS10a","SBS10b","SBS21",  # CRC
            "SBS2","SBS3","SBS4","SBS8","SBS13","SBS15","SBS17a", "SBS17b", "SBS18",  # Lung
              "SBS3","SBS8","SBS10","SBS13","SBS17","SBS18","SBS20","SBS26")  # Breast
sbs.names = unique(sbs.names)
catalogue = COSMIC_filtered[intersect(sbs.names, rownames(COSMIC_filtered)),]
reference_cat = catalogue[shared,]
# 2 -> sbs private to a group and with high frequency (common)
## apobec -> SBS2; tobacco smoking -> SBS4; chemotherapy -> SBS25
private_common = setdiff(rownames(catalogue), shared)
# 3 -> sbs private to a group and with low frequency (rare)
## defective DNA mmr -> SBS26; defective DNA mmr breast cancer -> SBS44; apobec -> SBS13;
## somatic hypermutation in lymphoid cells -> SBS9; exposure to e. coli in CRC -> SBS88
private_rare = c()
denovo_cat = catalogue[c(private_common, private_rare),]


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

# x = single_dataset(1000, 5, 100:1000, reference_cat, denovo_cat, reference_cosine, denovo_cosine,
#                    private_sigs=list("rare"=private_rare,"common"=private_common),
#                    private_fracs=list("rare"=0.05,"common"=0.3), cosine_limit=cosine_limit,
#                    seed=23, out_path="./script_test/simulations/")
# plot_alpha(x)


## Load Data ####
x = readRDS("~/GitHub/simbasilica/script_test/simulations/simul.N100.G2.s23.Rds"); N=100; G=2


dn_c = ggsci::pal_nejm()(x$exp_denovo[[1]] %>% nrow)
names(dn_c) = x$exp_denovo[[1]] %>% rownames()
sc_c = ggsci::pal_simpsons()(x$exp_fixed[[1]] %>% nrow)
names(sc_c) = x$exp_fixed[[1]] %>% rownames()
cls = c(sc_c, dn_c)

my_plot_signatures(x, cls=cls)
my_plot_exposure(x, cls=cls)
plot_simulated_data(x)

reticulate::use_condaenv("basilica-env"); py = NULL
py = reticulate::import_from_path(module="pybasilica", path="~/GitHub/pybasilica/")


## Run the inference ####

x.fit.noreg = fit(x=x$x[[1]], k=0:7, py=py, # groups = x$groups[[1]]-1,
                  reference_catalogue=basilica::COSMIC_catalogue[c(shared, private_common),],
                  # reference_catalogue=basilica::COSMIC_catalogue[c("SBS1","SBS5"),],
                  input_catalogue=basilica::COSMIC_catalogue[c("SBS1","SBS5"),],
                  reg_weight=0.)

# saveRDS(x.fit.noreg, "./script_test/simulations/fit_noreg_nogroups_subcosmic.N100.G2.s23.Rds")

x.fit.noreg %>% plot_signatures()
x.fit.noreg %>% plot_exposure()
x.fit.noreg %>% plot_similarity_reference(reference=basilica::COSMIC_catalogue[c(shared, private_common),])


## with regularization
x.fit.reg = fit(x=x$x[[1]], k=0:7, py=py, # groups = x$groups[[1]]-1,
                reference_catalogue=basilica::COSMIC_catalogue[c(shared, private_common),],
                # reference_catalogue=basilica::COSMIC_catalogue[c("SBS1","SBS5"),],
                input_catalogue=basilica::COSMIC_catalogue[c("SBS1","SBS5"),],
                reg_weight=1.)

x.fit.reg %>% plot_signatures()
x.fit.reg %>% plot_exposure()
x.fit.reg %>% plot_similarity_reference(reference=basilica::COSMIC_catalogue[c(shared, private_common),])

saveRDS(x.fit.reg, "./script_test/simulations/fit_reg_groups_subcosmic.N100.G2.s23.Rds")


