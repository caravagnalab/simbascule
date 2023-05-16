library(lsa)
library(dplyr)
library(ggplot2)
source("~/GitHub/simbasilica/nobuild/script_test/helper_fns.R")
py = reticulate::import_from_path("pybasilica","~/GitHub/pybasilica/")
devtools::load_all("~/GitHub/basilica/")
devtools::load_all()
# source("./R/generate_data.R")
# source("./R/signatures.R")
# source("./R/input.R")
# source("./R/exposure.R")
# source("./R/theta.R")

# Generate inputs ####

# cosmic = read.csv("./script_test/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t") %>%
#   tibble::column_to_rownames(var="Type") %>% t()

cosmic = COSMIC_filt_merged
cosine_limit = .8
n_fixed = 2 # n of fixed signatures

## De Novo signatures ####

# 1 -> sbs shared by ALL samples
# shared = c("SBS1","SBS5")
shared = c("SBS1", "SBS40 SBS3 SBS5")

# 2 -> sbs private to a group and with high frequency (common)
## apobec -> SBS2; tobacco smoking -> SBS4; chemotherapy -> SBS25
private_common = c("SBS17b", "SBS2", "SBS20")

# 3 -> sbs private to a group and with low frequency (rare)
## defective DNA mmr -> SBS26; defective DNA mmr breast cancer -> SBS44; apobec -> SBS13;
## somatic hypermutation in lymphoid cells -> SBS9; exposure to e. coli in CRC -> SBS88
private_rare = c("SBS8 SBS4")

# keep = lapply(rownames(cosmic), function(sname) {
#   snames = strsplit(sname, " ")[[1]]
#   if (any(snames %in% c(private_common, private_rare))) return(sname)
# }) %>% unlist()

denovo_catalogue = cosmic[c(private_common, private_rare),]

## Refrence signatures ####

# reference catalogue is all cosmic - the sigs selected as private common/rare and shared
# reference_cat = cosmic[!rownames(cosmic) %in% keep,]
# reference_cat = select_fixed_sbs(reference_cat, n_fixed=n_fixed, cosine_limit=cosine_limit)
reference_catalogue = cosmic[shared,]


## Cosine similarity ####
denovo_cosine = lsa::cosine(denovo_catalogue %>% t())
reference_cosine = lsa::cosine(reference_catalogue %>% t())


## Generate dataset ####

# x = single_dataset(350, 2, 150:350, reference_catalogue, denovo_catalogue, reference_cosine, denovo_cosine,
#                    private_sigs=list("rare"=private_rare,"common"=private_common),
#                    private_fracs=list("rare"=0.05,"common"=0.3),
#                    cosine_limit, seed=23, out_path="./nobuild/script_test/simulations/")

x.simul = readRDS("./nobuild/script_test/simulations/simul.N350.G2.s23.Rds")
xx = create_basilica_obj_simul(x.simul)
# plot_mutations(xx)
# plot_exposures(xx, sort_by="SBS1")
# plot_signatures(xx)

x.fit = fit(x.simul$x[[1]], k=0:3, py=py, input_catalogue=NULL,
            reference_catalogue=COSMIC_filt_merged,
            reg_bic=TRUE, filtered_cat=TRUE, regularizer = "KL")
plot_mutations(x.fit)
plot_exposures(x.fit, sort_by="SBS1")
plot_signatures(x.fit)


x.fit.g = fit(x$x[[1]], k=1:7, py=py, input_catalogue=NULL,
            reference_catalogue=COSMIC_filt_merged,
            reg_bic=TRUE, filtered_cat=TRUE, groups=x$groups[[1]]-1,
            regularizer="KL")
x.fit.g$groups = x$groups[[1]]-1
plot_mutations(x.fit.g)
plot_exposures(x.fit.g, sort_by="SBS1")
plot_signatures(x.fit.g)


two_steps_inference(x.simul$x[[1]],
                    k=0:7,
                    input_catalogue = COSMIC_filt_merged[c("SBS1", "SBS40 SBS3 SBS5", "SBS2", "SBS17b", "SBS20"),],
                    enforce_sparsity1=TRUE,
                    enforce_sparsity2=FALSE,
                    py=py)

x.fit.new = two_steps_inference(x$x[[1]],
                                k=0:7,
                                enforce_sparsity1=TRUE,
                                enforce_sparsity2=FALSE,
                                py=py)
x.new = x.fit.new$tot
plot_mutations(x.new)
plot_exposures(x.new, sort_by="SBS1")
plot_signatures(x.new)


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



