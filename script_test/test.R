library(lsa)
library(dplyr)
library(ggplot2)
source("./script_test/helper_fns.R")
source("./R/generate_data.R")
source("./R/signatures.R")
source("./R/input.R")
source("./R/exposure.R")
source("./R/theta.R")

# Generate inputs ####

cosmic = read.csv("./script_test/COSMIC_v3.3.1_SBS_GRCh38.txt", sep="\t") %>%
  tibble::column_to_rownames(var="Type") %>% t()

cosine_limit = .8
n_fixed = 4 # n of fixed signatures

## De Novo signatures ####

# 1 -> sbs shared by ALL samples
shared = c("SBS1","SBS5")

# 2 -> sbs private to a group and with high frequency (common)
## apobec -> SBS2; tobacco smoking -> SBS4; chemotherapy -> SBS25
private_common = c("SBS2", "SBS4", "SBS25")

# 3 -> sbs private to a group and with low frequency (rare)
## defective DNA mmr -> SBS26; defective DNA mmr breast cancer -> SBS44; apobec -> SBS13;
## somatic hypermutation in lymphoid cells -> SBS9; exposure to e. coli in CRC -> SBS88
private_rare = c("SBS44", "SBS9", "SBS13", "SBS88")

denovo_cat = cosmic[c(shared, private_common, private_rare),]

## Refrence signatures ####

# reference catalogue is all cosmic - the sigs selected as private common/rare and shared
reference_cat = cosmic[!rownames(cosmic) %in% c(shared, private_common, private_rare),]
reference_cat = select_fixed_sbs(reference_cat, n_fixed=n_fixed, cosine_limit=cosine_limit)



## Cosine similarity ####
denovo_cosine = lsa::cosine(denovo_cat %>% t())
reference_cosine = lsa::cosine(reference_cat %>% t())


## Generate dataset ####

x = single_dataset(1000, 5, 150:1000, reference_cat, denovo_cat, reference_cosine, denovo_cosine,
                   private_sigs=list("rare"=private_rare,"common"=private_common),
                   private_fracs=list("rare"=0.05,"common"=0.3),
                   cosine_limit, seed=23, out_path="./script_test/simulations/")


saveRDS(x, file = "./script_test/simulations/simul.N1000.G5.s23.Rds")


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



