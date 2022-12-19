# test synthetic data

#M_path <- "/home/azad/Documents/thesis/pybasilica/pybasilica/data/real/data_sigphylo.csv"
#input_catalogue_path <- "/home/azad/Documents/thesis/pybasilica/pybasilica/data/real/beta_aging.csv"
ref_path <- "/home/azad/Documents/thesis/pybasilica/pybasilica/data/cosmic/cosmic_catalogue.csv"
#M <- read.table(M_path, sep = ",", header = TRUE, check.names = FALSE)
#reference_catalogue <- read.table(reference_catalogue_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
#input_catalogue <- read.table(input_catalogue_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
#setwd("/home/azad/Documents/thesis/basilica")


s <- split.reference(reference_path=ref_path, ratio=0.7, seed=45)
reference_catalogue <- s$reference
denovo_catalogue <- s$denovo

theta <- generate.theta(mut_range=100:50, num_samples=3, seed=34)

reference_cosine <- cosine.matrix(reference_catalogue, reference_catalogue)
denovo_cosine <- cosine.matrix(denovo_catalogue, denovo_catalogue)
signatures <- generate.signatures(
  reference_catalogue=reference_catalogue,
  denovo_catalogue=denovo_catalogue,
  reference_cosine=reference_cosine, # cosine similarity matrix of reference signatures (SBS1 excluded)
  denovo_cosine=denovo_cosine,       # cosine similarity matrix of denovo signatures
  complexity=c(4,3),
  similarity_limit=0.5,
  seed=88
)

alpha <- generate.exposure(beta=signatures, groups=rep(1,10), seed=NULL)

