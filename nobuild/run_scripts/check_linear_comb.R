fitname = "~/Dropbox/shared/2022. Basilica/simulations/fits_dn.clust.nonparametric.nonsparsity.noreg.old_hier.0208/fit_clust.N1000.G1.s1.6.Rds"
keep_sigs = c("SBS1","SBS5")
x = readRDS(fitname)
lc = filter_signatures_QP(sign1=get_denovo_signatures(x), sign2=COSMIC_filt, return_weights=FALSE)
new_ref = COSMIC_filt[unique(c(keep_sigs, unlist(lc))),]
