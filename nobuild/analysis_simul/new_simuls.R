
gen_data = function(N, G, sbs_catalogue, alpha_range, alpha_sigma, seed,
                    shared_sbs=c("SBS1","SBS5")) {
  sbs = sample(setdiff(rownames(COSMIC_filt), shared_sbs), G + floor(G/2), replace=F)
  private_shared_sbs = sample(sbs, floor(G/2))
  private_sbs = setdiff(sbs, private_shared_sbs)

  return(
    generate_simulation_dataset(G=G, N=N, alpha_sigma=alpha_sigma,
                                sbs_catalogue=COSMIC_filt, private_sbs=private_sbs,
                                private_shared_sbs=private_shared_sbs,
                                shared_sbs=shared_sbs, alpha_range = alpha_range,
                                py=py, seed=seed)
  )
}


sbs_shared = c("SBS1","SBS5")
sbs_private = c("SBS4","SBS13","SBS10b","SBS7c","SBS7d",
                "SBS11","SBS90","SBS31","SBS10a","SBS22")

m = lsa::cosine(t(COSMIC_filt[c(sbs_shared, sbs_private),]))
diag(m) = 0
pheatmap::pheatmap(m)

data1 = gen_data(N=150, G=6, sbs_catalogue=COSMIC_filt[sbs_private,], alpha_range=c(.15,0.2), alpha_sigma=0.1, seed=1)
data2 = gen_data(N=100, G=3, sbs_catalogue=COSMIC_filt, alpha_range=c(.15,0.2), alpha_sigma=0.1, seed=1)

data1$alpha_plot
data2$alpha_plot


test = two_steps_inference(x=data1$counts[[1]] %>% dplyr::select(-sample, -groupid),
                           k=4:10, clusters=1:10,
                           enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                           lr=0.005, n_steps=2000, py=py,
                           hyperparameters=list("alpha_sigma"=0.1),
                           seed_list=c(10), reg_weight=1., regularizer="KL",
                           new_hier=FALSE, do_initial_fit=FALSE, verbose = T, regul_fixed = FALSE,
                           save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)

test.cosine = two_steps_inference(x=data1$counts[[1]] %>% dplyr::select(-sample, -groupid),
                           k=4:10, clusters=1:10,
                           enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                           lr=0.005, n_steps=2000, py=py,
                           hyperparameters=list("alpha_sigma"=0.1),
                           seed_list=c(10), reg_weight=1., regularizer="cosine",
                           new_hier=FALSE, do_initial_fit=TRUE, verbose = T, regul_fixed = FALSE,
                           save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)

test2 = two_steps_inference(x=data1$counts[[1]] %>% dplyr::select(-sample, -groupid),
                                  k=4:10, clusters=1:10,
                                  enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                                  lr=0.005, n_steps=2000, py=py,
                                  hyperparameters=list("alpha_sigma"=0.05),
                                  seed_list=c(10, 33, 92), reg_weight=0., regularizer="cosine",
                                  new_hier=FALSE, do_initial_fit=TRUE, verbose = T, regul_fixed = FALSE,
                                  save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)

test3 = two_steps_inference(x=data1$counts[[1]] %>% dplyr::select(-sample, -groupid),
                            k=4:10, clusters=1:10,
                            enforce_sparsity2=F, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                            lr=0.005, n_steps=2000, py=py,
                            hyperparameters=list("alpha_sigma"=0.05),
                            seed_list=c(10, 33, 92), reg_weight=0., regularizer="cosine",
                            new_hier=FALSE, do_initial_fit=TRUE, verbose = T, regul_fixed = FALSE,
                            save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)

get_assigned_missing(test2, reference_cat=data1$beta[[1]], cutoff=0.6)
get_assigned_missing(test3, reference_cat=data1$beta[[1]])


ump = umap::umap(get_exposure(test3))
ump2 = umap::umap(data1$alpha[[1]] %>% dplyr::select())
plot.umap(ump, test3$groups, colors=gen_palette(n=length(unique(test3$groups))))

plot.umap = function(x, labels,
                    main="A UMAP visualization of the Iris dataset",
                    colors=c("#ff7f00", "#e377c2", "#17becf"),
                    pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
                    cex.main=1, cex.legend=0.85) {
  layout <- x
  if (is(x, "umap")) {
   layout <- x$layout
  }

  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
   par(mar=c(0.2,0.7,1.2,0.7), ps=10)
   plot(xylim, xylim, type="n", axes=F, frame=F)
   rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
        cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)

  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
   legend.pos <- "bottomleft"
   legend.text <- paste(as.character(labels.u), legend.suffix)
  }

  legend(legend.pos, legend=legend.text, inset=0.03,
        col=colors[as.integer(labels.u)],
      bty="n", pch=pch, cex=cex.legend)
}


cls2 = gen_palette(nrow(data1$beta[[1]]) + 3) %>% setNames(c(rownames(data1$beta[[1]]), "D4","D8","D9"))
convert_sigs_names(test2, reference_cat=data1$beta[[1]], cutoff=0.6) %>%
  plot_exposures(cls=cls2, add_centroid = T) %>% patchwork::wrap_plots(
    data1$alpha_plot[[1]] + scale_fill_manual(values=cls2), ncol=1
  )

cls3 = gen_palette(nrow(data1$beta[[1]])) %>% setNames(rownames(data1$beta[[1]]))
convert_sigs_names(test3, reference_cat=data1$beta[[1]]) %>%
  plot_exposures(cls=cls3, add_centroid = T) %>% patchwork::wrap_plots(
    data1$alpha_plot[[1]] + scale_fill_manual(values=cls3), ncol=1
  )

grps3 = test3$groups
grps_true = data1$counts[[1]]$groupid
aricode::NMI(grps3, grps_true)
aricode::ARI(grps3, grps_true)

compute.mse(m_inf=get_exposure(convert_sigs_names(test3, reference_cat=data1$beta[[1]])),
            m_true=data1$alpha[[1]] %>% dplyr::select(-sample, -groupid))
compute.cosine(m1=get_exposure(test3), m2=data1$alpha[[1]] %>% dplyr::select(-sample, -groupid),
               assigned_missing=get_assigned_missing(test3, reference_cat=data1$beta[[1]]),
               what="expos")


test2 %>% plot_similarity_reference(reference = data1$beta[[1]])
ggsave("./tmp.pdf", height = 20, width = 10)

generate_simulation_dataset(G=4, N=100, alpha_sigma=0.1,
                            sbs_catalogue=COSMIC_filt, private_sbs=private_sbs,
                            private_shared_sbs=private_shared_sbs, min_alpha=.2, py = py)

alpha %>%
  # tibble::rownames_to_column(var="sample") %>%
  reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=sample, y=value, fill=variable), stat="identity")

counts %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=sample, y=value, fill=variable), stat="identity")



alpha_prior = c(0.6,0.3,0.1)
beta = reticulate::r_to_py(COSMIC_filt[c("SBS1","SBS5","SBS17b"),])
n_muts = reticulate::r_to_py(sample(1000:5000,N))
data_alpha1 = py$generate_model(reticulate::r_to_py(alpha_prior), beta, n_muts, N=N, alpha_sigma=0.4,
                                seed=as.integer(15), use_normal=F)
data_alpha1$alpha %>% tibble::rownames_to_column() %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=rowname, y=value, fill=variable), stat="identity") +
  geom_hline(yintercept=1-cumsum(reticulate::py_to_r(alpha_prior)))


data_alpha2 = py$generate_model(reticulate::r_to_py(rev(alpha_prior)), beta, n_muts, N=N,
                                seed=as.integer(15), use_normal=T, alpha_sigma=0.4)

data_alpha1$alpha %>% dplyr::add_row(data_alpha2$alpha) %>%
  tibble::rownames_to_column() %>% reshape2::melt() %>%
  ggplot() + geom_bar(aes(x=as.integer(rowname), y=value, fill=variable), stat="identity") +
  geom_hline(yintercept=1-cumsum(alpha_prior))

counts = data_alpha1$data %>% dplyr::add_row(data_alpha2$data)

test = two_steps_inference(x=counts %>% dplyr::select(-sample), k=0:6, clusters=1:5,
                           enforce_sparsity2=T, reference_catalogue=COSMIC_filt[c("SBS1","SBS5"),],
                           lr=0.005, n_steps=2000, py=py,
                           hyperparameters=list("alpha_sigma"=0.05),
                           seed_list=c(10, 33, 92), reg_weight=0., regularizer="noreg",
                           new_hier=FALSE, do_initial_fit=FALSE,
                           save_runs_seed=TRUE, save_all_fits=TRUE, nonparametric=TRUE)


