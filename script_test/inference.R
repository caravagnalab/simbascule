reticulate::use_condaenv("base")
pyb = reticulate::import("pybasilica")
pd = reticulate::import("pandas")

x = readRDS("./script_test/simulations/simul.N100.G5.s23.Rds")

counts = x$x[[1]]
groups_true = x$groups[[1]] -1
k_denovo_true = 9
k_fixed_true = 4
beta_fixed = pd$DataFrame(x$exp_fixed[[1]], index=rownames(x$exp_fixed[[1]]), columns=colnames(x$exp_fixed[[1]]))

# run.fit(data, k_list=k_denovo_true+k_fixed_true, groups=groups_true)
x.fit = pyb$fit(counts, k_list=as.integer(k_denovo_true+k_fixed_true), groups=groups_true)
x.fit = pyb$fit(counts, k_list=as.integer(k_denovo_true), groups=groups_true,
                 beta_fixed=beta_fixed)

N = x.fit$x %>% nrow()

beta = x.fit$beta_denovo %>%
  rbind(x.fit$beta_fixed)
beta_true = x$exp_fixed[[1]] %>% rbind(x$exp_denovo[[1]]) %>% as.data.frame()

alpha = x.fit$alpha$numpy() %>% as.data.frame(); colnames(alpha) = c(rownames(x.fit$beta_fixed), rownames(x.fit$beta_denovo))
alpha$group = groups_true +1
alpha_true = x$exp_exposure[[1]] %>% dplyr::mutate(group=x$groups[[1]])





plot_alpha(alpha)
plot_alpha(alpha_true)

plot_beta(beta)
plot_beta(beta_true)

simil = data.frame() %>%
  dplyr::mutate("sbs1"=as.character(NA), "sbs2"=as.character(NA), "simil"=as.numeric(NA))
for (b1 in rownames(beta))
  for (b2 in rownames(beta_true))
    simil = simil %>% dplyr::add_row(
      "sbs1"=b1, "sbs2"=b2, "simil"=lsa::cosine(unlist(beta[b1,]), unlist(beta_true[b2,]))[[1]]
    )

beta_pairs = simil %>%
  dplyr::group_by(sbs1) %>%
  dplyr::mutate(highest_cosine=max(simil)) %>%
  dplyr::filter(simil==highest_cosine) %>%
  dplyr::arrange(sbs2) %>%
  dplyr::mutate(pairs=paste0(sbs1,".",sbs2))

beta_pairs_unq = beta_pairs %>%
  dplyr::group_by(sbs2) %>%
  dplyr::mutate(nn=dplyr::n()) %>%
  dplyr::filter(nn==1) %>%
  dplyr::select(dplyr::contains("sbs"), "pairs") %>%
  dplyr::ungroup()


beta_pairs_dup = beta_pairs %>%
  dplyr::group_by(sbs2) %>%
  dplyr::mutate(nn=dplyr::n()) %>%
  dplyr::filter(nn>1) %>%
  dplyr::select(dplyr::contains("sbs"), "pairs") %>%
  dplyr::ungroup()


beta_all = rbind(beta %>% dplyr::mutate(type="estimated"), beta_true %>% dplyr::mutate(type="true"))

beta_all %>%
  tibble::rownames_to_column(var="sbs") %>%
  reshape2::melt(id=c("sbs", "type"), variable.name="context", value.name="beta") %>%
  dplyr::full_join(beta_pairs_unq %>% dplyr::select("sbs1","pairs"), by=c("sbs" = "sbs1")) %>%
  dplyr::full_join(beta_pairs_unq %>% dplyr::select("sbs2","pairs"), by=c("sbs" = "sbs2")) %>%
  dplyr::mutate(pairs=ifelse(is.na(pairs.x),pairs.y,pairs.x)) %>%
  dplyr::select(-pairs.x, -pairs.y) %>% tidyr::drop_na() %>%

  ggplot() +
  geom_bar(aes(x=context, y=beta), stat="identity") +
  facet_grid(pairs~type)


beta_all %>%
  tibble::rownames_to_column(var="sbs") %>%
  reshape2::melt(id=c("sbs", "type"), variable.name="context", value.name="beta") %>%
  dplyr::filter(sbs %in% c(beta_pairs_dup$sbs1, beta_pairs_dup$sbs2)) %>%

  ggplot() +
  geom_bar(aes(x=context, y=beta), stat="identity") +
  facet_wrap(sbs~type)





