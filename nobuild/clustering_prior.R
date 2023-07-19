x = abs(rnorm(50, sd=.05))

tmp = lapply(x, function(i) {
  data.frame(
    y = rnorm(50, mean=.9, sd=i),
    x = i
  )
}) %>% Reduce(f=bind_rows)

tmp %>% ggplot() + geom_point(aes(x=x, y=y)) + geom_hline(yintercept = .9 + c(-0.05, 0.05))



alpha_p = fit_cl.old$fit$params$alpha_prior / rowSums(fit_cl.old$fit$params$alpha_prior)
mm = alpha_p[1,1]

q.05 = mm - 0.15*mm
q.95 = mm + 0.15*mm
sd = (q.95-q.05) / (2*qnorm(0.95, mean=mm, sd=1))

sd_vals = c(sd, seq(from=0.001, to=0.2, length.out=10))
tmp = lapply(sd_vals, function(ss) {
  data.frame(
    dd = rnorm(50, mean=mm, sd=ss),
    sd_v = ss,
    mean_v = mm
  )
}) %>% do.call(rbind, .) %>% tibble::tibble()




tmp %>%
  dplyr::filter(sd_v==sd_vals[1]) %>%
  ggplot() +
  geom_histogram(aes(x=dd)) +
  geom_vline(aes(xintercept=mm-0.1*mm)) +
  geom_vline(aes(xintercept=mm+0.1*mm)) +
  facet_grid(~as.factor(sd_v))

