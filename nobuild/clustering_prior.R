x = abs(rnorm(50, sd=.05))

tmp = lapply(x, function(i) {
  data.frame(
    y = rnorm(50, mean=.9, sd=i),
    x = i
  )
}) %>% Reduce(f=bind_rows)

tmp %>% ggplot() + geom_point(aes(x=x, y=y)) + geom_hline(yintercept = .9 + c(-0.05, 0.05))
