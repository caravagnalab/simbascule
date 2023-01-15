#----------------------------------------------------------------------QC:PASSED
run.cohort <- function(
    cohort,
    k,
    lr,
    steps,
    phi,
    delta,
    input=TRUE
) {

  results <- NULL
  for (i in 1:nrow(cohort)) {

    cat('============================================\n') # TEST
    cat('                 Data No.', i, '\n') # TEST
    cat('============================================\n') # TEST

    xx <- simbasilica:::run.data(
      data = cohort[i, ],
      k = k,
      lr = lr,
      steps = steps,
      phi = phi,
      delta = delta,
      input = input
    )
    results <- rbind(results, xx)
  }
  return(results)
}
