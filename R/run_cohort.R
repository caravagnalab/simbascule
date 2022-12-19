#----------------------------------------------------------------------QC:PASSED
run.cohort <- function(
    cohort,
    k,
    lr,
    steps,
    phi,
    delta,
    lambda_rate = NULL,
    sigma = FALSE,
    input=TRUE
) {
  
  results <- NULL
  for (i in 1:nrow(cohort)) {
    
    cat('============================================\n') # TEST
    cat('                 Data No.', i, '\n') # TEST
    cat('============================================\n') # TEST
    
    xx <- basilica:::run.data(
      data = cohort[i, ],
      k = k,
      lr = lr,
      steps = steps,
      phi = phi,
      delta = delta,
      lambda_rate = lambda_rate,
      sigma = sigma,
      input = input
    )
    results <- rbind(results, xx)
  }
  return(results)
}
