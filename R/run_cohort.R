#----------------------------------------------------------------------QC:PASSED
run.cohort <- function(
    x,
    k,
    cohort = "MyCohort",
    use_reference = TRUE,
    use_input = FALSE,
    lr = 0.01,
    steps = 500,
    max_iterations = 20,
    blacklist = NULL,
    phi = 0.05,
    delta = 0.9,
    filt_pi =0.1,
    groups = NULL,
    lambda_rate = NULL,
    sigma = FALSE,
    CUDA = FALSE,
    compile = TRUE,
    enforce_sparsity = FALSE
) {

  results <- NULL
  for (i in 1:nrow(x)) {

    cat('============================================\n') # TEST
    cat('                 Data No.', i, '\n') # TEST
    cat('============================================\n') # TEST

    xx <- simbasilica:::run.data(
      data=x[i, ],
      k=k,
      cohort = paste0(cohort, "-", i),
      use_reference = use_reference,
      use_input = use_input,
      lr = lr,
      steps = steps,
      max_iterations = max_iterations,
      blacklist = blacklist,
      phi = phi,
      delta = delta,
      filt_pi = filt_pi,
      groups = groups,
      lambda_rate = lambda_rate,
      sigma = sigma,
      CUDA = CUDA,
      compile = compile,
      enforce_sparsity = enforce_sparsity
    )
    results <- rbind(results, xx)
  }
  return(results)
}
