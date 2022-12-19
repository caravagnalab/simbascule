#----------------------------------------------------------------------QC:PASSED

evaluate.cohort <- function(x) {
  res <- NULL
  counter <- 1 # TEST
  for (i in 1:nrow(x)) {
    res <- rbind(res, evaluate.data(x = x[i, ]))
    #print(paste('counter:', counter)) # TEST
    counter <- counter + 1 # TEST
  }
  return(res)
}

