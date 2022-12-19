
#----------------------------------------------------------------------QC:PASSED
generate.theta <- function(mut_range, num_samples, seed=NULL) {
  # generate.theta(mut_range=1000:4000, num_samples=15, seed=NULL)
  if (!(is.integer(mut_range))) {
    stop("not valid mut_range (mutational range) value! e.g., mut_range=1000:4000")
  }
  
  if (!(is.numeric(num_samples))) {
    stop("not valid num_samples value! e.g., num_samples=45")
  }
  
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  
  theta = sample(mut_range, num_samples)  # integer
  return(theta)
}
