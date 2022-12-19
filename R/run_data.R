#----------------------------------------------------------------------QC:PASSED
run.data <- function(
    data,
    k,
    lr,
    steps,
    phi,
    delta,
    input=TRUE
) {
  
  x <- data$x[[1]]
  reference <- data$ref_cat[[1]]
  if (input) {
    input <- data$input_cat[[1]]
  } else {
    input = NULL
  }
  
  obj <- basilica::fit(
    x=x,
    reference_catalogue = reference,
    k = k,
    lr = lr,
    steps = steps,
    phi = phi,
    delta = delta,
    groups = NULL,
    input_catalogue = input
  )
  
  simulation.fit.obj <- tibble::add_column(
    data,
    #inf_exposure = list(obj$exposure),
    #inf_denovo = list(obj$denovo_signatures),
    #inf_fixed = list(obj$catalogue_signatures),
    #bic = obj$bic,
    #losses = list(obj$losses),
    fit = list(obj)
  )
  return(simulation.fit.obj)
}
