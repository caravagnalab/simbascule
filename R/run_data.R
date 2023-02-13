#----------------------------------------------------------------------QC:

run.data <- function(
    data,
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

  # catalogue data
  x <- data$x[[1]]

  # reference catalogue
  if (use_reference) {
    reference <- data$ref_cat[[1]]
  } else {
    reference <- NULL
  }

  # input catalogue
  if (use_input) {
    input_catalogue <- data$ref_cat[[1]]
  } else {
    input_catalogue = NULL #basilica::COSMIC_catalogue["SBS1", ]
  }

  # RUN START ------------------------------------------------------------------
  obj <- basilica::fit(
    x=x,
    k,
    reference_catalogue = reference,
    input_catalogue = input_catalogue,
    cohort,
    lr,
    steps,
    max_iterations,
    blacklist,
    phi,
    delta,
    filt_pi,
    groups,
    lambda_rate,
    sigma,
    CUDA,
    compile,
    enforce_sparsity
  )
  # RUN END --------------------------------------------------------------------

  fit.obj <- tibble::add_column(
    data,
    fit = list(obj)
  )
  return(fit.obj)
}
