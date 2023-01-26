#----------------------------------------------------------------------QC:PASSED

#' Title
#'
#' @param reference_path e.g., '/home/file.csv', path for reference catalogue file in csv format
#' @param ratio e.g., 0.7, ratio between reference and denovo catalogues, both derived from splitting main reference catalogue in two parts
#' @param targetX e.g., c(a,b), a: No. of fixed signatures and b: No. of denovo signatures in data
#' @param inputX e.g., c(a,b), a: No. of fixed signatures and b: No. of denovo signatures in input
#' @param similarity_limit similarity between signatures in catalogues (trying to include non-similar signatures)
#' @param groups e.g., rep(1, n), n is number of patients in data (1st one is always one)
#' @param mut_range # range of mutations counts in patients
#' @param seed obvious
#' @param num_data e.g., 50, No. of the data (rows) in cohort. each data is a mutational catalogues comprising several patients
#'
#' @return a cohort includes several data
#' @export
#'
#' @examples
generate.cohort <- function(
    reference,
    ratio,
    targetX,
    inputX,
    similarity_limit,
    groups,
    mut_range,
    seed = NULL,
    num_data
) {

  reference_denovo <- simbasilica:::split.reference(reference=reference, ratio=ratio, seed=seed)

  reference_catalogue <- reference_denovo$reference
  denovo_catalogue <- reference_denovo$denovo

  reference_cosine <- simbasilica:::cosine.matrix(reference_catalogue, reference_catalogue)
  denovo_cosine <- simbasilica:::cosine.matrix(denovo_catalogue, denovo_catalogue)

  data <- NULL

  if (is.null(seed)) {
    for (i in 1:num_data) {
      xx <- simbasilica:::generate.data(
        reference_catalogue = reference_catalogue,
        denovo_catalogue = denovo_catalogue,
        reference_cosine = reference_cosine,
        denovo_cosine = denovo_cosine,
        targetX = targetX,
        inputX = inputX,
        similarity_limit = similarity_limit,
        groups = groups,
        mut_range = mut_range,
        seed = seed
      )
      data <- rbind(data, xx)
    }
  } else {
    for (i in 1:num_data) {
      xx <- simbasilica:::generate.data(
        reference_catalogue = reference_catalogue,
        denovo_catalogue = denovo_catalogue,
        reference_cosine = reference_cosine,
        denovo_cosine = denovo_cosine,
        targetX = targetX,
        inputX = inputX,
        similarity_limit = similarity_limit,
        groups = groups,
        mut_range = mut_range,
        seed = seed
      )
      data <- rbind(data, xx)
      seed <- seed + 1
    }
  }

  return(data)
}
