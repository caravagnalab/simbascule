#----------------------------------------------------------------------QC:PASSED
#' @import dplyr
evaluate.data <- function(x) {
  #--------------------------
  reference <- x$ref_cat[[1]]
  input <- x$input_cat[[1]]
  expected_fixed <- x$exp_fixed[[1]]
  #inferred_fixed <- x$inf_fixed[[1]]
  inferred_fixed <- x$fit[[1]]$catalogue_signatures
  a <- simbasilica:::fixed.accuracy(reference, input, expected_fixed, inferred_fixed)
  TP <- a$TP
  FP <- a$FP
  TN <- a$TN
  FN <- a$FN
  accuracy <- (TP + TN) / (TP + TN + FP + FN)
  #--------------------------
  m <- x$x[[1]]
  #alpha <- x$inf_exposure[[1]]
  alpha <- x$fit[[1]]$exposure
  beta <- rbind(x$fit[[1]]$catalogue_signatures, x$fit[[1]]$denovo_signatures)
  mr <- simbasilica:::reconstruct.count(m, alpha, beta)
  mae <- simbasilica:::compute.mae(m, mr)
  mse <- simbasilica:::compute.mse(m, mr)
  #--------------------------
  b <- simbasilica:::denovo.similarity(x$exp_denovo[[1]], x$fit[[1]]$denovo_signatures)
  denovo_similarity <- b$similarity_average  # numeric
  denovo_match <- b$match_df                 # data.frame
  #--------------------------
  denovo_ratio <- simbasilica:::denovo.ratio(x$exp_denovo[[1]], x$fit[[1]]$denovo_signatures)
  #--------------------------

  # CREATE TIBBLE ----------------------------
  obj <- tibble::tibble(

    targetX = x$targetX,
    inputX = x$inputX,
    num_samples = nrow(m),

    mae = mae,
    mse = mse,
    fixed_acc = accuracy,
    denovo_ratio = denovo_ratio,
    denovo_sim = denovo_similarity,
    denovo_match = list(denovo_match),
  )

  return(obj)
}
