#' ## Wrapping fns ####
#'
#' #----------------------------------------------------------------------QC:PASSED
#'
#' #' Title
#' #'
#' #' @param reference_path e.g., '/home/file.csv', path for reference catalogue file in csv format
#' #' @param ratio e.g., 0.7, ratio between reference and denovo catalogues, both derived from splitting main reference catalogue in two parts
#' #' @param targetX e.g., c(a,b), a: No. of fixed signatures and b: No. of denovo signatures in data
#' #' @param inputX e.g., c(a,b), a: No. of fixed signatures and b: No. of denovo signatures in input
#' #' @param similarity_limit similarity between signatures in catalogues (trying to include non-similar signatures)
#' #' @param groups e.g., rep(1, n), n is number of patients in data (1st one is always one)
#' #' @param mut_range # range of mutations counts in patients
#' #' @param seed obvious
#' #' @param num_data e.g., 50, No. of the data (rows) in cohort. each data is a mutational catalogues comprising several patients
#' #'
#' #' @return a cohort includes several data
#' #' @export
#' #'
#' generate.cohort <- function(
#'     reference,
#'     ratio,
#'     targetX,
#'     inputX,
#'     similarity_limit,
#'     groups,
#'     mut_range,
#'     seed = NULL,
#'     num_data
#' ) {
#'
#'   reference_denovo <- simbasilica:::split.reference(reference=reference, ratio=ratio, seed=seed)
#'
#'   reference_catalogue <- reference_denovo$reference
#'   denovo_catalogue <- reference_denovo$denovo
#'
#'   reference_cosine <- simbasilica:::cosine.matrix(reference_catalogue, reference_catalogue)
#'   denovo_cosine <- simbasilica:::cosine.matrix(denovo_catalogue, denovo_catalogue)
#'
#'   data <- NULL
#'
#'   if (is.null(seed)) {
#'     for (i in 1:num_data) {
#'       xx <- simbasilica:::generate.data(
#'         reference_catalogue = reference_catalogue,
#'         denovo_catalogue = denovo_catalogue,
#'         reference_cosine = reference_cosine,
#'         denovo_cosine = denovo_cosine,
#'         targetX = targetX,
#'         inputX = inputX,
#'         similarity_limit = similarity_limit,
#'         groups = groups,
#'         mut_range = mut_range,
#'         seed = seed
#'       )
#'       data <- rbind(data, xx)
#'     }
#'   } else {
#'     for (i in 1:num_data) {
#'       xx <- simbasilica:::generate.data(
#'         reference_catalogue = reference_catalogue,
#'         denovo_catalogue = denovo_catalogue,
#'         reference_cosine = reference_cosine,
#'         denovo_cosine = denovo_cosine,
#'         targetX = targetX,
#'         inputX = inputX,
#'         similarity_limit = similarity_limit,
#'         groups = groups,
#'         mut_range = mut_range,
#'         seed = seed
#'       )
#'       data <- rbind(data, xx)
#'       seed <- seed + 1
#'     }
#'   }
#'
#'   return(data)
#' }
#'
#'
#' #----------------------------------------------------------------------QC:PASSED
#' generate.data <- function(
#'     reference_catalogue,
#'     denovo_catalogue,
#'     targetX,
#'     inputX,
#'     similarity_limit,
#'     groups,
#'     mut_range,
#'     reference_cosine = NULL,
#'     denovo_cosine = NULL,
#'     private_sigs=list("rare"=c(), "common"=c()),
#'     private_fracs=list("rare"=1., "common"=0.),
#'     thr = 0.1,
#'     seed=NULL
#' ) {
#'
#'   if (!is.null(seed)) set.seed(seed = seed)
#'
#'   # SIGNATURES --------------------------------
#'   signatures <- generate.signatures(
#'     reference_catalogue = reference_catalogue,
#'     denovo_catalogue  = denovo_catalogue,
#'     reference_cosine = reference_cosine,
#'     denovo_cosine = denovo_cosine,
#'     complexity = targetX,
#'     similarity_limit = similarity_limit,
#'     seed = seed
#'   )
#'   beta <- rbind(signatures$fixed, signatures$denovo)
#'
#'   # INPUT ------------------------------------
#'   input <- generate.input(
#'     reference_catalogue = reference_catalogue,
#'     beta_fixed = signatures$fixed,
#'     complexity = inputX,
#'     seed = seed
#'   )
#'
#'   # EXPOSURE ---------------------------------
#'   alpha <- generate.exposure(beta=beta,
#'                              groups=groups,
#'                              private_sigs=private_sigs,
#'                              private_fracs=private_fracs,
#'                              seed=seed,
#'                              thr=thr) # include group column
#'
#'   # removing group column
#'   if (!is.null(alpha$group)) {
#'     groups = alpha$group
#'     alpha <- subset(alpha, select = -c(group))
#'   }
#'
#'   # apply exposure limit (<0.05) to one fixed and one de-novo signature
#'   # alpha <- simbasilica:::edit.exposure(alpha = alpha)
#'
#'   cat("ALPHA DONE\n")
#'
#'   # THETA ------------------------------------
#'   num_samples <- length(groups)
#'   theta <- generate.theta(mut_range=mut_range, num_samples=num_samples, seed=seed)
#'
#'   cat("THETA DONE\n")
#'
#'   # COUNT MATRIX -----------------------------
#'   # m <- generate.counts(alpha=alpha, beta=beta, theta=theta, seed=seed)
#'
#'   # M <- as.data.frame(round(as.matrix(alpha*theta) %*% as.matrix(beta), digits = 0))
#'   rate = round(as.matrix(alpha*theta) %*% as.matrix(beta), digits = 0)
#'   M = sapply(colnames(beta), function(s)
#'     sapply(1:num_samples, function(n)
#'       rpois(1, rate[n, s]))) %>% as.data.frame()
#'
#'   # print(rate)
#'   # print(rowSums(M) < mut_range[1])
#'   #
#'   # while (any(rowSums(M) < mut_range[1])) {
#'   #   M = sapply(colnames(beta), function(s)
#'   #     sapply(1:num_samples, function(n)
#'   #       rpois(1, rate[n, s]))) %>% as.data.frame()
#'   # }
#'
#'   rownames(M) <- rownames(alpha)
#'   colnames(M) <- colnames(beta)
#'
#'   cat("COUNTS DONE\n")
#'
#'   # MODIFY COMPLEXITY VALUES -----------------
#'   if (is.numeric(targetX)) {
#'     targetX <- paste(targetX, collapse = "|")
#'   }
#'   if (is.numeric(inputX)) {
#'     inputX <- paste(inputX, collapse = "|")
#'   }
#'
#'   # CREATE TIBBLE ----------------------------
#'   obj <- tibble::tibble(
#'     x = list(M),
#'     input_cat = list(input),
#'     ref_cat = list(reference_catalogue),
#'     exp_exposure = list(alpha),
#'     exp_fixed = list(signatures$fixed),
#'     exp_denovo = list(signatures$denovo),
#'     targetX = targetX,
#'     inputX = inputX,
#'     groups = list(groups)
#'   )
#'   return(obj)
#' }
#'
#'
#' #----------------------------------------------------------------------QC:PASSED
#'
#' generate.input <- function(
#'     reference_catalogue,
#'     beta_fixed,
#'     complexity,
#'     seed=NULL
#' ) {
#'
#'   if (!is.null(seed)) {
#'     set.seed(seed = seed)
#'   }
#'
#'   k_fixed <- nrow(beta_fixed)
#'
#'   if (is.null(complexity)) {
#'     return(NULL)
#'   } else if (is.numeric(complexity) & length(complexity)==2) {
#'     n_overlap <- complexity[1]
#'     n_extra <- complexity[2]
#'     if (n_overlap > k_fixed) {
#'       stop("overlap signatures are more than fixed signatures!")
#'     }
#'   } else if (is.character(complexity)) {
#'     if (complexity=='low') {
#'       n_overlap <- sample(1:k_fixed, 1)
#'       n_extra <- 0
#'     }
#'     else if (complexity=='medium') {
#'       n_overlap <- sample(1:k_fixed, 1)
#'       n_extra <- sample(1:k_fixed, 1)
#'     }
#'     else if (complexity=='high') {
#'       n_overlap <- 0
#'       n_extra <- sample(1:(k_fixed), 1)
#'     }
#'     else {
#'       stop("complexity argument should be selected from {'low', 'medium', 'high'}")
#'     }
#'   } else {
#'     stop("wrong complexity argument!")
#'   }
#'
#'   extra_ref <- reference_catalogue[setdiff(rownames(reference_catalogue), rownames(beta_fixed)), ]
#'
#'   if (n_overlap > 0) {
#'     overlap <- sample(rownames(beta_fixed))[1:n_overlap]
#'   } else {
#'     overlap <- NULL
#'   }
#'
#'   if (n_extra > 0) {
#'     extra <- sample(rownames(extra_ref))[1:n_extra]
#'   } else {
#'     extra <- NULL
#'   }
#'
#'   df <- reference_catalogue[c(overlap, extra), ]
#'
#'   return(df)
#' }
#'
#'
#'
#' ## Run ####
#'
#' #----------------------------------------------------------------------QC:PASSED
#' run.cohort <- function(
#'     x,
#'     k,
#'     cohort = "MyCohort",
#'     use_reference = TRUE,
#'     use_input = FALSE,
#'     lr = 0.01,
#'     steps = 500,
#'     max_iterations = 20,
#'     blacklist = NULL,
#'     phi = 0.05,
#'     delta = 0.9,
#'     filt_pi =0.1,
#'     groups = NULL,
#'     lambda_rate = NULL,
#'     sigma = FALSE,
#'     CUDA = FALSE,
#'     compile = TRUE,
#'     enforce_sparsity = FALSE
#' ) {
#'
#'   results <- NULL
#'   for (i in 1:nrow(x)) {
#'
#'     cat('============================================\n') # TEST
#'     cat('                 Data No.', i, '\n') # TEST
#'     cat('============================================\n') # TEST
#'
#'     xx <- simbasilica:::run.data(
#'       data=x[i, ],
#'       k=k,
#'       cohort = paste0(cohort, "-", i),
#'       use_reference = use_reference,
#'       use_input = use_input,
#'       lr = lr,
#'       steps = steps,
#'       max_iterations = max_iterations,
#'       blacklist = blacklist,
#'       phi = phi,
#'       delta = delta,
#'       filt_pi = filt_pi,
#'       groups = groups,
#'       lambda_rate = lambda_rate,
#'       sigma = sigma,
#'       CUDA = CUDA,
#'       compile = compile,
#'       enforce_sparsity = enforce_sparsity
#'     )
#'     results <- rbind(results, xx)
#'   }
#'   return(results)
#' }
#'
#'
#' #----------------------------------------------------------------------QC:
#'
#' run.data <- function(
#'     data,
#'     k,
#'     cohort = "MyCohort",
#'     use_reference = TRUE,
#'     use_input = FALSE,
#'     lr = 0.01,
#'     steps = 500,
#'     max_iterations = 20,
#'     blacklist = NULL,
#'     phi = 0.05,
#'     delta = 0.9,
#'     filt_pi =0.1,
#'     groups = NULL,
#'     lambda_rate = NULL,
#'     sigma = FALSE,
#'     CUDA = FALSE,
#'     compile = TRUE,
#'     enforce_sparsity = FALSE
#' ) {
#'
#'   # catalogue data
#'   x <- data$x[[1]]
#'
#'   # reference catalogue
#'   if (use_reference) {
#'     reference <- data$ref_cat[[1]]
#'   } else {
#'     reference <- NULL
#'   }
#'
#'   # input catalogue
#'   if (use_input) {
#'     input_catalogue <- data$ref_cat[[1]]
#'   } else {
#'     input_catalogue = NULL #basilica::COSMIC_catalogue["SBS1", ]
#'   }
#'
#'   # RUN START ------------------------------------------------------------------
#'   obj <- basilica::fit(
#'     x=x,
#'     k,
#'     reference_catalogue = reference,
#'     input_catalogue = input_catalogue,
#'     cohort,
#'     lr,
#'     steps,
#'     max_iterations,
#'     blacklist,
#'     phi,
#'     delta,
#'     filt_pi,
#'     groups,
#'     lambda_rate,
#'     sigma,
#'     CUDA,
#'     compile,
#'     enforce_sparsity
#'   )
#'   # RUN END --------------------------------------------------------------------
#'
#'   fit.obj <- tibble::add_column(
#'     data,
#'     fit = list(obj)
#'   )
#'   return(fit.obj)
#' }
#'
#'
#'
#'
#' ## Counts ####
#'
#' # may use theta*alpha*beta or generate.counts (should be discussed to )
#' generate.counts <- function(alpha, beta, theta, seed=NULL) {
#'
#'   if (!is.null(seed)) {
#'     set.seed(seed = seed)
#'   }
#'
#'   if (!is.null(alpha$group)) {
#'     alpha <- subset(alpha, select = -c(group))
#'   }
#'
#'   alpha <- alpha[, order(colnames(alpha))]
#'   beta <- beta[order(rownames(beta)), ]
#'   if (!identical(colnames(alpha), rownames(beta))) {
#'     print(colnames(alpha))
#'     print(rownames(beta))
#'     stop("alpha and beta are NOT valid!")
#'   }
#'
#'   num_samples <- nrow(alpha)
#'
#'   M <- matrix(rep(0, num_samples*96) , ncol = 96)
#'
#'   # iterate over samples
#'   for (sample in 1:num_samples) {
#'
#'     p <- as.numeric(alpha[sample, ]) # select sample i
#'
#'     # iterate over number of mutations in sample i
#'     for (j in 1:theta[sample]) {
#'
#'       # sample signature profile index from categorical data
#'       signature_idx <- extraDistr::rcat(1, p)
#'       signature <- beta[signature_idx, ]
#'
#'       # sample mutation feature index for corresponding signature from categorical data
#'       mutation_idx <- extraDistr::rcat(1, as.numeric(signature))
#'
#'       # add +1 to the mutation feature in position j in branch i
#'       M[sample, mutation_idx] <- M[sample, mutation_idx] + 1
#'
#'     }
#'   }
#'   M <- as.data.frame(M)
#'   colnames(M) <- colnames(beta)
#'   rownames(M) <- rownames(alpha)
#'   return(M)
#' }
#'
#'
#' ## Signatures ####
#' #----------------------------------------------------------------------QC:PASSED
#' # generate signatures which includes:
#' #   * fixed signatures (SBS1 included)
#' #   * denovo signatures
#' # similarity between catalogue signatures are less than threshold
#' # similarity between denovo signatures are less than threshold
#' # but similarity between catalogue and denovo signatures are not taken to the account (future work)
#' generate.signatures <- function(
#'     reference_catalogue,
#'     denovo_catalogue,
#'     reference_cosine = NULL, # cosine similarity matrix of reference signatures (SBS1 excluded)
#'     denovo_cosine = NULL,    # cosine similarity matrix of denovo signatures
#'     complexity = NULL,       # c(3,2) | "low" | "medium" | "high" |
#'     similarity_limit = NULL,
#'     seed = NULL
#' ) {
#'
#'   if (complexity == -1)
#'     return(list(fixed = reference_catalogue, denovo = denovo_catalogue))
#'
#'   if (!is.null(seed)) {
#'     set.seed(seed = seed)
#'   }
#'
#'   if (is.numeric(complexity) & length(complexity)==2) {
#'     fixed_num <- complexity[1]
#'     denovo_num <- complexity[2]
#'   } else if (is.character(complexity)) {
#'     if (complexity=='low') {
#'       fixed_num <- sample(3:5, 1)
#'       denovo_num <- sample(0:2, 1)
#'     }
#'     else if (complexity=='medium') {
#'       fixed_num <- sample(1:2, 1)
#'       denovo_num <- sample(3:5, 1)
#'     }
#'     else if (complexity=='high') {
#'       fixed_num <- sample(3:5, 1)
#'       denovo_num <- sample(3:5, 1)
#'     }
#'     else {
#'       stop("complexity argument should be selected from {'low', 'medium', 'high'}")
#'     }
#'   } else {
#'     stop("wrong complexity argument!")
#'   }
#'
#'   SBS1 <- reference_catalogue['SBS1', ] # save SBS1 (data.frame)
#'   reference <- reference_catalogue[!(rownames(reference_catalogue) %in% c("SBS1")), ] # excludes SBS1
#'
#'   # catalogue signatures -------------------------------------------------------
#'
#'   if (fixed_num > 1) {
#'
#'     while (TRUE) {
#'       shuffled_reference = reference[sample(1:nrow(reference)), ]
#'       signatures <- rownames(shuffled_reference[1:(fixed_num-1), ])
#'       cos_matrix <- reference_cosine[c("SBS1", signatures), c("SBS1", signatures)]
#'       for (i in 1:nrow(cos_matrix)) {
#'         cos_matrix[i, i] <- 0
#'       }
#'       max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
#'       if (cos_matrix[max][1] < similarity_limit) {
#'         fixed_df <- rbind(SBS1, reference[signatures, ])
#'         break
#'       }
#'     }
#'   }
#'   else {
#'     fixed_df <- SBS1
#'   }
#'
#'   # denovo signatures ----------------------------------------------------------
#'
#'   if (denovo_num > 1) {
#'
#'     while (TRUE) {
#'       shuffled_denovo = denovo_catalogue[sample(1:nrow(denovo_catalogue)), ]
#'       signatures <- rownames(shuffled_denovo[1:denovo_num, ])
#'       cos_matrix <- denovo_cosine[signatures, signatures]
#'       for (i in 1:nrow(cos_matrix)) {
#'         cos_matrix[i, i] <- 0
#'       }
#'       max = which(cos_matrix == max(cos_matrix), arr.ind = TRUE)
#'       if (cos_matrix[max][1] < similarity_limit) {
#'         denovo_df <- denovo_catalogue[signatures, ]
#'         rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
#'         break
#'       }
#'     }
#'   }
#'   else if (denovo_num==1) {
#'     shuffled_denovo = denovo_catalogue[sample(1:nrow(denovo_catalogue)), ]
#'     denovo_df <- shuffled_denovo[1, ]
#'     rownames(denovo_df) <- paste0(rownames(denovo_df), "_D")
#'   }
#'   else {
#'     denovo_df <- NULL
#'   }
#'
#'   obj <- list(fixed = fixed_df, denovo = denovo_df)
#'   return(obj)
#' }
#'
#'
#'
#'
#' ## Exposures ####
#' library(gtools)
#' library(fdrtool)
#'
#' generate.exposure <- function(beta, groups, private_sigs, private_fracs,
#'                               n_priv_per_group=1, seed=NULL, thr=0.1) {
#'
#'   if (is.null(n_priv_per_group)) n_priv_per_group = round(length(private_sigs.all)/2)
#'
#'   signatures <- rownames(beta)
#'   if (!('SBS1' %in% signatures)) stop('Wrong signatures! SBS1 not included!')
#'
#'   if (!is.null(seed)) set.seed(seed = seed)
#'
#'   if (length(signatures) < 2) stop("not valid! there are not enough signatures!")
#'
#'   df_list <- list()
#'
#'   private_sigs.all = private_sigs.start = c(private_sigs$common, private_sigs$rare)
#'   shared_sigs = setdiff(signatures, private_sigs.all)
#'
#'   unq = unique(groups) %>% sort()
#'
#'   data = data.frame(matrix(ncol=length(signatures)+1, nrow=0))
#'   colnames(data) = c(signatures, "group")
#'
#'   for (i in 1:length(unq)) {
#'     group = unq[i]
#'     # print(group)
#'     # print(private_sigs.all)
#'
#'     if (length(unique(groups))==1) {
#'       sigNums <- length(signatures)
#'       sigNames <- signatures
#'     } else {
#'
#'       sigNames_priv = c()
#'
#'       # number of sigs to assign to "group"
#'       if (length(private_sigs.all)>0) {
#'         if (i == length(unq))
#'           sigNums_priv = length(private_sigs.all) else
#'             sigNums_priv = base::sample(1:n_priv_per_group, 1)
#'
#'           sigNames_priv = base::sample(private_sigs.all, sigNums_priv)
#'
#'           sigNames_priv = c(sigNames_priv, private_sigs.start[sigNames_priv],
#'                             private_sigs.start[private_sigs.start==sigNames_priv] %>% names) %>%
#'             purrr::discard(function(x) x=="" || is.na(x)) %>% unique()
#'
#'           private_sigs.all = private_sigs.all %>% purrr::discard(function(x) x %in% sigNames_priv)
#'
#'           sigNames = c(shared_sigs, sigNames_priv)
#'           sigNums = length(sigNames)
#'       }
#'     }
#'
#'     num_samples <- length(groups[groups==group])
#'
#'     mean_prior = fdrtool::rhalfnorm(sigNums, 1) %>% setNames(sigNames)
#'     alpha = sapply(1:sigNums, function(s)
#'       sapply(1:num_samples, function(x) sample_positive_norm(mean=mean_prior[s], sd=1))
#'     )
#'     alpha = apply(alpha, 1, function(x) x/sum(x)) %>% t() %>% as.data.frame()
#'     colnames(alpha) <- sigNames
#'
#'     ## change values for rare/common
#'
#'     alpha = adjust_frequency(alpha, N=length(groups),
#'                              columns=intersect(private_sigs$rare, colnames(alpha)),
#'                              frac=private_fracs$rare, check=">=",
#'                              # mean=rep(0, length.out=length(intersect(private_sigs$rare, colnames(alpha)))) %>%
#'                              #   setNames(intersect(private_sigs$rare, colnames(alpha))),
#'                              mean=rep(0, length.out=ncol(alpha)) %>% setNames(colnames(alpha)),
#'                              sd=0, min_exp=0.3)
#'
#'     alpha = adjust_frequency(alpha, N=length(groups),
#'                              columns=intersect(private_sigs$common, colnames(alpha)),
#'                              frac=private_fracs$common, check="<",
#'                              mean=mean_prior, sd=1, thr=thr, min_exp=0.1)
#'
#'     alpha = adjust_frequency(alpha, N=length(groups),
#'                              columns=shared_sigs,
#'                              frac=nrow(alpha), check="<",
#'                              mean=mean_prior,
#'                              sd=1, thr=thr, min_exp=0.1)
#'
#'     alpha = apply(alpha, 1, function(x) x/sum(x)) %>% t() %>% as.data.frame()
#'
#'     # alpha = alpha / rowSums(alpha). ## check if it works
#'     alpha$group <- group
#'
#'     data = data %>%
#'       dplyr::add_row(alpha)
#'
#'   }
#'
#'   # merge all different group exposure matrices
#'   # data <- Reduce(function(x, y) { print(c(x,y)); merge(x, y, all=TRUE) }, df_list)
#'
#'   # sort columns
#'   column_names <- colnames(data)
#'   #column_names <- column_names[order(column_names)]
#'   column_names <- append(setdiff(column_names, "group"), "group")
#'   #column_names[length(column_names)+1] <- "group"
#'   data <- data[, column_names]
#'
#'   data[is.na(data)] <- 0    # convert 'NA' to zero
#'   data[order(data$group), ] # sort rows by group column
#'
#'   return(data)
#' }
#'
#'
#' adjust_frequency = function(alpha, N, columns, frac, check, mean, sd, thr=0., min_exp=0.2) {
#'   if (!is.integer(frac)) { nn = min(nrow(alpha), round(N * frac)) } else { nn = frac}
#'   for (colname in columns) {
#'     alpha = alpha %>% tibble::as_tibble() %>% dplyr::mutate(!!colname:=replace( .[[colname]], {{colname}}<thr, 0 ))
#'
#'     if ( !match.fun(check)(sum(alpha[,colname]>0), nn) ) next
#'     samples.tmp = sample(1:nrow(alpha), size=nn)
#'     alpha[-samples.tmp, colname] = sample_positive_norm(mean[colname], sd)
#'     while (any(alpha[samples.tmp, colname] < min_exp))
#'       alpha[samples.tmp, colname] = sample_positive_norm(mean[colname], sd=1)
#'   }
#'
#'   return(alpha)
#' }
#'
#'
#' sample_positive_norm = function(mean, sd) {
#'   while(TRUE) {
#'     n = rnorm(1, mean, sd)
#'     if (n >= 0)
#'       return(n)
#'   }
#' }
#'
#'
#' ## Theta ####
#'
#' #----------------------------------------------------------------------QC:PASSED
#' generate.theta <- function(mut_range, num_samples, seed=NULL) {
#'   # generate.theta(mut_range=1000:4000, num_samples=15, seed=NULL)
#'   if (!(is.integer(mut_range))) {
#'     stop("not valid mut_range (mutational range) value! e.g., mut_range=1000:4000")
#'   }
#'
#'   if (!(is.numeric(num_samples))) {
#'     stop("not valid num_samples value! e.g., num_samples=45")
#'   }
#'
#'   if (!is.null(seed)) {
#'     set.seed(seed = seed)
#'   }
#'
#'   theta = sample(mut_range, num_samples)  # integer
#'   return(theta)
#' }
#'
#'
#'
#' ## Other ####
#'
#' #----------------------------------------------------------------------QC:PASSED
#' # split reference catalogue to 2 sub catalogue:
#' # reference catalogue (SBS1 included) + denovo catalogue
#' split.reference <- function(reference, ratio, seed=NULL) {
#'
#'   if (!is.null(seed)) {
#'     set.seed(seed = seed)
#'   }
#'
#'   # read csv file as data.frame
#'   #reference <- read.table(reference_path, sep = ",", row.names = 1, header = TRUE, check.names = FALSE)
#'   num_ref <- round(ratio * nrow(reference))
#'
#'   SBS1 <- reference['SBS1', ]   # save SBS1 (data.frame)
#'   reference <- reference[!(rownames(reference) %in% c("SBS1")), ] # excludes SBS1
#'
#'   # shuffle the reference catalogue
#'   shuffled_reference = reference[sample(1:nrow(reference)), ]
#'
#'   ref <- shuffled_reference[1:(num_ref-1), ]
#'   ref <- ref[order(rownames(ref)), ]
#'   ref <- rbind(SBS1, ref) # includes SBS1
#'
#'   denovo <- shuffled_reference[num_ref:nrow(shuffled_reference), ]
#'   denovo <- denovo[order(rownames(denovo)), ]
#'
#'   obj <- list(reference=ref, denovo=denovo)
#'   return(obj)
#' }
#'
#'
#'
#'
#' ## Evaluation ####
#'
#' evaluate.cohort <- function(x) {
#'   res <- NULL
#'   counter <- 1 # TEST
#'   for (i in 1:nrow(x)) {
#'     res <- rbind(res, evaluate.data(x = x[i, ]))
#'     #print(paste('counter:', counter)) # TEST
#'     counter <- counter + 1 # TEST
#'   }
#'   return(res)
#' }
#'
#' #' @import dplyr
#' evaluate.data <- function(x) {
#'   #--------------------------
#'   reference <- x$ref_cat[[1]]
#'   input <- x$input_cat[[1]]
#'   expected_fixed <- x$exp_fixed[[1]]
#'   #inferred_fixed <- x$inf_fixed[[1]]
#'   inferred_fixed <- x$fit[[1]]$catalogue_signatures
#'   a <- simbasilica:::fixed.accuracy(reference, input, expected_fixed, inferred_fixed)
#'   TP <- a$TP
#'   FP <- a$FP
#'   TN <- a$TN
#'   FN <- a$FN
#'   accuracy <- (TP + TN) / (TP + TN + FP + FN)
#'   #--------------------------
#'   m <- x$x[[1]]
#'   #alpha <- x$inf_exposure[[1]]
#'   alpha <- x$fit[[1]]$exposure
#'   beta <- rbind(x$fit[[1]]$catalogue_signatures, x$fit[[1]]$denovo_signatures)
#'   mr <- simbasilica:::reconstruct.count(m, alpha, beta)
#'   mae <- simbasilica:::compute.mae(m, mr)
#'   mse <- simbasilica:::compute.mse(m, mr)
#'   #--------------------------
#'   b <- simbasilica:::denovo.similarity(x$exp_denovo[[1]], x$fit[[1]]$denovo_signatures)
#'   denovo_similarity <- b$similarity_average  # numeric
#'   denovo_match <- b$match_df                 # data.frame
#'   #--------------------------
#'   denovo_ratio <- simbasilica:::denovo.ratio(x$exp_denovo[[1]], x$fit[[1]]$denovo_signatures)
#'   #--------------------------
#'
#'   # CREATE TIBBLE ----------------------------
#'   obj <- tibble::tibble(
#'
#'     targetX = x$targetX,
#'     inputX = x$inputX,
#'     num_samples = nrow(m),
#'
#'     mae = mae,
#'     mse = mse,
#'     fixed_acc = accuracy,
#'     denovo_ratio = denovo_ratio,
#'     denovo_sim = denovo_similarity,
#'     denovo_match = list(denovo_match),
#'   )
#'
#'   return(obj)
#' }
#'
#'
#' ## Utils ####
#' cosine.vector = function(a, b) {
#'
#'   if (!identical(colnames(a), colnames(b))) {
#'     a = a[names(b)]
#'   }
#'
#'   numerator = sum(a * b)
#'   denominator = sqrt(sum(a^2)) * sqrt(sum(b^2))
#'   return(numerator / denominator)
#' }
#'
#'
#' cosine.matrix = function(a, b) {
#'   # a and b are data.frame
#'
#'   df = data.frame(matrix(0, nrow(a), nrow(b)))
#'   rownames(df) = rownames(a)
#'   colnames(df) = rownames(b)
#'
#'   for (i in 1:nrow(a)) {
#'     denovo = a[i, ]
#'     for (j in 1:nrow(b)) {
#'       ref = b[j, ]
#'
#'       score = cosine.vector(denovo, ref)
#'       df[i,j] = score
#'     }
#'   }
#'
#'   return(df)
#' }
#'
#'
#' edit.exposure = function(alpha) {
#'
#'   n = nrow(alpha)
#'   val = 0.05
#'
#'   while ( TRUE ) {
#'
#'     exp = alpha
#'
#'     # fixed
#'     newFixed = runif(n, 0, val)
#'     fixed_diff = exp[, 4] - newFixed
#'     exp[, 4] = newFixed
#'     exp[, 3] = exp[, 3] + fixed_diff
#'
#'     # de-novo
#'     newDenovo = runif(n, 0, val)
#'     denovo_diff = exp[, 7] - newDenovo
#'     exp[, 7] = newDenovo
#'     exp[, 6] = exp[, 6] + denovo_diff
#'
#'     if ( sum(rowSums(exp))==n & all(apply(exp, 1, function(x) all(x > 0))) ) {
#'       return(exp)
#'     }
#'   }
#' }
#'
#'
#' fixed.accuracy = function(reference, input, expected_fixed, inferred_fixed) {
#'   ref_list = rownames(reference)
#'   if (is.null(expected_fixed)) {exp_list = c()} else {exp_list = rownames(expected_fixed)}
#'   if (is.null(inferred_fixed)) {inf_list = c()} else {inf_list = rownames(inferred_fixed)}
#'   if (is.null(input)) {input_list = c()} else {input_list = rownames(input)}
#'
#'   TP = length(intersect(inf_list, exp_list))
#'   FP = length(setdiff(inf_list, exp_list))
#'   #TN = length( setdiff( setdiff(ref_list, exp_list), inf_list) )
#'   TN = length( setdiff(setdiff(input_list, exp_list), inf_list)  )
#'   FN = length(setdiff(exp_list, inf_list))
#'
#'   accuracy = list(TP=TP, FP=FP, TN=TN, FN=FN)
#'
#'   #accuracy = (TP + TN) / (TP + TN + FP + FN)
#'   return(accuracy)
#' }
#'
#'
#' reconstruct.count = function(m, alpha, beta) {
#'   # all args are data.frame
#'   theta = diag(rowSums(m))               # matrix
#'   alpha = theta %*% as.matrix(alpha)     # matrix
#'   beta = as.matrix(beta)                 # matrix
#'
#'   mr_matrix = alpha %*% beta
#'   mr = round(as.data.frame(mr_matrix))
#'   rownames(mr) = rownames(m)
#'   return(mr)
#' }
#'
#'
#'
#' denovo.similarity = function(expected_denovo, inferred_denovo) {
#'
#'   if (length(expected_denovo)==0 | length(inferred_denovo)==0) {
#'     return(NULL)
#'   } else {
#'     df = data.frame(matrix(nrow = nrow(inferred_denovo), ncol = nrow(expected_denovo)))
#'     colnames(df) = rownames(expected_denovo)
#'     rownames(df) = rownames(inferred_denovo)
#'
#'     for (i in 1:nrow(inferred_denovo)) {
#'       inferred = inferred_denovo[i,]
#'       inferred_name = rownames(inferred)
#'       for (j in 1:nrow(expected_denovo)) {
#'         target = expected_denovo[j, ]
#'         target_name = rownames(target)
#'         score = cosine.vector(inferred, target)
#'         df[inferred_name, target_name] = score
#'       }
#'     }
#'
#'     match_df = data.frame(matrix(nrow = nrow(inferred_denovo), ncol = 2))
#'     colnames(match_df) = c("match", "similarity")
#'     rownames(match_df) = rownames(inferred_denovo)
#'
#'     similarity = 0
#'     iter = min(nrow(inferred_denovo), nrow(expected_denovo))
#'     for (i in 1:iter) {
#'
#'       max = which(df == max(df), arr.ind = TRUE)
#'       similarity = similarity + df[max]
#'
#'       row = row.names(df[max[,1],])
#'       column = names(df[max[,2]])
#'
#'       match_df[row, 'match'] = column
#'       match_df[row, 'similarity'] = df[max]
#'
#'       df[row, column] = 0
#'     }
#'
#'     return( list( similarity_average=(similarity / iter), match_df=match_df ) )
#'   }
#' }
#'
#'
#' denovo.ratio = function(expected_denovo, inferred_denovo) {
#'
#'   if (is.null(expected_denovo)) {n_exp = 0} else {n_exp = nrow(expected_denovo)}
#'   if (is.null(inferred_denovo)) {n_inf = 0} else {n_inf = nrow(inferred_denovo)}
#'
#'   denovo_ratio = (n_inf + 1) / (n_exp + 1)
#'
#'   return(denovo_ratio)
#' }
#'
