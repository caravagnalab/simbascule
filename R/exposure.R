library(gtools)
library(fdrtool)

generate.exposure <- function(beta, groups, private_sigs, private_fracs, seed=NULL) {

  signatures <- rownames(beta)
  if (!('SBS1' %in% signatures)) {
    stop('Wrong signatures! SBS1 not included!')
  }

  if (!is.null(seed)) {
    set.seed(seed = seed)
  }

  if (length(signatures) < 2) {
    stop("not valid! there are not enough signatures!")
  }

  df_list <- list()

  private_sigs.all = c(private_sigs$common, private_sigs$rare)
  shared_sigs = setdiff(signatures, private_sigs.all)

  unq = unique(groups)

  data = data.frame(matrix(ncol=length(signatures)+1, nrow=0))
  colnames(data) = c(signatures, "group")

  for (i in 1:length(unq)) {
    group = unq[i]

    if (length(unique(groups))==1) {
      sigNums <- length(signatures)
      sigNames <- signatures
    } else {

      sigNames_priv = c()

      # number of sigs to assign to "group"
      if (length(private_sigs.all)>0) {
        if (i == length(unq))
          sigNums_priv = length(private_sigs.all) else
          sigNums_priv = base::sample(1:round(length(private_sigs.all)/2), 1)

        sigNames_priv = base::sample(private_sigs.all, sigNums_priv)
        private_sigs.all = setdiff(private_sigs.all, sigNames_priv)

      sigNames = c(shared_sigs, sigNames_priv)
      sigNums = length(sigNames)
      }
    }

    num_samples <- length(groups[groups==group])

    # TEST - start
    mean_prior = fdrtool::rhalfnorm(sigNums, 1) %>% setNames(sigNames)
    alpha = sapply(1:sigNums, function(s)
      sapply(1:num_samples, function(x) sample_positive_norm(mean=mean_prior[s], sd=1))
      )
    alpha = apply(alpha, 1, function(x) x/sum(x)) %>% t() %>% as.data.frame()
    colnames(alpha) <- sigNames
    # TEST - end

    ## change values for rare/common

    alpha = adjust_frequency(alpha,
                     columns=intersect(private_sigs$rare, colnames(alpha)),
                     frac=private_fracs$rare, check=">=", mean=rep(0,length.out=sigNums) %>% setNames(sigNames), sd=0)
    alpha = adjust_frequency(alpha,
                     columns=intersect(private_sigs$common, colnames(alpha)),
                     frac=private_fracs$common, check="<=", mean=mean_prior, sd=1)

    # alpha = alpha / rowSums(alpha). ## check if it works
    alpha$group <- rep(group, num_samples)

    data = data %>%
      dplyr::add_row(alpha)

  }

  # merge all different group exposure matrices
  # data <- Reduce(function(x, y) { print(c(x,y)); merge(x, y, all=TRUE) }, df_list)

  # sort columns
  column_names <- colnames(data)
  #column_names <- column_names[order(column_names)]
  column_names <- append(setdiff(column_names, "group"), "group")
  #column_names[length(column_names)+1] <- "group"
  data <- data[, column_names]

  data[is.na(data)] <- 0    # convert 'NA' to zero
  data[order(data$group), ] # sort rows by group column

  return(data)
  }


adjust_frequency = function(alpha, columns, frac, check, mean, sd) {
  nn = round(nrow(alpha) * frac)
  for (colname in columns) {
    if (!match.fun(check)(sum(alpha[,colname]>0), nn) ) next
    samples.tmp = sample(1:nrow(alpha), size=nn)
    alpha[-samples.tmp, colname] = sample_positive_norm(mean[colname], sd)
  }

  return(alpha)
}


sample_positive_norm = function(mean, sd) {
  while(TRUE) {
    n = rnorm(1, mean, sd)
    if (n >= 0)
      return(n)
  }
}
