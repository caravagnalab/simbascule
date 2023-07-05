library(gtools)
library(fdrtool)

generate.exposure <- function(beta, groups, private_sigs, private_fracs,
                              n_priv_per_group=1, seed=NULL, thr=0.1) {

  if (is.null(n_priv_per_group)) n_priv_per_group = round(length(private_sigs.all)/2)

  signatures <- rownames(beta)
  if (!('SBS1' %in% signatures)) stop('Wrong signatures! SBS1 not included!')

  if (!is.null(seed)) set.seed(seed = seed)

  if (length(signatures) < 2) stop("not valid! there are not enough signatures!")

  df_list <- list()

  private_sigs.all = private_sigs.start = c(private_sigs$common, private_sigs$rare)
  shared_sigs = setdiff(signatures, private_sigs.all)

  unq = unique(groups) %>% sort()

  data = data.frame(matrix(ncol=length(signatures)+1, nrow=0))
  colnames(data) = c(signatures, "group")

  for (i in 1:length(unq)) {
    group = unq[i]
    # print(group)
    # print(private_sigs.all)

    if (length(unique(groups))==1) {
      sigNums <- length(signatures)
      sigNames <- signatures
    } else {

      sigNames_priv = c()

      # number of sigs to assign to "group"
      if (length(private_sigs.all)>0) {
        if (i == length(unq))
          sigNums_priv = length(private_sigs.all) else
          sigNums_priv = base::sample(1:n_priv_per_group, 1)

        sigNames_priv = base::sample(private_sigs.all, sigNums_priv)

        sigNames_priv = c(sigNames_priv, private_sigs.start[sigNames_priv],
                          private_sigs.start[private_sigs.start==sigNames_priv] %>% names) %>%
          purrr::discard(function(x) x=="" || is.na(x)) %>% unique()

        private_sigs.all = private_sigs.all %>% purrr::discard(function(x) x %in% sigNames_priv)

        sigNames = c(shared_sigs, sigNames_priv)
        sigNums = length(sigNames)
      }
    }

    num_samples <- length(groups[groups==group])

    mean_prior = fdrtool::rhalfnorm(sigNums, 1) %>% setNames(sigNames)
    alpha = sapply(1:sigNums, function(s)
      sapply(1:num_samples, function(x) sample_positive_norm(mean=mean_prior[s], sd=1))
      )
    alpha = apply(alpha, 1, function(x) x/sum(x)) %>% t() %>% as.data.frame()
    colnames(alpha) <- sigNames

    ## change values for rare/common

    alpha = adjust_frequency(alpha, N=length(groups),
                     columns=intersect(private_sigs$rare, colnames(alpha)),
                     frac=private_fracs$rare, check=">=",
                     # mean=rep(0, length.out=length(intersect(private_sigs$rare, colnames(alpha)))) %>%
                     #   setNames(intersect(private_sigs$rare, colnames(alpha))),
                     mean=rep(0, length.out=ncol(alpha)) %>% setNames(colnames(alpha)),
                     sd=0, min_exp=0.3)

    alpha = adjust_frequency(alpha, N=length(groups),
                     columns=intersect(private_sigs$common, colnames(alpha)),
                     frac=private_fracs$common, check="<",
                     mean=mean_prior, sd=1, thr=thr, min_exp=0.1)

    alpha = adjust_frequency(alpha, N=length(groups),
                             columns=shared_sigs,
                             frac=nrow(alpha), check="<",
                             mean=mean_prior,
                             sd=1, thr=thr, min_exp=0.1)

    alpha = apply(alpha, 1, function(x) x/sum(x)) %>% t() %>% as.data.frame()

    # alpha = alpha / rowSums(alpha). ## check if it works
    alpha$group <- group

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


adjust_frequency = function(alpha, N, columns, frac, check, mean, sd, thr=0., min_exp=0.2) {
  if (!is.integer(frac)) { nn = min(nrow(alpha), round(N * frac)) } else { nn = frac}
  for (colname in columns) {
    alpha = alpha %>% tibble::as_tibble() %>% dplyr::mutate(!!colname:=replace( .[[colname]], {{colname}}<thr, 0 ))

    if ( !match.fun(check)(sum(alpha[,colname]>0), nn) ) next
    samples.tmp = sample(1:nrow(alpha), size=nn)
    alpha[-samples.tmp, colname] = sample_positive_norm(mean[colname], sd)
    while (any(alpha[samples.tmp, colname] < min_exp))
      alpha[samples.tmp, colname] = sample_positive_norm(mean[colname], sd=1)
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
