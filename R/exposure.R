library(gtools)

generate.exposure <- function(beta, groups, seed=NULL) {

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

  signatures <- signatures[! signatures %in% c('SBS1')] # excludes SBS1

  for (group in unique(groups)) {

    if (length(unique(groups))==1) {
      sigNums <- length(signatures)
      sigNames <- c('SBS1', signatures)
    } else {
      sigNums <- sample(1:length(signatures), 1)
      sigNames <- c('SBS1', sample(signatures, sigNums))
    }
    #sigNums <- sample(1:length(signatures), 1)
    #sigNames <- c('SBS1', sample(signatures, sigNums))

    num_samples <- length(groups[groups==group])

    #print(paste("group", group, "has", sigNums+1, "signatures, and", num_samples, "samples"))

    # TEST - start
    alpha <- data.frame(matrix(ncol = sigNums+1, nrow = 0))
    for (i in 1:num_samples) {
      s <- rdirichlet(1, alpha = sample(1:100, sigNums+1, replace=FALSE))
      alpha <- rbind(alpha, c(s))
    }
    colnames(alpha) <- sigNames
    alpha$group <- rep(group, num_samples)
    # TEST - end


    #x <- matrix( runif(num_samples * (sigNums+1), 0, 1), ncol = sigNums+1 )
    #alpha <- x / rowSums(x)
    #alpha <- as.data.frame(alpha)
    #colnames(alpha) <- sigNames
    #alpha$group <- rep(group, num_samples)
    #print(alpha)


    df_list[length(df_list)+1] <- list(alpha)
  }

  # merge all different group exposure matrices
  data <- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

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
