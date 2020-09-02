# Function to return the id of the first selected variable
etree_sel <- function(response,
                      covariates,
                      case.weights = NULL,
                      minbucket = 5,
                      alpha = 0.05,
                      R = 1000,
                      split.type = 'coeff',
                      coef.split.type = 'test',
                      p.adjust.method = 'fdr',
                      nb = 5) {

  # Check whether covariates is a list
  if(!is.list(covariates)) stop("Argument 'covariates' must be provided as a list")

  # Number of covariates
  n.var = length(covariates)

  # If the case weights are not provided, they are all initialized as 1
  if(is.null(case.weights))
    case.weights <- rep(1L, as.numeric(length(response)))

  # New list of covariates (needed here to build the df used by party)
  newcovariates = lapply(covariates, function(j){

    if(class(j) == 'fdata'){

      if(split.type == "coeff"){

        foo <- fda.usc::optim.basis(j, numbasis = nb)
        fd3 <- fda.usc::fdata2fd(foo$fdata.est,
                                 type.basis = "bspline",
                                 nbasis = foo$numbasis.opt)
        foo <- t(fd3$coefs)

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      attr(foo, 'cov.type') <- 'fdata'
      return(foo)

    } else if(class(j) == 'list' & all(sapply(j, function(x) attributes(x)$names) == 'diagram')){

      foo <- as.factor(1:length(response))
      attr(foo, 'cov.type') <- 'diagram'
      return(foo)

    } else if(class(j) == 'list' &
              all(sapply(j, class) == 'igraph')){

      if(split.type == "coeff"){

        foo <- graph.shell(j)

      } else if(split.type == "cluster"){

        foo <- as.factor(1:length(response))

      }

      attr(foo, 'cov.type') <- 'graph'
      return(foo)

    } else {

      return(j)

    }
  })

  # Covariates name
  if(!is.null(names(covariates))){
    names(newcovariates) <- names(covariates)
    #control if any name is void, i.e. if it is ''
    no_name <- which(sapply(names(covariates), function(n) n == '', USE.NAMES = FALSE))
    names(newcovariates) <- replace(names(newcovariates),
                                    no_name,
                                    as.factor(1:length(no_name)))
  } else {
    warning('No names available for covariates. Numbers are used instead.')
    names(newcovariates) <- 1:length(newcovariates)
  }

  # Distances
  cov.distance <- lapply(covariates, compute.dissimilarity)

  # Large list with covariates, newcovariates and distances
  covariates.large = list('cov' = covariates, 'newcov' = newcovariates, 'dist' = cov.distance)

  # Number of original covariates
  n.cov = length(covariates.large$cov)

  # Independence test between the response and each covariate
  p = lapply(covariates.large$dist,
             function(cov.dist) {
               ct <- energy::dcor.test(cov.dist,
                                       compute.dissimilarity(response),
                                       R = R)
               if (!is.na(ct$statistic)) {
                 return(c(ct$statistic, ct$p.value))
               } else{
                 c(NA, NA)
               }
             }
  )

  # Create matrix with test stats and p-values; check if p-values are all NULL
  p = t(matrix(unlist(p), ncol = 2, byrow = T))
  rownames(p) <- c("statistic", "p-value")
  if (all(is.na(p[2,]))) return(NULL)

  # Multiple testing correction
  adj_p <- p.adjust(p[2,], method = p.adjust.method)

  # Stop criterion
  if (min(adj_p) > alpha) return(NULL)

  # Variable selection (based on original p-values)
  if (length(which(p[2,] == min(p[2,], na.rm = T))) > 1) {
    xselect <- which.max(p[1,])
    #in case of multiple minima, take that with the highest test statistic
  } else{
    xselect <- which.min(p[2,])
  }

  return(xselect)
}
