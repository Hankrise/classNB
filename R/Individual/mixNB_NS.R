#' Perform inference of mixture of negative binomial models using numerical solution
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of initial value of mean (optional).
#' @param var a vector of initial value of variance (optional).
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_NS algorithm (optional).
#' @param EPSILON a value for the threshold used for the stopping criteria for the mixNB_NS algorithm (optional).
#' @return A list of 3 objets.
#' \describe{
#' \item{\code{parameter}}{ a list with the inferred parameters: Size, Probability, mean, variance, weight.}
#' \item{\code{emisPr}}{a matrix with emission probability associated to each states in column and each position in row.}
#' \item{\code{postPr}}{a matrix with posterior probability associated to each series in column and each position in row.}
#' }
#' @export
#' @examples
#' data(toyexample)
#' # Perform inference of mixture of negative binomial models using numerical solution
#' resmixNB_NS <- mixNB_NS(X = toydata, nbS = 3)
#' @seealso \code{\link{mixNB_CF}}, \code{\link{mixNB_CF}}

mixNB_NS <- function(X, nbS, avg = NULL, var = NULL, iterMax = 15000, EPSILON = 1e-7) {
  eWeight <- rep(1, nbS) / nbS # the weight in NB mixture model
  if (is.null(c(avg, var))) {
    # if initial value of mean & variance are not given, use kmeans algorithm to generate
    esAvg <- 1
    esVar <- 0 # ensure mean < variance
    while (any(esAvg > esVar)) {
      X <- as.matrix(X)
      clust <- kmeans(x = as.vector(t(X)), centers = nbS) # apply kmeans
      esAvg <- as.vector(clust$centers) # calculate initial mean
      ordAvg <- order(esAvg) # the order of initial mean
      esAvg <- esAvg[ordAvg] # rearrange mean in increasing order
      esVar <- clust$withinss / (clust$size - 1) # calculate initial variance
      esVar <- esVar[ordAvg] # rearrange variance in increasing order of mean
      eSize <- esAvg^2 / (esVar - esAvg)
    }
  } else {
    # if initial value is given, use them.
    esAvg <- avg
    esVar <- var
    eSize <- esAvg^2 / (esVar - esAvg)
  }

  for (iter in 1:iterMax) {
    prevParam <- c(esAvg, eSize, eWeight)
    # cat("iter=",iter,"\n")
    ## Estep ————
    emisPr <- t(sapply(1:nbS, function(s) dnbinom(X, size = eSize[s], mu = esAvg[s])))
    rawPost <- emisPr * eWeight
    postPr <- t(t(rawPost) / colSums(rawPost))

    ## Mstep ————
    ## BFGS
    for (k in 1:nbS) {
      par <- c(eSize[k], esAvg[k])
      objFunc <- function(par) {
        -sum(dnbinom(X, size = par[1], mu = par[2], log = TRUE) * postPr[k, ])
      }
      devFunc <- function(par) {
        -as.matrix(c(sum(postPr[k, ] * (digamma(X + par[1]) - digamma(par[1] + log(par[2])))), sum(postPr[k, ] * (par[1] / par[2] - X / (1 - par[2])))))
      }

      eParam <- optim(par, objFunc, devFunc, method = "BFGS")$par

      eSize[k] <- eParam[1]
      esAvg[k] <- eParam[2]
    }

    sumWeight <- rowSums(postPr)
    eWeight <- sumWeight / sum(sumWeight)

    ordMean <- order(esAvg)
    esAvg <- sort(esAvg)
    eSize <- eSize[ordMean]
    eWeight <- eWeight[ordMean]

    newParam <- c(esAvg, eSize, eWeight)
    if (any(is.na(esAvg)) || any(is.na(eSize))) {
      warning("mixNB_NB_NS - NA is appeared!")
      esAvg <- prevParam[1:nbS]
      eSize <- prevParam[(nbS + 1):(2 * nbS)]
      break
    }

    D <- sum((newParam - prevParam)^2)
    if (D < EPSILON) {
      break
    }
  }
  cat("iterStop = ", iter, "\n")
  eProb <- eSize / (eSize + esAvg)
  esVar <- (esAvg^2 + eSize * esAvg) / eSize
  parameter <- list()
  parameter[["eSize"]] <- eSize
  parameter[["eProb"]] <- eProb
  parameter[["esAvg"]] <- esAvg
  parameter[["esVar"]] <- esVar
  parameter[["eWeight"]] <- eWeight
  structure(list(parameter = parameter, emisPr = emisPr, postPr = postPr), class = "mixNB_NS")
}
