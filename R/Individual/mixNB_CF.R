#' Perform inference of mixture of negative binomial models using closed form algorithm
#' @param X a vector of observations.
#' @param nbS a integer specifying the numbers of states.
#' @param avg a vector of initial value of mean (optional).
#' @param var a vector of initial value of variance (optional).
#' @param iterMax an integer specifying the maximal number of iterations for the mixNB_CF algorithm (optional).
#' @param EPSILON a value for the threshold used for the stopping criteria for the mixNB_CF algorithm (optional).
#' @return A list of 3 objets.
#' \describe{
#' \item{\code{parameter}}{ a list with the inferred parameters: Size, Probability, mean, variance, weight.}
#' \item{\code{emisPr}}{a matrix with emission probability associated to each states in column and each position in row.}
#' \item{\code{postPr}}{a matrix with posterior probability associated to each series in column and each position in row.}
#' }
#' @export
#' @examples
#' data(toyexample)
#' # Perform inference of mixture of negative binomial models using closed form algorithm
#' resmixNB_CF <- mixNB_CF(X = toydata, nbS = 3)
#' @seealso \code{\link{mixNB_NS}}, \code{\link{mixNB_NS}}
mixNB_CF <- function(X, nbS, avg = NULL, var = NULL, iterMax = 15000, EPSILON = 1e-7) {
  eWeight <- rep(1, nbS) / nbS # the initial value of weight
  nbT <- length(X) # amount of samples

  if (is.null(c(avg, var))) {
    # if initial value of mean & variance are not given, use kmeans algorithm to generate
    esAvg <- 1
    esVar <- 0 # ensure mean < variance
    while (any(esAvg > esVar)) {
      clust <- kmeans(x = as.vector(t(X)), centers = nbS) # apply kmeans
      esAvg <- as.vector(clust$centers) # calculate initial mean
      ordAvg <- order(esAvg) # the order of initial mean
      esAvg <- esAvg[ordAvg] # rearrange mean in increasing order
      esVar <- clust$withinss / (clust$size - 1) # calculate initial variance
      esVar <- esVar[ordAvg] # rearrange variance in increasing order of mean
    }
  } else {
    # if initial value is given, use them.
    esAvg <- avg
    esVar <- var
  }
  eSize <- esAvg^2 / (esVar - esAvg) # calculate initial value of size
  eProb <- esAvg / esVar # calculate initial value of probability
  ones_nbS <- as.matrix(rep(1, nbS)) # a vector with nbS elements equals 1
  ones_nbT <- as.matrix(rep(1, nbT)) # a vector with nbT elements equals 1
  ones_nbTS <- ones_nbT %*% t(ones_nbS) # a matrix with nbT*nbS elements equals 1
  ones_nbST <- ones_nbS %*% t(ones_nbT) # a matrix with nbS*nbT elements equals 1
  ones_nbSS <- ones_nbS %*% t(ones_nbS) # a matrix with nbS*nbS elements equals 1
  # iteration starts
  for (iter in 1:iterMax) {
    prevParam <- c(esAvg, eSize, eWeight) # define iterative parameters in an arrray
    Nmat <- eSize %*% t(ones_nbT) # matrix of size
    Mmat <- esAvg %*% t(ones_nbT) # matrix of mean
    Xmat <- ones_nbS %*% t(X) # matrix of sample X
    ## Estep --------
    Pmat <- dnbinom(Xmat, size = Nmat, mu = Mmat) # matrix of emission Probability
    rawPost <- Pmat * eWeight
    postPr <- rawPost / (ones_nbSS %*% rawPost) # matrix of post Probability
    beta <- 1 - 1 / (1 - eProb) - 1 / log(eProb) # value of beta
    Bmat <- beta %*% t(ones_nbT) # matrix of beta
    Dmat <- Nmat * (digamma(Nmat + Xmat) - digamma(Nmat)) # matrix of delta
    ## Mstep ---------
    eProb <- beta * diag(postPr %*% t(Dmat)) / diag(postPr %*% t(Xmat - ((ones_nbST - Bmat) * Dmat))) # new value of probability
    eSize <- -diag(postPr %*% t(Dmat)) / diag(postPr %*% ones_nbTS) / log(eProb) # new value of size
    sumWeight <- diag(postPr %*% ones_nbTS)
    eWeight <- sumWeight / sum(sumWeight) # new value of weight
    esAvg <- eSize * (1 - eProb) / eProb # new value of mean
    esVar <- eSize * (1 - eProb) / eProb^2 # new value of variance
    ordAvg <- order(esAvg) # order of mean
    esAvg <- esAvg[ordAvg] # rearrange mean with increasing order
    eSize <- eSize[ordAvg] # rearrange size in increasing order of mean
    eWeight <- eWeight[ordAvg] # rearrange weight in increasing order of mean
    eProb <- eProb[ordAvg] # rearrange probability in increasing order of mean
    esVar <- esVar[ordAvg] # rearrange variance in increasing order of mean
    newParam <- c(esAvg, eSize, eWeight) # update iterative parameters in an arrray
    # throw warning when any element in the iterative parameters can't be calculated
    if (any(is.na(esAvg)) || any(is.na(esAvg))) {
      warning("mixNB_NB_CF - NA is appeared!")
      esAvg <- prevParam[1:nbS]
      eSize <- prevParam[(nbS + 1):(2 * nbS)]
      break
    }
    D <- sum((newParam - prevParam)^2) # when precision reach the desired level, stop iteration, otherwise reach the maximum limit
    if (D < EPSILON) {
      break
    }
  }
  cat("iterStop = ", iter, "\n")
  # return with a list
  parameter <- list()
  parameter[["eSize"]] <- eSize
  parameter[["eProb"]] <- eProb
  parameter[["esAvg"]] <- esAvg
  parameter[["esVar"]] <- esVar
  parameter[["eWeight"]] <- eWeight
  structure(list(parameter = parameter, emisPr = Pmat, postPr = postPr), class = "mixNB_CF")
}
