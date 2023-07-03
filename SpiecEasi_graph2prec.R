
library(SpiecEasi)


graph2prec <- function(Graph, posThetaLims=c(2,3), negThetaLims=-posThetaLims, targetCondition=100, epsBin=1e-2,
                        numBinSearch=100) {

    if (class(Graph) != 'graph') stop('input is not a graph')

    Graph = SpiecEasi::make_graph(D = 20, e = 10, method = "erdos_renyi")

    n <- ncol(Graph)


    posThetaLims <- sort(posThetaLims)
    negThetaLims <- sort(negThetaLims)

    #add diagonal
    Theta <- Graph + diag(n)
    # node degree
    degVec <- colSums(Graph)

    # add random edge weights uniformly from theta_min to theta_max
    # this selects the upper triangular matrix without the diagonal
    utri  <- SpiecEasi::triu(Theta, k=1)
    # non-zero entries
    nzind <- which(utri != 0)

    if (length(posThetaLims) > 2 || length(negThetaLims) > 2)
        stop("theta_max and theta_min should be a numeric vector of length 1 or 2")

    # generate random numbers, those are on the range 0, 1 now
    rands <- runif(length(nzind))

    # map to range
    mapToRange <- function(x, lim) {
        span <- diff(sort(lim))
        min(lim) + (x * span)
    }

    # randomly assign positive and negative values
    boolind <- sample(c(TRUE, FALSE), length(nzind), replace=TRUE)


    # wird dann auf daen positiven oder den negativen theta support gemappt
    rands[boolind]  <- mapToRange(rands[boolind], posThetaLims)
    rands[!boolind] <- mapToRange(rands[!boolind], negThetaLims)


    utri[nzind] <- rands
    #da wird dann eine symmetrische matrix draus
    Theta <- triu2diag(utri, 1)

    # find smallest eigenvalue such that Theta is invertible
    eigVals <- eigen(Theta)$values
    minEig  <- min(eigVals)
    maxEig  <- max(eigVals)

    # if smallest eigenvalue is negative, add the absolute value to the diagonal, otherwise do nithing (???)
    if (minEig < 1e-2) Theta <- Theta + abs(minEig)*diag(n)

    # hier muss ich noch verstehen was genau passiert
    diagConst <- .binSearchCond(Theta, targetCondition, numBinSearch, epsBin)
    Theta <- Theta + diagConst*diag(n)
    return(Theta)
}



#' @noRd
.binSearchCond <- function(Theta, condTheta, numBinSearch, epsBin) {
# Internal function that determines the constant in the diagonal to satisfy the
# condition constraint on the Precision/Covariance matrix

    n <- nrow(Theta)
    currCondTheta <- kappa(Theta)
    if (currCondTheta < condTheta) {
        # Max entry in the diagonal (lower bound)
        currLB   <- -max(diag(Theta))
        stepSize <- currLB+.Machine$double.eps

        while (currCondTheta < condTheta) {
            currCondTheta <- kappa(Theta+stepSize*diag(n))
            stepSize      <- stepSize/2
        }
        currUB <- stepSize
    } else {
        currLB <- 0
        stepSize = 0.1

        while (currCondTheta > condTheta) {
            currCondTheta <- kappa(Theta + stepSize*diag(n))
            stepSize      <- 2*stepSize
        }
        currUB <- stepSize
    }

    for (i in 1:numBinSearch) {
        diagConst <- (currUB+currLB)/2
        currCondTheta <- kappa(Theta+diagConst*diag(n))

        if (currCondTheta < condTheta) currUB <- diagConst
        else currLB <- diagConst

        if (abs(currCondTheta-condTheta)<epsBin) break
    }
    diagConst
}
