#' Performs a simulation study to evaluate the ability of a design
#' to identify active main effects using best-subsets regression
#' and the AICc criterion.
#' @param X A design matrix containing the mixed-level factors.
#' @param SNR A numeric value specifying the signal-to-noise ratio.
#' @param n.active An integer specifying the number of active factors.
#' @param nrep An integer specifying the number of simulation replications.
#' Default is 1000.
#' @return A data frame containing the SNR, number of active factors,
#' power, and Type I error rate.
#' @importFrom leaps regsubsets
#' @importFrom stats rexp rnorm coef
#' @examples
#' \dontrun{
#' X <- matrix(rnorm(20 * 5), nrow = 20, ncol = 5)
#' simulation(X, SNR = 1, n.active = 2, nrep = 100)
#' }
#' @export
simulation <- function(X, SNR, n.active, nrep = 1000) {
  X <- as.matrix(X)
  if (!is.numeric(X)) {
    stop("X must be a numeric design matrix.")
  }
  p <- ncol(X)
  if (n.active < 1 || n.active > p) {
    stop("n.active must be between 1 and the number of columns of X.")
  }
  if (length(SNR) != 1 || !is.numeric(SNR)) {
    stop("SNR must be a single numeric value.")
  }
  if (length(nrep) != 1 || nrep < 1 || nrep != as.integer(nrep)) {
    stop("nrep must be a positive integer.")
  }
  quad <- X^2
  colnames(quad) <- paste0(colnames(X), "_Q")
  TP.total <- 0
  FP.total <- 0
  Inactive.total <- 0
  for (rep in seq_len(nrep)) {
    active <- sample(seq_len(p), n.active)
    beta <- rep(0, p)
    beta[active] <-
      (SNR - log(2) + rexp(n.active, rate = 1)) *
      sample(c(-1, 1), n.active, replace = TRUE)
    y <- as.vector(X %*% beta + rnorm(nrow(X)))
    dat <- data.frame(y = y, X, quad)
    fit <- leaps::regsubsets(
      y ~ .,
      data = dat,
      nvmax = ncol(dat) - 1
    )
    s <- summary(fit)
    rss <- s$rss
    k <- sapply(
      seq_along(rss),
      function(i) length(coef(fit, i))
    )
    n <- nrow(dat)
    AIC <- n * log(rss / n) + 2 * k
    AICc <- AIC +
      (2 * k * (k + 1)) / (n - k - 1)
    best <- which.min(AICc)
    coef.best <- names(coef(fit, best))
    selected <- match(coef.best, colnames(X))
    selected <- selected[!is.na(selected)]
    TP <- length(intersect(active, selected))
    inactive <- setdiff(seq_len(p), active)
    FP <- length(intersect(inactive, selected))
    TP.total <- TP.total + TP
    FP.total <- FP.total + FP
    Inactive.total <- Inactive.total + length(inactive)
  }
  Power <- TP.total / (nrep * n.active)
  Type1 <- FP.total / Inactive.total
  results <- data.frame(
    SNR = SNR,
    Active = n.active,
    Power = Power,
    Type1 = Type1
  )
  return(results)
}
