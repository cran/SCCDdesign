#' Constructs a D-optimal two-level design with \code{m} runs
#' and \code{m - 1} factors. The factor levels are coded as
#' -1 and +1.
#' @param m A positive integer specifying the number of runs.
#'   The minimum value is 3.
#' @return A numeric matrix containing the D-optimal two-level
#'   design.
#' @details
#' The candidate set is generated using all possible combinations
#' of -1 and +1 for the factors. The D-optimal design is selected
#' using \code{AlgDesign::optFederov()}.
#' @examples
#' dOptimalTwoLevel(4)
#' dOptimalTwoLevel(8)
#' @export
dOptimalTwoLevel <- function(m) {

  if (!is.numeric(m) ||
      length(m) != 1 ||
      is.na(m) ||
      !is.finite(m)) {

    stop("m must be a single finite numeric value.")
  }

  if (m != as.integer(m))
    stop("m must be an integer.")

  m <- as.integer(m)

  if (m < 3)
    stop("m must be at least 3.")

  k <- m - 1

  candidate <- expand.grid(
    rep(
      list(c(-1, 1)),
      k
    )
  )

candidate <- as.data.frame(candidate)

names(candidate) <- paste0(
    "Z",
    seq_len(k)
  )

  result <- AlgDesign::optFederov(
    ~ .,
    data = candidate,
    nTrials = m,
    approximate = FALSE
  )

  T <- as.matrix(
    result$design
  )

  storage.mode(T) <- "numeric"

  if ("(Intercept)" %in% colnames(T)) {

    T <- T[
      ,
      colnames(T) != "(Intercept)",
      drop = FALSE
    ]
  }

  if (
    nrow(T) != m ||
    ncol(T) != k
  ) {

    stop(
      "D-optimal search did not return the required dimensions."
    )
  }

  if (!all(T %in% c(-1, 1))) {

    stop(
      "D-optimal matrix must contain only -1 and +1."
    )
  }

  return(T)
}
