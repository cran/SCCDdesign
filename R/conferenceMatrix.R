#' Constructs a conference matrix of a specified order using a direct
#' construction for order 2 and Paley-based constructions for supported
#' odd prime or prime-power values of \eqn{q = n - 1}.
#' @param n A single finite integer specifying the order of the conference
#'   matrix. The value must be at least 2.
#' @return A numeric \eqn{n \times n} conference matrix.
#' @details
#' For \code{n = 2}, the conference matrix is constructed directly.
#' For larger orders, the function sets \eqn{q = n - 1} and determines
#' whether \eqn{q} is a supported odd prime or odd prime power.
#' Depending on the value of \eqn{q} modulo 4, the corresponding
#' Paley construction is used.
#' An error is returned when the supplied order does not satisfy the
#' requirements of the implemented construction.
#' @examples
#' conferenceMatrix(6)
#' @export
conferenceMatrix <- function(n) {
  if (missing(n))
    stop("Please specify the order n.")
  if (!is.numeric(n))
    stop("n must be numeric.")
  if (length(n) != 1)
    stop("n must be a single value.")
  if (is.na(n) || !is.finite(n))
    stop("n must be a finite value.")
  if (n != as.integer(n))
    stop("n must be an integer.")
  if (n < 2)
    stop("n must be at least 2.")
  n <- as.integer(n)
  if (n == 2) {
    C <- matrix(
      c(
        0, 1,
        1, 0
      ),
      nrow = 2,
      byrow = TRUE
    )

    return(C)
  }

  q <- n - 1

  info <- primePowerInfo(q)

  if (is.null(info)) {

    stop(
      paste0(
        "q = n - 1 = ",
        q,
        " is neither a prime nor a prime power supported ",
        "by the Paley construction implemented in this package."
      )
    )
  }

  if (q %% 2 == 0) {

    stop(
      paste0(
        "q = n - 1 = ",
        q,
        " is even. The implemented Paley construction ",
        "requires an odd prime or odd prime power."
      )
    )
  }

  if (q %% 4 == 1) {

    paleyType <- "Paley Type II"

  } else if (q %% 4 == 3) {

    paleyType <- "Paley Type I"

  } else {

    stop("Paley construction is not applicable.")
  }

  if (info$m == 1) {

    C <- paleyPrime(n)

  } else {

    # Check that GF(p^m) is implemented
    irreduciblePolynomial(
      info$p,
      info$m
    )

    C <- paleyPrimePower(n)
  }

  return(C)
}
