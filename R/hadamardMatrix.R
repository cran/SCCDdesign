#' Constructs a Hadamard matrix of order \code{n} using
#' \code{pracma::hadamard()}.
#' @param n A positive integer specifying the order of the matrix.
#' @return A Hadamard matrix of order \code{n}.
#' @details
#' A Hadamard matrix contains only +1 and -1, and its rows
#' are mutually orthogonal. The supported orders are of the
#' forms \eqn{2^e}, \eqn{12(2^e)}, and \eqn{20(2^e)}.
#' @examples
#' hadamardMatrix(4)
#' hadamardMatrix(8)
#' hadamardMatrix(12)
#' @export
hadamardMatrix <- function(n) {

  if (missing(n)) {
    stop("Please specify the order n.")
  }
  if (!is.numeric(n) ||
      length(n) != 1L ||
      is.na(n) ||
      !is.finite(n)) {
    stop("n must be a single finite numeric value.")
  }
   if (n != as.integer(n)) {
    stop("n must be an integer.")
  }
  n <- as.integer(n)
   if (n <= 1L) {
    stop("n must be greater than 1.")
  }
# Construct using pracma
  H <- tryCatch(
    pracma::hadamard(n),
    error = function(e) {
      stop(
        paste0(
          "A Hadamard matrix of order ", n,
          " cannot be constructed using pracma::hadamard()."
        ),
        call. = FALSE
      )
    }
  )
  if (!isTRUE(
    all.equal(
      H %*% t(H),
      n * diag(n),
      tolerance = 1e-10
    )
  )) {
    stop("The constructed matrix failed Hadamard verification.")
  }
  H
}
