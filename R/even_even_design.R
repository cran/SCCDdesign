#' Constructs an orthogonal design for two types of even-level factors using
#' a Hadamard matrix. The first column of the Hadamard matrix is removed,
#' and scaled copies of the resulting matrix are combined to obtain the
#' required factor levels.
#' @param a Number of levels for the first factor. Must be even.
#' @param b Number of levels for the second factor. Must be even.
#' @param m Order of the Hadamard matrix.Must be a multiple of 4. The constructed design contains
#'   `m - 1` factors of each of the two factor types.
#' @return A matrix containing the constructed design. The first
#'   `m - 1` columns correspond to the first type factor, and the remaining
#'   `m - 1` columns correspond to the second type factor.
#' @keywords internal
even_even_design <- function(a, b, m) {
  if (length(a) != 1 || !is.numeric(a) || a <= 0 ||
      a != floor(a)) {
    stop("a must be a positive integer.")
  }
  if (length(b) != 1 || !is.numeric(b) || b <= 0 ||
      b != floor(b)) {
    stop("b must be a positive integer.")
  }
  if (length(m) != 1 || !is.numeric(m) || m <= 0 ||
      m != floor(m)) {
    stop("m must be a positive integer.")
  }
  if (a %% 2 != 0 || b %% 2 != 0) {
    stop("a and b must both be even.")
  }
  k1 <- a / 2
  k2 <- b / 2
  if (k1 <= k2) {
    stop("a/2 must be greater than b/2.")
  }
  if (k1 %% k2 != 0) {
    stop("Construction is not feasible: b/2 must be a factor of a/2.")
  }
  if (m %% 4 != 0) {
    stop("m must be a multiple of 4.")
  }
  H <- hadamardMatrix(m)
  Hstar <- H[, -1, drop = FALSE]
  levels1 <- seq(-(2 * k1 - 1),
                 (2 * k1 - 1),
                 by = 2)
  positive1 <- levels1[levels1 > 0]
  H1 <- lapply(positive1, function(x) {
    x * Hstar
  })
  X1 <- do.call(rbind, H1)
  X1_foldover <- rbind(X1, -X1)
  levels2 <- seq(-(2 * k2 - 1),
                 (2 * k2 - 1),
                 by = 2)

  positive2 <- levels2[levels2 > 0]

  H2 <- lapply(positive2, function(x) {
    x * Hstar
  })
  X2 <- do.call(rbind, H2)
  n_matrix_1 <- nrow(X1_foldover)
  n_matrix_2 <- nrow(X2)
  n_rep <- n_matrix_1 / n_matrix_2

  if (n_rep != floor(n_rep)) {
    stop("Equal number of matrices cannot be obtained by replication.")
  }
  X2_replicate <- X2[
    rep(seq_len(n_matrix_2), times = n_rep),
    ,
    drop = FALSE
  ]
  design <- cbind(X1_foldover, X2_replicate)

  return(design)
}
