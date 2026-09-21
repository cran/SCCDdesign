#' Constructs an orthogonal design for two types of odd-level factors using
#' conference matrices.
#' @param a Number of levels for the first odd-level factor.
#' @param b Number of levels for the second odd-level factor.
#' @param m Order of the conference matrix. Must be a multiple of 4.The design contains
#'   `m` factors of the first odd-level type and `m` factors of the second
#'   odd-level type.
#' @return A matrix containing the constructed design. The first `m`
#'   columns correspond to the first odd-level factor and the remaining
#'   `m` columns correspond to the second odd-level factor. Both sets of
#'   columns are constructed using a conference matrix.
#' @keywords internal
odd_odd_design <- function(a, b, m) {
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
  if (a %% 2 == 0) {
    stop("a must be odd.")
  }
  if (b %% 2 == 0) {
    stop("b must be odd.")
  }
  if (m %% 4 != 0) {
    stop("m must be a multiple of 4.")
  }
  k1 <- (a - 1) / 2
  k2 <- (b - 1) / 2
  if (k1 <= k2) {
    stop("Construction is not feasible: k1 must be greater than k2.")
  }
  if (k1 %% k2 != 0) {
    stop("Construction is not feasible: k2 must be a factor of k1.")
  }
  C <- conferenceMatrix(m)
  levels1 <- if (k1 == 1) {
    1
  } else {
    seq(1, k1)
  }
  C1 <- lapply(levels1, function(x) {
    x * C
  })
  X1 <- do.call(rbind, C1)
  X1_foldover <- rbind(
    X1,
    -X1
  )
  levels2 <- if (k2 == 1) {
    1
  } else {
    seq(1, k2)
  }
  C2 <- lapply(levels2, function(x) {
    x * C
  })
  X2 <- do.call(rbind, C2)
  n_matrix_1 <- nrow(X1_foldover)
  n_matrix_2 <- nrow(X2)
  n_rep <- (n_matrix_1 / 2) / n_matrix_2
  if (n_rep != floor(n_rep)) {
    stop("Equal number of matrices cannot be obtained by replication.")
  }

  X2_positive_replicate <- X2[
    rep(seq_len(n_matrix_2), times = n_rep),
    ,
    drop = FALSE
  ]

  X2_replicate <- rbind(
    -X2_positive_replicate,
    X2_positive_replicate
  )

  design <- cbind(
    X1_foldover,
    X2_replicate
  )
  return(design)
}
