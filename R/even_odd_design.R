#' Constructs an orthogonal design for an odd-level factor and an
#' even-level factor using a conference matrix and a Hadamard matrix,
#' respectively.
#' @param a Number of levels for the odd-level factor.
#' @param b Number of levels for the even-level factor.
#' @param m Order of the conference and Hadamard matrices. Must be a
#'   multiple of 4.The design contains`m` factors of the odd-level type and `m - 1` factors of the
#'   even-level type.
#' @return A matrix containing the constructed design. The first `m`
#'   columns correspond to the odd-level factors and are constructed
#'   using a conference matrix, while the remaining `m - 1` columns
#'   correspond to the even-level factors and are constructed using
#'   a Hadamard matrix with its first column removed.
#' @keywords internal
odd_even_design <- function(a, b, m) {
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
  if (b %% 2 != 0) {
    stop("b must be even.")
  }
  if (m %% 4 != 0) {
    stop("m must be a multiple of 4.")
  }
  k1 <- (a - 1) / 2
  k2 <- b / 2
  C <- conferenceMatrix(m)
  H <- hadamardMatrix(m)
  Hstar <- H[, -1, drop = FALSE]
  if (k1 >= k2) {
    if (k1 %% k2 != 0) {
      stop("Construction is not feasible: k2 must be a factor of k1.")
    }
    levels1 <- if (k1 == 1) {
      1
    } else {
      seq(1, k1)
    }
    C1 <- lapply(levels1, function(x) {
      x * C
    })
    X1 <- do.call(rbind, C1)
    X1_foldover <- rbind(X1, -X1)
    levels2 <- seq(1, 2 * k2 - 1, by = 2)
    H2 <- lapply(levels2, function(x) {
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
   design <- cbind(
      X1_foldover,
      X2_replicate
    )
  } else {
    if (k2 %% k1 != 0) {
      stop("Construction is not feasible: k1 must be a factor of k2.")
    }
    levels2 <- seq(1, 2 * k2 - 1, by = 2)
    H2 <- lapply(levels2, function(x) {
      x * Hstar
    })
    X2 <- do.call(rbind, H2)
    n_rep_H <- k2 / k1
    X2_replicate <- X2[
      rep(seq_len(nrow(X2)), times = n_rep_H),
      ,
      drop = FALSE
    ]
    levels1 <- if (k1 == 1) {
      1
    } else {
      seq(1, k1, by = 1)
    }
    C1 <- lapply(levels1, function(x) {
      x * C
    })
    X1_positive <- do.call(rbind, C1)
    half_hadamard <- nrow(X2_replicate) / 2
    n_rep_C <- half_hadamard / nrow(X1_positive)
    if (n_rep_C != floor(n_rep_C)) {
      stop("Conference component cannot be replicated to the required size.")
    }
    X1_replicate <- X1_positive[
      rep(seq_len(nrow(X1_positive)), times = n_rep_C),
      ,
      drop = FALSE
    ]
    X1 <- rbind(
      X1_replicate,
      -X1_replicate
    )

    design <- cbind(
      X1,
      X2_replicate
    )
  }

  return(design)
}
