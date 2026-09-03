#' Quadratic Residues Modulo a Prime
#' Computes the nonzero quadratic residues modulo a prime number.
#' @param q A prime number.
#' @return A vector containing the nonzero quadratic residues modulo q.
#' @keywords internal
quadraticResiduesPrime <- function(q) {

  if (!is_prime(q))
    stop("q must be prime.")

  QR <- unique(
    ((1:(q - 1))^2) %% q
  )

  QR <- QR[QR != 0]

  return(sort(QR))
}
#' Quadratic Character Modulo a Prime
#' Computes the quadratic character of an integer modulo a prime.
#' It returns 0 for zero, 1 for a quadratic residue, and -1 otherwise.
#' @param x An integer.
#' @param q A prime number.
#' @param QR A vector containing the nonzero quadratic residues modulo q.
#' @return The quadratic character: 0, 1, or -1.
#' @keywords internal
quadraticCharacterPrime <- function(x, q, QR) {

  x <- x %% q

  if (x == 0)
    return(0)

  if (x %in% QR)
    return(1)

  return(-1)
}
#' Jacobsthal Matrix for a Prime
#' Constructs a Jacobsthal matrix when q is prime using
#' the quadratic character modulo q.
#' @param q A prime number specifying the order of the
#'   Jacobsthal matrix.
#' @return A numeric Jacobsthal matrix of order q.
#' @keywords internal
jacobsthalPrime <- function(q) {

  if (!is_prime(q))
    stop("q must be prime.")

  QR <- quadraticResiduesPrime(q)

  J <- matrix(
    0,
    nrow = q,
    ncol = q
  )

  for (i in 0:(q - 1)) {

    for (j in 0:(q - 1)) {

      difference <- (i - j) %% q

      J[i + 1, j + 1] <-
        quadraticCharacterPrime(
          difference,
          q,
          QR
        )
    }
  }

  return(J)
}
#' Paley Conference Matrix for a Prime
#' Constructs a Paley conference matrix of order n when
#' q = n - 1 is an odd prime.
#' @param n A positive integer specifying the order of the
#'   conference matrix.
#' @return A numeric conference matrix of order n.
#' @keywords internal
# Construct a Paley conference matrix when q = n - 1 is prime
paleyPrime <- function(n) {

  q <- n - 1

  if (!is_prime(q))
    stop("q = n - 1 must be prime.")

  if (q %% 2 == 0)
    stop("q must be an odd prime.")

  J <- jacobsthalPrime(q)

  C <- matrix(
    0,
    nrow = n,
    ncol = n
  )

  # First row
  C[1, 2:n] <- 1

  # First column
  C[2:n, 1] <- 1

  # Jacobsthal block
  C[2:n, 2:n] <- J

  # Conference matrix has zero diagonal
  diag(C) <- 0

  return(C)
}
#' Jacobsthal Matrix for a Prime Power
#' Constructs a Jacobsthal matrix over the finite field GF(p^m)
#' when q = p^m is an odd prime power.
#' @param q An odd prime power specifying the order of the
#'   Jacobsthal matrix.
#' @return A numeric Jacobsthal matrix of order q.
#' @keywords internal
jacobsthalPrimePower <- function(q) {

  info <- primePowerInfo(q)

  if (is.null(info))
    stop("q must be a prime power.")

  p <- info$p
  m <- info$m

  if (m == 1)
    stop("q is prime. Use the prime construction.")

  if (p == 2)
    stop("This construction requires an odd prime power.")

  # Check that the finite field is implemented
  irreduciblePolynomial(p, m)

  # Generate GF(p^m)
  GF <- generateGF(p, m)

  # Find nonzero quadratic residues
  QR <- quadraticResiduesGF(p, m)

  J <- matrix(
    0,
    nrow = q,
    ncol = q
  )

  for (i in seq_len(q)) {

    for (j in seq_len(q)) {

      difference <- gfSubtract(
        GF[i, ],
        GF[j, ],
        p
      )

      J[i, j] <- quadraticCharacterGF(
        difference,
        QR
      )
    }
  }

  return(J)
}
#' Paley Conference Matrix for a Prime Power
#' Constructs a Paley conference matrix of order n when
#' q = n - 1 is an odd prime power.
#' @param n A positive integer specifying the order of the
#'   conference matrix.
#' @return A numeric conference matrix of order n.
#' @keywords internal
paleyPrimePower <- function(n) {

  q <- n - 1

  info <- primePowerInfo(q)

  if (is.null(info))
    stop("q = n - 1 must be a prime power.")

  if (info$m == 1)
    stop("q is prime. Use the prime construction.")

  if (info$p == 2)
    stop("This construction requires an odd prime power.")

  J <- jacobsthalPrimePower(q)

  C <- matrix(
    0,
    nrow = n,
    ncol = n
  )

  # First row
  C[1, 2:n] <- 1

  # First column
  C[2:n, 1] <- 1

  # Jacobsthal block
  C[2:n, 2:n] <- J

  # Zero diagonal
  diag(C) <- 0

  return(C)
}
