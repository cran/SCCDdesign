#' Irreducible Polynomial
#' Returns a predefined irreducible polynomial used for
#' constructing the finite field GF(p^m).
#' @param p A prime number.
#' @param m A positive integer specifying the field extension degree.
#' @return A numeric vector containing the polynomial coefficients.
#' @keywords internal
irreduciblePolynomial <- function(p, m) {

  # GF(3^2): x^2 + 1
  if (p == 3 && m == 2)
    return(c(1, 0, 1))

  # GF(5^2): x^2 + 2
  if (p == 5 && m == 2)
    return(c(2, 0, 1))

  # GF(3^3): x^3 + 2x + 1
  if (p == 3 && m == 3)
    return(c(1, 2, 0, 1))

  stop(
    paste0(
      "Finite-field construction for GF(",
      p, "^", m,
      ") is not currently implemented."
    )
  )
}
#' Convert an Integer to a Finite-Field Element
#' Converts an integer to its coefficient-vector representation
#' in GF(p^m).
#' @param x A non-negative integer.
#' @param p A prime number.
#' @param m A positive integer specifying the field extension degree.
#' @return An integer vector representing an element of GF(p^m).
#' @keywords internal
intToGF <- function(x, p, m) {

  coefficients <- integer(m)

  for (i in seq_len(m)) {

    coefficients[i] <- x %% p
    x <- x %/% p
  }

  return(coefficients)
}
#' Generate Elements of a Finite Field
#' Generates the coefficient-vector representations of all
#' elements of GF(p^m).
#' @param p A prime number.
#' @param m A positive integer specifying the field extension degree.
#' @return A matrix containing all elements of GF(p^m).
#' @keywords internal
generateGF <- function(p, m) {

  q <- p^m

  GF <- matrix(
    0,
    nrow = q,
    ncol = m
  )

  for (x in 0:(q - 1)) {

    GF[x + 1, ] <- intToGF(
      x,
      p,
      m
    )
  }

  return(GF)
}
#' Finite-Field Addition
#' Adds two finite-field elements componentwise modulo p.
#' @param a A finite-field element.
#' @param b A finite-field element.
#' @param p A prime number.
#' @return The sum of the two elements modulo p.
#' @keywords internal
gfAdd <- function(a, b, p) {

  return((a + b) %% p)
}
#' Finite-Field Subtractio
#' Subtracts two finite-field elements componentwise modulo p.
#' @param a A finite-field element.
#' @param b A finite-field element.
#' @param p A prime number.
#' @return The difference of the two elements modulo p.
#' @keywords internal
gfSubtract <- function(a, b, p) {

  return((a - b) %% p)
}
#' Polynomial Multiplication over GF(p
#' Multiplies two polynomials with coefficients in GF(p).
#' @param a A vector containing the coefficients of the first polynomial.
#' @param b A vector containing the coefficients of the second polynomial.
#' @param p A prime number.
#' @return A vector containing the coefficients of the product modulo p.
#' @keywords internal
polyMultiply <- function(a, b, p) {

  result <- rep(
    0,
    length(a) + length(b) - 1
  )

  for (i in seq_along(a)) {

    for (j in seq_along(b)) {

      result[i + j - 1] <-
        result[i + j - 1] +
        a[i] * b[j]
    }
  }

  return(result %% p)
}
#' Polynomial Reduction over GF(p)
#' Reduces a polynomial modulo a specified irreducible polynomial
#' over GF(p).
#' @param poly A vector containing polynomial coefficients.
#' @param modulus A vector containing the coefficients of the
#'   irreducible polynomial.
#' @param p A prime number.
#' @return The reduced polynomial as a coefficient vector.
#' @keywords internal
polyReduce <- function(poly, modulus, p) {

  poly <- poly %% p

  m <- length(modulus) - 1

  while (length(poly) > m) {

    lead <- poly[length(poly)] %% p

    if (lead != 0) {

      shift <-
        length(poly) -
        length(modulus)

      for (i in seq_along(modulus)) {

        position <- i + shift

        poly[position] <-
          (
            poly[position] -
              lead * modulus[i]
          ) %% p
      }
    }

    poly <- poly[-length(poly)]
  }

  if (length(poly) < m) {

    poly <- c(
      poly,
      rep(0, m - length(poly))
    )
  }

  return(poly %% p)
}
#' Finite-Field Multiplication
#' Multiplies two elements of GF(p^m) using polynomial
#' multiplication and reduction.
#' @param a A finite-field element.
#' @param b A finite-field element.
#' @param p A prime number.
#' @param modulus The irreducible polynomial used to define the field.
#' @return The product of the two finite-field elements.
#' @keywords internal
gfMultiply <- function(a, b, p, modulus) {

  product <- polyMultiply(
    a,
    b,
    p
  )

  result <- polyReduce(
    product,
    modulus,
    p
  )

  return(result)
}
#' Remove Duplicate Finite-Field Elements
#' Removes duplicate rows from a matrix of finite-field elements.
#' @param M A matrix containing finite-field elements as rows.
#' @return A matrix containing the unique finite-field elements.
#' @keywords internal
uniqueGFElements <- function(M) {

  keep <- !duplicated(
    as.data.frame(M)
  )

  return(
    M[
      keep,
      ,
      drop = FALSE
    ]
  )
}
#' Quadratic Residues in a Finite Field
#' Finds the nonzero quadratic residues in GF(p^m).
#' @param p A prime number.
#' @param m A positive integer specifying the field extension degree.
#' @return A matrix containing the nonzero quadratic residues
#'   of GF(p^m).
#' @keywords internal
quadraticResiduesGF <- function(p, m) {

  q <- p^m

  GF <- generateGF(p, m)

  modulus <- irreduciblePolynomial(
    p,
    m
  )

  squares <- matrix(
    0,
    nrow = q - 1,
    ncol = m
  )

  for (i in 2:q) {

    element <- GF[i, ]

    squares[i - 1, ] <- gfMultiply(
      element,
      element,
      p,
      modulus
    )
  }

  QR <- uniqueGFElements(squares)

  zeroRows <- apply(
    QR,
    1,
    function(x) all(x == 0)
  )

  QR <- QR[
    !zeroRows,
    ,
    drop = FALSE
  ]

  return(QR)
}
#' Quadratic Character in a Finite Field
#' Computes the quadratic character of an element of GF(p^m).
#' It returns 0 for the zero element, 1 for a nonzero quadratic
#' residue, and -1 otherwise.
#' @param x A finite-field element.
#' @param QR A matrix containing the nonzero quadratic residues
#'   of the finite field.
#' @return The quadratic character: 0, 1, or -1.
#' @keywords internal
quadraticCharacterGF <- function(x, QR) {

  # chi(0) = 0
  if (all(x == 0))
    return(0)

  # chi(x) = 1 if x is a nonzero square
  for (i in seq_len(nrow(QR))) {

    if (all(x == QR[i, ]))
      return(1)
  }

  # chi(x) = -1 otherwise
  return(-1)
}
