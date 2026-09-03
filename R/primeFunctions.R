#' Check for a Prime Number.
#' Checks whether a given number is prime.
#' @param q A numeric value to be checked.
#' @return TRUE if q is prime and FALSE otherwise.
#' @keywords internal
is_prime <- function(q) {

  if (length(q) != 1 || is.na(q))
    return(FALSE)

  if (q <= 1)
    return(FALSE)

  if (q == 2)
    return(TRUE)

  if (q %% 2 == 0)
    return(FALSE)

  limit <- floor(sqrt(q))

  if (limit >= 3) {

    for (i in seq(3, limit, by = 2)) {

      if (q %% i == 0)
        return(FALSE)
    }
  }

  return(TRUE)
}
#' Identify a Prime Power
#' Determines whether q can be expressed as p^m, where p
#' is a prime number and m is a positive integer.
#' @param q A positive integer to be checked.
#' @return A list containing p and m if q is a prime or
#' prime power; otherwise NULL.
#' @keywords internal
primePowerInfo <- function(q) {

  if (q <= 1)
    return(NULL)

  # Prime case: q = q^1
  if (is_prime(q)) {

    return(
      list(
        p = q,
        m = 1
      )
    )
  }
max_m <- floor(log(q, base = 2))

  for (m in 2:max_m) {

    p0 <- round(q^(1 / m))
    candidates <- unique(
      c(p0 - 1, p0, p0 + 1)
    )

    candidates <- candidates[candidates >= 2]

    for (p in candidates) {

      if (is_prime(p) && p^m == q) {

        return(
          list(
            p = p,
            m = m
          )
        )
      }
    }
  }

  return(NULL)
}
#' Check for a Prime Power.
#' Checks whether q is a prime or a power of a prime number.
#' @param q A positive integer to be checked.
#' @return TRUE if q is a prime or prime power and FALSE otherwise.
#' @keywords internal
is_prime_power <- function(q) {

  !is.null(primePowerInfo(q))
}
