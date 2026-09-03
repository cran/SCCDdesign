#' Constructs a pseudo-conference matrix of order \code{n}
#' using a determinant-based search procedure.
#' @param n A positive integer specifying the order of the matrix.
#' @param starts Number of random starting matrices used in the
#' search. Default is 1000.
#' @param max.iter Maximum number of iterations for each starting
#' matrix. Default is 100.
#' @return A list containing the pseudo-conference matrix and
#' the corresponding log-determinant value.
#' @details
#' The matrix has zero diagonal elements and off-diagonal
#' elements equal to -1 or +1. The search selects a matrix
#' by maximizing the determinant criterion.
#' @examples
#' result <- pseudoConference(6, starts = 10, max.iter = 10)
#' result$Matrix
#' result$LogDet
#' @export
pseudoConference <- function(n,
                             starts = 1000,
                             max.iter = 100) {

  if (missing(n))
    stop("Please specify the order n.")

  if (!is.numeric(n) ||
      length(n) != 1 ||
      is.na(n) ||
      !is.finite(n)) {

    stop("n must be a single finite numeric value.")
  }

  if (n != as.integer(n))
    stop("n must be an integer.")

  n <- as.integer(n)

  if (n < 2)
    stop("n must be at least 2.")


  if (!is.numeric(starts) ||
      length(starts) != 1 ||
      is.na(starts) ||
      starts < 1 ||
      starts != as.integer(starts)) {

    stop("starts must be a positive integer.")
  }

  starts <- as.integer(starts)


  if (!is.numeric(max.iter) ||
      length(max.iter) != 1 ||
      is.na(max.iter) ||
      max.iter < 1 ||
      max.iter != as.integer(max.iter)) {

    stop("max.iter must be a positive integer.")
  }

  max.iter <- as.integer(max.iter)

  best.det <- -Inf

  best.C <- NULL

  for (start in seq_len(starts)) {

    C <- matrix(
      sample(
        c(-1, 1),
        n * n,
        replace = TRUE
      ),
      nrow = n,
      ncol = n
    )

    diag(C) <- 0

    improved <- TRUE

    iter <- 1


    while (improved && iter <= max.iter) {

      improved <- FALSE

      for (i in seq_len(n)) {

        for (j in seq_len(n)) {

          if (i != j) {

            old <- C[i, j]


            det.old <- determinant(
              crossprod(C),
              logarithm = TRUE
            )$modulus

            det.old <- as.numeric(det.old)


            C[i, j] <- 1

            det.plus <- determinant(
              crossprod(C),
              logarithm = TRUE
            )$modulus

            det.plus <- as.numeric(det.plus)

            C[i, j] <- -1

            det.minus <- determinant(
              crossprod(C),
              logarithm = TRUE
            )$modulus

            det.minus <- as.numeric(det.minus)

            if (
              det.plus >= det.minus &&
              det.plus > det.old
            ) {

              C[i, j] <- 1

              improved <- TRUE

            } else if (
              det.minus > det.old
            ) {

              C[i, j] <- -1

              improved <- TRUE

            } else {

              C[i, j] <- old
            }
          }
        }
      }

      iter <- iter + 1
    }

    det.final <- determinant(
      crossprod(C),
      logarithm = TRUE
    )$modulus

    det.final <- as.numeric(det.final)

    if (det.final > best.det) {

      best.det <- det.final

      best.C <- C
    }
  }

  if (is.null(best.C))
    stop("Pseudo-conference matrix construction failed.")


  if (!all(diag(best.C) == 0))
    stop("Diagonal elements of the matrix are not zero.")


  if (!all(best.C %in% c(-1, 0, 1)))
    stop("Pseudo-conference matrix contains invalid entries.")


  colnames(best.C) <- paste0(
    "X",
    seq_len(n)
  )

  return(
    list(
      Matrix = best.C,
      LogDet = best.det
    )
  )
}

