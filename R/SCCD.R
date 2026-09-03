#' Construct Saturated screening designs for continous and categorical factor (SCCD)
#' Using Method I
#' Constructs an SCCD when the number of runs (n=2m) is divisible by 8,
#' using a conference matrix and a Hadamard matrix.
#' The number of continuous factor is m and the number of categorical factors is m-1)
#' @param n A positive integer specifying the number of runs.
#' @return A numeric SCCD matrix.
#' @keywords internal
SCCD_method1 <- function(n) {

  if (n %% 8 != 0) {

    stop(
      "Method I requires the number of runs to be divisible by 8."
    )
  }
m <- n / 2
C <- conferenceMatrix(m)
H <- hadamardMatrix(m)
Hstar <- H[, -1, drop = FALSE]
top <- cbind(
    C,
    Hstar
  )
bottom <- cbind(
    -C,
    Hstar
  )
D <- rbind(
    top,
    bottom
  )
continuousNames <- paste0(
    "X",
    seq_len(m)
  )
categoricalNames <- paste0(
    "Z",
    seq_len(m - 1)
  )
colnames(D) <- c(
    continuousNames,
    categoricalNames
  )
  return(D)
}
#' Construct SCCD Using Method II
#' Constructs an SCCD when the number of runs (n=2m) is divisible by 4
#' but not divisible by 8, using a conference matrix and a
#' D-optimal two-level design.
#' The number of continuous factor is m and number of categorical factors is m-1
#' @param n A positive integer specifying the number of runs.
#' @return A numeric SCCD matrix.
#' @keywords internal
SCCD_method2 <- function(n) {

  if (n %% 4 != 0 || n %% 8 == 0) {

    stop(
      paste0(
        "Method II requires n to be divisible by 4 ",
        "but not divisible by 8."
      )
    )
  }

  m <- n / 2

  C <- conferenceMatrix(m)

  T <- dOptimalTwoLevel(m)

  if (nrow(C) != m || ncol(C) != m) {

    stop(
      "Conference matrix has incorrect dimensions."
    )
  }

  if (nrow(T) != m || ncol(T) != (m - 1)) {

    stop(
      "D-optimal two-level matrix has incorrect dimensions."
    )
  }


  top <- cbind(
    C,
    T
  )

  bottom <- cbind(
    -C,
    T
  )

  D <- rbind(
    top,
    bottom
  )


  continuousNames <- paste0(
    "X",
    seq_len(m)
  )

  categoricalNames <- paste0(
    "Z",
    seq_len(m - 1)
  )

  colnames(D) <- c(
    continuousNames,
    categoricalNames
  )
  return(D)
}
#' Construct SCCD Using Method III
#' Constructs an SCCD when the number of runs (n=2m) is even but not
#' divisible by 4, using a pseudo-conference matrix and a
#' D-optimal two-level design.
#' The number of continuous factor is m and the number of categorical factors is m-1)
#' @param n A positive integer specifying the number of runs.
#' @param starts Number of random starting matrices used in the
#' pseudo-conference matrix search. Default is 1000.
#' @param max.iter Maximum number of iterations for each starting
#' matrix. Default is 100.
#' @return A numeric SCCD matrix.
#' @keywords internal
SCCD_method3 <- function(n,
                         starts = 1000,
                         max.iter = 100) {


  if (n %% 2 != 0 || n %% 4 == 0) {

    stop(
      paste0(
        "Method III requires n to be even ",
        "but not divisible by 4."
      )
    )
  }
  m <- n / 2

  PC <- pseudoConference(
    n = m,
    starts = starts,
    max.iter = max.iter
  )

  Cstar <- PC$Matrix

  T <- dOptimalTwoLevel(m)

  Cstar <- as.matrix(Cstar)
  T <- as.matrix(T)

  storage.mode(Cstar) <- "numeric"
  storage.mode(T) <- "numeric"

  if (
    nrow(Cstar) != m ||
    ncol(Cstar) != m
  ) {

    stop(
      "Pseudo-conference matrix has incorrect dimensions."
    )
  }


  if (
    nrow(T) != m ||
    ncol(T) != (m - 1)
  ) {

    stop(
      "D-optimal two-level matrix has incorrect dimensions."
    )
  }
  top <- cbind(
    Cstar,
    T
  )

  bottom <- cbind(
    -Cstar,
    T
  )

  D <- rbind(
    top,
    bottom
  )

  continuousNames <- paste0(
    "X",
    seq_len(m)
  )

  categoricalNames <- paste0(
    "Z",
    seq_len(m - 1)
  )

  colnames(D) <- c(
    continuousNames,
    categoricalNames
  )
  return(D)
}
#' Constructs an SCCD for a specified number of runs (n=2m).
#' The construction method is selected automatically according
#' to the number of runs.
#' @param n A positive even integer specifying the number of runs.
#' @return A numeric SCCD matrix containing m continuous and m-1
#' categorical factors.
#' @details
#' Method I is used when n is divisible by 8.
#' Method II is used when n is divisible by 4 but not by 8.
#' Method III is used when n is even but not divisible by 4.
#' @examples
#' SCCD(16)
#' @export
SCCD <- function(n) {

  if (missing(n))
    stop("Please specify the number of runs n.")

  if (
    !is.numeric(n) ||
    length(n) != 1 ||
    is.na(n) ||
    !is.finite(n)
  ) {

    stop(
      "n must be a single finite numeric value."
    )
  }


  if (n != as.integer(n))
    stop("n must be an integer.")


  n <- as.integer(n)


  if (n < 8)
    stop("n must be at least 8.")


  if (n %% 2 != 0)
    stop("n must be even.")

  if (n %% 8 == 0) {

    return(
      SCCD_method1(n)
    )
  }

  if (n %% 4 == 0) {

    return(
      SCCD_method2(n)
    )
  }
  return(
    SCCD_method3(n)
  )
}
