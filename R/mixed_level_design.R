#' Constructs an orthogonal design for two types of factors with different
#' numbers of levels. Three construction cases are considered: even-even level factors,
#' odd-even level factors and odd-odd level factors. The appropriate construction is selected
#' automatically according to the parity of the numbers of levels of the
#' two factors. Hadamard matrices are used for even-level factors, while
#' conference matrices are used for odd-level factors. The order of both
#' the Hadamard and conference matrices is `m`. For the even-even case,
#' the design contains `m - 1` factors of each type. For the odd-even case,
#' the design contains `m` odd-level factors and `m - 1` even-level factors.
#' For the odd-odd case, the design contains `m` factors of each type.
#' @param a Number of levels of the first factor.
#' @param b Number of levels of the second factor.
#' @param m Order of the Hadamard and conference matrices used in the
#'   construction. Must be a multiple of 4.
#' @return A matrix containing the constructed orthogonal design. The total
#'   number of factors depends on the construction case: `2(m - 1)` for
#'   the even-even case, `2m - 1` for the odd-even case, and `2m` for
#'   the odd-odd case.
#' @examples
#' mixed_level_design(a = 4, b = 2, m = 4)
#' mixed_level_design(a = 3, b = 2, m = 8)
#' mixed_level_design(a = 5, b = 3, m = 4)
#' @export
mixed_level_design <- function(a, b, m) {

  if (a %% 2 == 0 && b %% 2 == 0) {

    design <- even_even_design(a, b, m)

  } else if (a %% 2 != 0 && b %% 2 != 0) {

    design <- odd_odd_design(a, b, m)

  } else {

    design <- odd_even_design(a, b, m)
  }

  return(design)
}
