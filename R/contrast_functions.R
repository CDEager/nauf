
#' Create named sum contrasts for a factor.
#'
#' \code{named_contr_sum} creates sum contrasts for a factor which are named
#' with the levels of the factor rather than with number (e.g. if a factor
#' \code{f1} has levels \code{A}, \code{B}, and \code{C}, then rather than
#' creating contrast columns \code{f11} and \code{f12}, it creates columns
#' \code{f1A} and \code{f1B}).  The absolute value of the non-zero elements
#' of the matrix can also be specified.
#'
#' First, \code{x} is coerced to character, and its unique non-\code{NA} values
#' are sorted alphabetically.  If there are two unique values, and they are
#' equal to (ignoring case) "F" and "T", "FALSE" and "TRUE", "N" and "Y",
#' "NO" and "YES", or "0" and "1", then their order is reversed (this makes it
#' so the positive level gets the dummy coefficient rather than the negative
#' level, which has a more intuitive interpretation).  Then
#' \code{\link[stats]{contr.sum}} is called, and the column names of the
#' resulting contrast matrix are set using the character vector of unique values
#' (excluding the final element that gets coded as \code{-1} for all dummy
#' variables).  This entire matrix is then multiplied by \code{sumcoef}; with
#' the default value of \code{1}, this does not change the matrix; if, for
#' example, \code{sumcoef = 0.5}, then rather than each column containing
#' \code{-1, 0, 1}, each column would contain \code{-0.5, 0, 0.5}.
#' If \code{return_contr = TRUE}, then this contrast matrix is
#' returned.  If \code{return_contr = FALSE}, then \code{x} is converted to an
#' unordered factor with the named sum contrats and returned. \code{NA} is never
#' assigned as a level in the contrast matrix or in the factor returned by the
#' function, but \code{NA} values in \code{x} are not removed in the factor
#' returned when \code{return_contr = FALSE}. See the examples.
#'
#' @param x A factor or vector, with at least two unique non-\code{NA} values,
#'   which can be coerced to character with \code{\link[base]{as.character}}.
#' @param sumcoef A positive number indicating by which the entire contrast
#'   matrix returned by \code{\link[stats]{contr.sum}} is multiplied.  See
#'   'Details'.
#' @param return_contr A logical. If \code{TRUE} (the default), a contrast
#'   matrix is returned. If \code{FALSE}, \code{x} is converted to an unordered
#'   factor with the contrast matrix applied, and the factor is returned.
#'
#' @return If \code{return_contr = TRUE}, a contrast matrix obtained from
#'   \code{\link[stats]{contr.sum}} with named columns rather than numbered
#'   columns.  If \code{return_contr = FALSE}, then \code{x} is returned
#'   as an unordered factor with the named sum contrasts applied.
#'
#' @seealso \code{\link{nauf_interaction}} for interactions involving unordered
#'   factors where at least one has \code{NA} values,
#'   \code{\link{nauf_model_frame}} for automatic assignment of unordered factor
#'   contrasts, and \code{\link{nauf_model_matrix}} for the treatment of
#'   \code{NAs} in regressions.
#'
#' @examples
#' f <- factor(rep(c("a", "b", "c", NA), 2), levels = c("b", "c", "a"))
#' f <- addNA(f)
#' levels(f)  # NA listed as factor
#' contrasts(f)  # NA included in contrast matrix
#' named_contr_sum(f)  # named sum contrasts (NA dropped; levels alphabetized)
#' named_contr_sum(levels(f))  # same output
#' named_contr_sum(f, FALSE)  # f (values unchanged) with named sum contrasts
#'
#' f <- c(TRUE, FALSE, FALSE, TRUE)
#' class(f)  # logical
#' named_contr_sum(f)  # TRUE gets the dummy variable
#' f <- named_contr_sum(f, FALSE)
#' class(f)  # factor
#'
#' named_contr_sum(letters[1:5])  # character argument
#' named_contr_sum(rep(letters[1:5], 2), FALSE)  # creates factor
#'
#' named_contr_sum(c(1, 6, 3))  # numeric argument
#' named_contr_sum(rep(c(1, 6, 3), 2), FALSE)  # creates factor
#'
#' # ordered factors are converted to unordered factors, so use with caution
#' f <- factor(rep(1:3, 2), ordered = TRUE)
#' is.ordered(f)
#' f
#' f <- named_contr_sum(f, FALSE)
#' is.ordered(f)
#' f
#'
#' \dontrun{
#' # error from stats::contr.sum because only one unique non-NA value
#' named_contr_sum(5)
#' named_contr_sum(rep(c("a", NA), 3), FALSE)
#' }
#'
#' @export
named_contr_sum <- function(x, sumcoef = 1, return_contr = TRUE) {
  if (length(dim(x))) stop("'x' must be a vector or factor")
  if (!is.scalar(sumcoef, 1)) {
    stop("'sumcoef' must be a single positive number")
  }
  
  x <- as.character(x)
  n <- length(lvs <- reorder_ft(sort(unique(x))))
  contr <- stats::contr.sum(lvs) * sumcoef
  colnames(contr) <- lvs[-n]

  if (return_contr) return(contr)

  x <- factor(x, ordered = FALSE, levels = lvs)
  contrasts(x) <- contr
  return(x)
}


#' Unordered factor interactions involving \code{NAs}.
#'
#' The function determines whether unordered factor levels need to be dropped
#' in the computation of an interaction term due to \code{NA} values, treating
#' characters and logicals as unordered factors.
#'
#' If there are not at least two unordered factors in the specified columns,
#' an error is thrown. Otherwise, the unordered factors involved in the
#' interaction are subsetted.  Then a table of all their unique combinations in
#' \code{x} is made, excluding any combination which contains \code{NAs}.
#' Then levels which are not used in the no-\code{NA} subset are dropped.
#' Then, if any level in a factor is collinear with a level in another factor,
#' the level is dropped.  If after doing this, there are still empty cells
#' in the full contingency table of the factors, a warning is issued
#' indicating that there is collinearity in the interaction term (the term
#' shouldn't be included and will cause \code{NA} regression coefficients).
#' The levels for each factor which are applicable in the interaction are
#' returned along with a logical indicating if any levels were dropped. Note
#' that only the highest order interaction of the unordered factors is
#' considered.  If unordered factors \code{A}, \code{B}, and \code{C} are
#' involved in an interaction, it may be that \code{A} needs levels dropped in
#' the third-order interaction, but not in the interaction \code{A:B}, for
#' example.  The function would need to be called separately with
#' \code{cols = c("A", "B")} to determine the levels for the second-order
#' interaction.
#'
#' This function is implemented in \code{\link{nauf_model_frame}}, and is
#' necessary in cases where \code{NA} values are collinear with non-\code{NA}
#' values in other unordered factors.  For example, if factor \code{f1}
#' has levels \code{A}, \code{B}, and \code{C}, and factor \code{f2} has
#' levels \code{D} and \code{E} when \code{f1} is \code{A} or \code{B}, but
#' \code{f2} is always \code{NA} when \code{f1} is \code{C}, then the main
#' effects for these factors are automatically assigned the following contrasts.
#'
#' Main effect contrasts for \code{f1}.
#' \tabular{lrr}{
#' f1 \tab  f1A \tab  f1B \cr
#' A  \tab  1   \tab    0 \cr
#' B  \tab  0   \tab    1 \cr
#' C  \tab -1   \tab   -1
#' }
#'
#' Main effect contrasts for \code{f2}.
#' \tabular{lr}{
#' f2 \tab  f2D \cr
#' D  \tab    1 \cr
#' E  \tab   -1 \cr
#' NA \tab    0
#' }
#'
#' This setup allows the regression coefficient \code{f2D} to
#' only apply when \code{f1} is either \code{A} or \code{B}, provided that
#' whenever \code{f1} is \code{C}, \code{f2} is \code{NA}.
#' Provided that all covariates have been put on unit scale (see
#' \code{\link[base]{scale}}) and orthogonal polynomial contrasts have been set
#' for all ordered factors (see \code{\link[stats]{contr.poly}}), both of which
#' are good ideas for most regressions anyway, then these types of contrasts
#' also result in the intercept being the corrected mean, and the output is
#' easier to interpret.  If we are interested in the
#' interaction \code{f1 * f2}, rather than just the main effects, then merely
#' multiplying the main effect contrasts to obtain interaction contrasts
#' yields the following undesirable contrasts.
#'
#' Undesirable contrasts for \code{f1:f2} interaction.
#' \tabular{llrr}{
#' f1 \tab f2 \tab f1A:f2D \tab f1B:f2D \cr
#' A  \tab D  \tab       1 \tab       0 \cr
#' A  \tab E  \tab      -1 \tab       0 \cr
#' B  \tab D  \tab       0 \tab       1 \cr
#' B  \tab E  \tab       0 \tab      -1 \cr
#' C  \tab NA \tab       0 \tab       0
#' }
#'
#' These contrasts are undesirable because the same interaction term could be
#' expressed with one column with some releveling, and in cases more complicated
#' than this one, there could also be collinearity.  So, instead, we should use
#' the following contrasts:
#'
#' Contrasts for \code{f1:f2} as implemented in \code{nauf}.
#' \tabular{llrr}{
#' f1 \tab f2 \tab f1A:f2D \cr
#' A  \tab D  \tab       1 \cr
#' A  \tab E  \tab      -1 \cr
#' B  \tab D  \tab      -1 \cr
#' B  \tab E  \tab       1 \cr
#' C  \tab NA \tab       0
#' }
#'
#' For this example \code{nauf_interaction} determines that \code{C} is
#' not applicable in the interaction term, and so it drops it from the levels
#' of \code{f1}, while keeping \code{f2} the same.  In this case,
#' \code{nauf_interaction} would return a list with an element \code{levels}
#' which has sub-elements \code{f1} and \code{f2}, each of which is
#' a character vector (\code{c("A", "B")} and \code{c("D", "E")},
#' respectively), along with a second element \code{changed = TRUE} to indicate
#' that at least one level was dropped from at least one factor.  See
#' \code{\link{nauf_model_matrix}} for implementation of interaction contrasts
#' which don't match main effect contrasts. See the examples section for code
#' implementing this example.
#'
#' @param x A data.frame.
#' @param cols A vector specifying columns in \code{x} involved in an
#'   interaction.  At least two must be unordered factors, and at least one must
#'   be an unordered factor with \code{NA} values.  Defaults to all columns in
#'   \code{x}.
#'
#' @return A list with two elements.  The first, \code{levels}, is a named list
#'   with one entry for each unordered factor in \code{cols}.  Each entry is a
#'   character vector of the factor levels which should be coded for in the full
#'   interaction term (with sum contrasts).  The second element, \code{changed},
#'   is a logical indicating whether or not \code{levels} is the same as the
#'   levels in the main effects.
#'
#' @seealso \code{\link{named_contr_sum}} for the contrasts that \code{nauf}
#'   assigns to unordered factors,
#'   \code{\link{nauf_model_frame}} for automatic assignment of unordered factor
#'   contrasts, and \code{\link{nauf_model_matrix}} for the treatment of
#'   \code{NAs} in regressions.
#'
#' @examples
#' dat <- data.frame(
#'   f1 = c("A", "A", "B", "B", "C"),
#'   f2 = c("D", "E", "D", "E", NA))
#' nauf_interaction(dat)  # drops C from f1
#'
#' dat <- as.data.frame(lapply(dat, function(n) rep(n, 10)))
#' dat$x <- rnorm(nrow(dat))
#' dat$o <- factor(rep(1:3, length.out = 50), ordered = TRUE)
#' nauf_interaction(dat)  # ignores x and o; same result as above
#'
#' \dontrun{
#' nauf_interaction(dat, "f1")  # error; only one column
#' nauf_interaction(dat, c("f1", "x"))  # error; only one unordered factor
#' }
#'
#' @export
nauf_interaction <- function(x, cols = colnames(x)) {
  if (!is.data.frame(x)) stop("'x' must be a data.frame'")
  if (!is.character(cols)) {
    if (is.numeric(cols) || is.logical(cols)) {
      cols <- colnames(x)[cols]
    } else {
      stop("'cols' must be a character, numeric, or logical vector")
    }
  }
  cols <- unique(cols)
  nm <- paste(cols, collapse = ":")

  # coerce char/logical/binary to unordered factor and record levels
  x <- droplevels(charlogbin_to_uf(x[, cols, drop = FALSE]))
  uf <- sapply(x, is.uf)
  if (sum(uf) < 2) {
    stop("interaction ", nm, " does not involve at least two unordered factors")
  }
  x <- x[, uf]
  mlvs <- lapply(x, function(n) reorder_ft(sort(levels(n))))

  # remove any unique combination which involves NAs
  x <- unique(x)
  x <- droplevels(x[!rowna(x), , drop = FALSE])
  if (!nrow(x)) {
    stop("No unique applicable combinations in ", nm)
  }
  
  # remove levels which are redundant
  # e.g. f1 in [a, b, c], f2 in [d, e, f, g]
  #    when f1 = a, f2 in [d, e, f]
  #    when f1 = b, f2 in    [e, f, g]
  #    when f1 = c, f2 = NA
  #    then at this point [c] has been dropped form f1,
  #    but we still need to drop [d, g] from f2
  if (empty_cells(x)) {
    torm <- vector("list", ncol(x))
    for (j in 1:length(torm)) {
      for (i in levels(x[, j])) {
        check <- droplevels(x[x[, j] == i, -j, drop = FALSE])
        if (any(sapply(check, nlevels) == 1)) {
          torm[[j]] <- c(torm[[j]], i)
        }
      }
    }
    for (j in 1:length(torm)) {
      if (length(torm[[j]])) {
        x[x[, j] %in% torm[[j]], j] <- NA
      }
    }
    x <- droplevels(x[!rowna(x), , drop = FALSE])
    if (!nrow(x)) {
      stop("No unique applicable combinations in ", nm)
    }
  }
  
  ilvs <- lapply(x, function(n) reorder_ft(sort(levels(n))))
  changed <- !isTRUE(all.equal(mlvs, ilvs))
  if (any(unlist(lapply(x, nlevels)) < 2)) {
    stop("At least one factor in ", nm,
      " has only one level when NAs are removed")
  }
  
  # there can still be collinearity, but in this case it isn't structural
  # so warning rather than error (i.e. a column will be dropped from
  # the model matrix if the interaction is included in a regression)
  if (empty_cells(x)) warning("Collinearity in the interaction ", nm)

  return(list(levels = ilvs, changed = changed))
}


#' Contrasts for a model fit with \code{\link{nauf_reg}}.
#'
#' \code{nauf_contrasts} lists all factor contrasts in a model fit with
#' \code{\link{nauf_reg}}, including changes in contrasts for interactions
#' involving \code{NAs}.
#'
#' @param object A model fit using \code{\link{nauf_reg}}, a model frame made
#'   by \code{\link{nauf_model_frame}}, or the \code{terms} object from either
#'   of these.
#'
#' @return A list with two elements: \code{mefc} (main effect factor contrasts)
#'   and \code{ccna} (contrast changes due to \code{NAs}).
#'
#'   \code{mefc} is a list with two elements \code{unordered} and \code{ordered},
#'   each of which is a named list with an entry for each unordered or ordered
#'   factor in the regression which contains the contrasts applied to the main
#'   effect.
#'
#'   \code{ccna} is a named list with an element for each combination of
#'   unordered factors which required a change in factor contrasts from the
#'   main effects, with an entry for each of the factors containing the
#'   contrasts applied to any term involving those unordered factors (and
#'   possibly also ordered factors and covariates), but no other unordered
#'   factors.
#'
#' @examples
#' \dontrun{
#' m <- nauf_reg(y ~ f1 * f2 * (f3 + x1 + x2), data = mydata)
#' nauf_contrasts(m)
#' }
#'
#' @export
nauf_contrasts <- function(object) {
  if (is.nauf.terms(object)) {
    vars <- attr(object, "nauf")$vars
  } else {
    mf <- model.frame(object)
    if (!is.nauf.frame(object)) {
      stop("must supply a nauf.terms, nauf.frame, or nauf model")
    }
    vars <- mf@nauf$vars
  }
  hasna <- vars$hasna
  uf <- mf@nauf$vars$uf
  of <- mf@nauf$vars$of
  
  contr <- list()
  for (v in names(uf)) {
    lvs <- uf[[v]][[1]]
    v.contr <- named_contr_sum(lvs)
    if (length(uf[[v]]) > 1) {
      for (cj in 2:length(uf[[v]])) {
        vcj.contr <- named_contr_sum(uf[[v]][[cj]])
        cnms <- paste(".c", cj, ".", colnames(vcj.contr), sep = "")
        colnames(vcj.contr) <- cnms
        z <- matrix(0, nrow(v.contr), ncol(vcj.contr))
        rownames(z) <- rownames(v.contr)
        colnames(z) <- colnames(vcj.contr)
        z[rownames(vcj.contr), ] <- vcj.contr
        v.contr <- cbind(v.cont, z)
      }
    }
    if (hasna[v]) {
      v.contr <- rbind(v.contr, 0)
      rownames(v.contr) <- c(rownames(v.contr), NA)
    }
    contr[[v]] <- v.contr
  }
  for (v in names(of)) {
    contr[[v]] <- of[[v]]$contrasts
  }
  
  if (!length(contr)) return(NULL)
  return(contr)
}





































