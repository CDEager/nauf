
#' Get the model frame for a model with unordered factor \code{NA} values.
#'
#' \code{nauf_model_frame} calls the default \code{\link[stats]{model.frame}}
#' method with \code{na.action = "na.pass"} and then sets sum contrasts for
#' unordered factors and checks for \code{NAs} that cause factor levels to
#' need to be dropped in interactions.
#'
#' After calling the default \code{\link[stats]{model.frame}} method,
#' the function converts any character and logical columns to unordered factors,
#' and any other column which has only two unique non-\code{NA} values,
#' regardless of their class, are also coerced to unordered factors (e.g.
#' a factor coded as ordered even though it only has two levels, or a binary
#' variable coded as integer \code{0 / 1}).  Then all unordered factors are
#' assigned named sum contrasts via \code{\link{named_contr_sum}}.
#'
#' After this, all columns are checked for \code{NA} values.  If any columns
#' have \code{NAs} but are not unordered factors, then an error is thrown
#' (these types of \code{NA} values are not the same as the ones dealt with in
#' this package; there are many options for imputing them, but this needs to be
#' done prior to calling any \code{nauf} functions).  If there are any unordered
#' factors with \code{NA} values which are involved in interactions with other
#' unordered factors, then for each of these interactions
#' \code{\link{nauf_interaction}} is called, and, if changes in contrasts are
#' required in the interaction term, they are recorded.
#'
#' @param formula A \code{\link[stats]{formula}}. \code{\link[stats]{terms}}
#'   attributes are ignored.
#' @param data A data.frame containing the variables for the model.
#' @param na.action Changes from the default value (\code{\link[stats]{na.pass}})
#'   are ignored (\code{NA} values in factors are not removed).
#' @param drop.unused.levels Changes from the default value (\code{FALSE}) are
#'   ignored (factor levels with no observations are dropped).
#' @param xlev Changes from the default value (\code{NULL}) are ignored. Levels
#'   for unordered factors are determined alphabetically so they are consistent
#'   when reordered in interactions (see \code{\link{nauf_interaction}}).
#' @param ... Other parameters to pass to the default
#'   \code{\link[stats]{model.frame}} method such as \code{subset},
#'   \code{weights} and \code{offset}.
#'
#' @return A model frame that inherits from \code{\link[=nauf-class]{nauf}} and
#'   whose \code{\link[stats]{terms}} attribute inherits from
#'   \code{\link[=nauf-class]{nauf}} and has three additional attributes:
#'   \code{mefc} (main effect factor contrasts) is a named list with an element
#'   for each factor in the model, containing its contrast matrix, levels,
#'   and a logical indicating whether the factor is ordered or not; \code{ccna}
#'   (contrasts changed due to \code{NA} values) is a named list with an entry
#'   for each unordered factor interaction which requires contrasts to be
#'   changed from those applied to the main effects, containing a named list
#'   of levels which are kept in the interaction (see
#'   \code{\link{nauf_interaction}}); and \code{hasna} is a named logical vector
#'   indicating which variables in the model frame have \code{NA} values.
#'
#' @examples
#' \dontrun{
#' fr <- nauf_model_frame(y ~ f1 * f2 * (f3 + x1), mydata)
#' }
#'
#' @export
nauf_model_frame <- function(formula, data = NULL, na.action = "na.pass",
                     drop.unused.levels = TRUE, xlev = NULL, ...) {

  if (na.action != "na.pass") {
    warning("'na.action' must be na.pass and changes to ",
      "the default are ignored.")
  }
  if (!drop.unused.levels) {
    warning("'drop.unused.levels' must be TRUE and changes to ",
      "the default are ignored.")
  }
  if (!is.null(xlev)) {
    warning("'xlev' must be NULL and changes to ",
      "the default are ignored.")
  }

  mf <- stats::model.frame(nauf_on(formula), data = data, na.action = "na.pass",
    drop.unused.levels = TRUE, xlev = NULL, ...)

  return(mf)
}


#' @describeIn nauf_model_frame S3 method for class 'nauf'
#'
#' @export
model.frame.nauf <- function(formula, data = NULL, na.action = "na.pass",
                    drop.unused.levels = TRUE, xlev = NULL, ...) {

  if (inherits(formula, "lm") && is.nauf(formula)) {
    return(formula$model)
  }
  data <- nauf_off(data)
  if (!inherits(formula, "formula")) {
    stop("must supply a formula or terms object")
  }
  if (!is.data.frame(data)) stop("must supply a data.frame")
  if ("ccna" %in% names(attributes(formula))) {
    formula <- nauf_off(formula)
    mf <- stats::model.frame(formula, data, na.action = "na.pass",
      drop.unused.levels = FALSE, xlev = NULL, ...)
    mefc <- attr(formula, "mefc")
    for (j in names(mefc)) {
      mf[, j] <- factor(mf[, j], ordered = mefc[[j]]$ordered,
        levels = mefc[[j]]$levels)
      contrasts(mf[, j]) <- mefc[[j]]$contrasts
    }
    attr(mf, "terms") <- nauf_on(formula)

  } else {
    attr(formula, "terms") <- NULL
    class(formula) <- "formula"
    mf <- stats::model.frame(formula, data, na.action = "na.pass",
      drop.unused.levels = FALSE, xlev = NULL, ...)
    if (stats::is.empty.model(mf)) {
      stop("model is empty")
    }

    mt <- attr(mf, "terms")
    mta <- attributes(mt)
    fmat <- mta$factors
    # gets rid of both response, offset, etc
    v <- rownames(fmat)[rowSums(fmat) > 0]
    nval <- unlist(lapply(mf[, v, drop = FALSE],
      function(n) length(sort(unique(n)))))
    dc <- mta$dataClasses[v]
    dc[dc == "logical" | nval == 2] <- "factor"
    uf <- dc == "factor"
    of <- dc == "ordered" & !uf
    hasna <- unlist(lapply(mf, anyNA)[v])
    if (any(hasna & !uf)) {
      stop("the following have NA values but are not unordered factors: ",
        v[hasna & !uf])
    }
    nauf <- hasna & uf
    if (any(fmat > 1)) {
      warning("nauf has not been tested for models that violate ",
        "the interaction hierarchy.")
    }
    if (!mta$intercept) {
      warning("nauf has not been tested for zero-intercept models.")
    }

    mefc <- list()
    for (j in v[uf]) {
      mf[, j] <- named_contr_sum(mf[, j], FALSE)
      mefc[[j]] <- list(ordered = FALSE, levels = levels(mf[, j]),
        contrasts = contrasts(mf[, j]))
      attr(attr(mf, "terms"), "dataClasses")[j] <- "factor"
    }
    ofc <- as.character(getOption("contrasts")[2])
    for (j in v[of]) {
      # if you don't do this then the attributes aren't equal
      # in testing for the ordered factors
      if (is.null(attr(mf[, j], "contrasts"))) {
        contrasts(mf[, j]) <- ofc
      } else if (any(!(xtabs(~ mf[, j])))) {
        mf[, j] <- droplevels(mf[, j])
        contrasts(mf[, j]) <- ofc
        warning("Dropping levels from ordered factor ", j)
      }
      attr(mf[, j], "contrasts") <- contrasts(mf[, j])
      rownames(attr(mf[, j], "contrasts")) <- levels(mf[, j])
      mefc[[j]] <- list(ordered = TRUE, levels = levels(mf[, j]),
        contrasts = contrasts(mf[, j]))
    }

    ccna <- list()
    if (any(nauf)) {
      vmat <- fmat[v, , drop = FALSE]
      ufmat <- vmat[v[uf], , drop = FALSE]
      naufmat <- ufmat[v[nauf], , drop = FALSE]
      has2uf <- colSums(ufmat) > 1
      hasnauf <- colSums(naufmat) > 0
      check_naufi <- which(has2uf & hasnauf)

      if (length(check_naufi)) {
        check_ccna <- list()
        for (j in 1:length(check_naufi)) {
          # loop since using apply() could return matrix instead of list
          check_ccna[[j]] <- sort(rownames(ufmat)[ufmat[, check_naufi[j]] > 0])
        }
        check_ccna <- unique(check_ccna)

        for (i in check_ccna) {
          checki <- nauf_interaction(mf, i)
          if (checki$changed) {
            nn <- paste(i, collapse = ":")
            ccna[[nn]] <- list(factors = list(),
              assign = which((colSums(ufmat) == length(i)) & apply(
                ufmat[i, , drop = FALSE], 2, all)))
            for (j in i) {
              ccna[[nn]]$factors[[j]] <- list(ordered = FALSE,
                levels = checki$levels[[j]],
                contrasts = named_contr_sum(checki$levels[[j]]))
            }
          }
        }
      }
    }

    attr(attr(mf, "terms"), "hasna") <- hasna
    attr(attr(mf, "terms"), "mefc") <- mefc
    attr(attr(mf, "terms"), "ccna") <- ccna
    attr(mf, "terms") <- nauf_on(attr(mf, "terms"))
  }

  return(nauf_on(mf))
}

