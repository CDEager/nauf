
#' Run a regression on data with \code{NA} values in the unordered factors
#'
#' \code{nauf_reg} calls either \code{\link[stats]{lm}},
#' \code{\link[stats]{glm}}, or \code{\link[MASS]{glm.nb}} depending on the
#' regression family, ensuring that the \code{nauf} methods are used to form
#' the model frame and model matrix, and also ensuring that generic functions
#' used on the fitted model object (predict, summary, etc.) use these methods.
#'
#' The \code{formula}, \code{family}, and \code{data} arguments are the same
#' as would be passed to \code{\link[stats]{lm}} and \code{\link[stats]{glm}},
#' except that \code{family = "negbin"} is allowed and fits a negative binomial
#' regression.  The regression fitting function that is called depends on the
#' value of \code{family}.  If \code{family = gaussian(link = identity)}, which
#' is the default, then \code{\link[stats]{lm}} is called.
#' If \code{family = "negbin"}, then \code{\link[MASS]{glm.nb}} is called.
#' Any other value for \code{family} leads to a call to \code{\link[stats]{glm}}.
#' See the help pages for these functions to see what additional arguments can
#' be passed to them through \code{...}, as they differ between the functions.
#'
#' The only other difference between a call to \code{nauf_reg} and a call to
#' the model-fitting functions that it calls are the following arguments whose
#' values cannot be changed (you can enter different values for them, but they
#' will be ignored): \code{model} and \code{x} are always \code{TRUE}, as
#' some generic functions will not work properly otherwise;  \code{na.action}
#' is always set to \code{\link[stats]{na.pass}}, since otherwise the
#' observations with \code{NA} values would be dropped rather than coded
#' using \code{nauf} methods; and \code{contrasts} is always set to \code{NULL}
#' (contrasts for unordered factors are forced to \code{\link{named_contr_sum}}
#' and contrasts for ordered factors with 3 or more levels are not altered; to
#' set custom contrasts for ordered factors, do so in the data frame prior
#' to calling \code{nauf_reg}).
#'
#' @param formula A \code{\link[stats]{formula}} describing the model to be fit.
#' @param family A regression \code{\link[stats]{family}} or, for negative
#'   binomial regressions, a character string \code{"negbin"}.  The default
#'   is \code{gaussian}.
#' @param data A \code{data.frame} with the variables in \code{formula}.
#' @param ... Additional parameters passed to the fitting function (depends on
#'   \code{family}.  See 'Details').
#'
#' @return A fitted model object which inherits from
#'   \code{\link[=nauf-class]{nauf}} and \code{lm}, and possibly also from
#'   \code{mlm}, \code{glm}, or \code{negbin} depending on \code{family}.
#'
#' @seealso For information about the elements in the fitted model object and
#'   about the generic functions that can be applied (predict, summary, etc.),
#'   see \code{\link[stats]{lm}} for \code{family = gaussian},
#'   \code{\link[MASS]{glm.nb}} for \code{family = "negbin"}, and
#'   \code{\link[stats]{glm}} for all other values of \code{family}.  For
#'   information about the handling of \code{NA} values, see
#'   \code{\link{named_contr_sum}} and \code{\link{nauf_interaction}} for
#'   unordered factor contrasts and \code{\link{nauf_model_frame}} and
#'   \code{\link{nauf_model_matrix}} for how these contrasts are implemented.
#'
#' @examples
#' \dontrun{
#' # linear model (calls lm)
#' m <- nauf_reg(y ~ f1 * f2 * (f3 + x), mydata)
#' # logistic model (calls glm)
#' m <- nauf_reg(y ~ f1 * f2 * (f3 + x), mydata, family = binomial)
#' # poisson model (calls glm)
#' m <- nauf_reg(y ~ f1 * f2 * (f3 + x), mydata, family = poisson)
#' # negative binomial model (calls glm.nb)
#' m <- nauf_reg(y ~ f1 * f2 * (f3 + x), mydata, family = "negbin")
#' }
#'
#' @export
nauf_reg <- function(formula, data, family = gaussian, ...) {
  mc <- match.call()
  mc$model <- TRUE
  mc$x <- TRUE
  mc$contrasts <- NULL
  mc$na.action <- "na.pass"

  mce <- mc
  mce$formula <- call("nauf_on", formula)

  if (is.function(family)) {
    family <- family()
    if (is.null(family$family)) {
      stop("'family' not recognized")
    }
    if (family$family == "gaussian" && family$link == "identity") {
      mce[[1]] <- quote(stats::lm)
      mce <- mce[-which(names(mce) == "family")]
    } else {
      mce[[1]] <- quote(stats::glm)
    }
  } else if (is.character(family)) {
    if (family == "gaussian") {
      mce[[1]] <- quote(stats::lm)
    } else if (family == "negbin") {
      mce[[1]] <- quote(MASS::glm.nb)
      mce$family <- NULL
    } else {
      stop("'family' not recognized")
    }
  } else {
    stop("'family' not recognized")
  }

  mod <- eval(mce, parent.frame())
  mod$call <- mc

  return(nauf_on(mod))
}

