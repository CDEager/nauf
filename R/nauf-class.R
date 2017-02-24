
#' S3 class for fitting models with non-applicable unordered factor levels.
#'
#' The \code{nauf} class is used to ensure that the proper S3 methods
#' for \code{link[stats]{model.frame}} and \code{\link[stats]{model.matrix}}
#' are called by already-existing model fitting functions.
#'
#' The \code{nauf} package mostly relies on two functions, both of which are
#' S3 generics: \code{model.frame.nauf} and \code{model.matrix.nauf} with
#' corresponding wrapper functions \code{\link{nauf_model_frame}} and
#' \code{\link{nauf_model_matrix}}.
#'
#' The \code{nauf_on} function adds \code{nauf} as the first
#' \code{\link[base]{class}} attribute to \code{object} if \code{object}
#' \code{\link[base]{inherits}} from \code{\link[stats]{formula}},
#' \code{\link[base]{data.frame}}, or \code{\link[stats]{lm}} and is not
#' an S4 object.  If the object does not meet these critera, it is returned
#' unaltered and a warning is issued.
#'
#' The \code{nauf_off} function removes \code{nauf} from the the class
#' of \code{object} if it is there and then returns the object.
#'
#' The \code{is.nauf} function returns \code{TRUE} if \code{object} inherits
#' from \code{nauf} and \code{FALSE} if it does not.
#'
#' @param object Any R object
#'
#' @name nauf-class
NULL


#' @rdname nauf-class
#' @export
is.nauf <- function(object) {
  return(inherits(object, "nauf"))
}


#' @rdname nauf-class
#' @export
nauf_on <- function(object) {
  if (!isS4(object) && inherits(object, c("formula", "lm", "data.frame",
  "list"))) {
    cl <- class(object)
    cl <- unique(c("nauf", cl))
    class(object) <- cl
  } else {
    warning("incompatible object; returned unaltered")
  }
  return(object)
}


#' @rdname nauf-class
#' @export
nauf_off <- function(object) {
  if (is.nauf(object)) {
    cl <- class(object)
    cl <- cl[cl != "nauf"]
    class(object) <- cl
  }
  return(object)
}

