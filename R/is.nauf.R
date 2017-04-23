

#' Determine if an object inherits from a \code{nauf} class.
#'
#' The \code{is.nauf} functions return \code{TRUE} if \code{object} inherits
#' from the class (or one of the classes) implied in the name of the funciton
#' and \code{FALSE} otherwise.
#'
#' The following inheritances will result in a value of \code{TRUE}:
#' \describe{
#'   \item{is.nauf.formula}{\code{\linkS4class{nauf.formula}}.}
#'   \item{is.nauf.terms}{\code{\linkS4class{nauf.terms}}.}
#'   \item{is.nauf.frame}{\code{\linkS4class{nauf.frame}}.}
#'   \item{is.nauf.glm}{\code{\linkS4class{nauf.glm}}.}
#'   \item{is.nauf.lmerMod}{\code{\linkS4class{nauf.lmerMod}}.}
#'   \item{is.nauf.glmerMod}{\code{\linkS4class{nauf.glmerMod}}.}
#'   \item{is.nauf.merMod}{\code{\linkS4class{nauf.lmerMod}} or
#'     \code{\linkS4class{nauf.glmerMod}}.}
#'   \item{is.nauf.model}{\code{\linkS4class{nauf.glm}},
#'     \code{\linkS4class{nauf.lmerMod}}, or
#'     \code{\linkS4class{nauf.glmerMod}}.}
#'   \item{is.nauf.ref.grid}{\linkS4class{nauf.ref.grid}.}
#' }
#'
#' @param object Any R object.
#'
#' @name is.nauf
NULL


#' @rdname is.nauf
#' @export
is.nauf.formula <- function(object) {
  return(inherits(object, "nauf.formula"))
}


#' @rdname is.nauf
#' @export
is.nauf.terms <- function(object) {
  return(inherits(object, "nauf.terms"))
}


#' @rdname is.nauf
#' @export
is.nauf.frame <- function(object) {
  return(inherits(object, "nauf.frame"))
}


#' @rdname is.nauf
#' @export
is.nauf.glm <- function(object) {
  return(inherits(object, "nauf.glm"))
}


#' @rdname is.nauf
#' @export
is.nauf.lmerMod <- function(object) {
  return(inherits(object, "nauf.lmerMod"))
}


#' @rdname is.nauf
#' @export
is.nauf.glmerMod <- function(object) {
  return(inherits(object, "nauf.glmerMod"))
}


#' @rdname is.nauf
#' @export
is.nauf.merMod <- function(object) {
  return(inherits(object, c("nauf.lmerMod", "nauf.glmerMod")))
}


#' @rdname is.nauf
#' @export
is.nauf.model <- function(object) {
  return(inherits(object, c("nauf.glm", "nauf.lmerMod", "nauf.glmerMod")))
}


#' @rdname is.nauf
#' @export
is.nauf.ref.grid <- function(object) {
  return(inherits(object, "nauf.ref.grid"))
}

