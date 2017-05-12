



nauf_info <- function(object, ...) {
  UseMethod("nauf_info")
}

nauf_info.nauf.terms <- function(object, ...) {
  return(attr(object, "nauf.info"))
}

nauf_info.nauf.frame <- function(object, ...) {
  return(attr(attr(object, "terms"), "nauf.info"))
}

nauf_info.nauf.lmerMod <- function(object, ...) {
  return(attr(attr(object@frame, "terms"), "nauf.info"))
}

nauf_info.nauf.glmerMod <- function(object, ...) {
  return(attr(attr(object@frame, "terms"), "nauf.info"))
}

nauf_info.nauf.glm <- function(object, ...) {
  return(attr(object$terms, "nauf.info"))
}

nauf_info.nauf.ref.grid <- function(object, ...) {
  return(attr(object$ref.grid@model.info$terms, "nauf.info"))
}

nauf_info.default <- function(object, ...) {
  stop("Cannot extract nauf.info from object with class '",
    class(object)[1],"'")
}



