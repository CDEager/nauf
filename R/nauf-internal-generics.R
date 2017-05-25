

###### get number of chains in a stan object
nchain <- function(x) {
  UseMethod("nchain")
}


nchain.default <- function(x) {
  return(0L)
}


nchain.stanfit <- function(x) {
  return(x@sim$chains)
}


nchain.stanreg <- function(x) {
  return(nchain(x$stanfit))
}



###### get number of iterations in a stan object
niter <- function(x) {
  UseMethod("niter")
}


niter.default <- function(x) {
  return(0L)
}


niter.stanfit <- function(x) {
  return(x@sim$n_save[1])
}


niter.stanreg <- function(x) {
  return(niter(x$stanfit))
}



###### get and alter the nauf.info attribute

nauf.info <- function(x) {
  UseMethod("nauf.info")
}


`.nauf.info<-` <- function(x, value) {
  UseMethod(".nauf.info<-")
}


`nauf.info<-` <- function(x, which, value) {
  if (missing(which)) {
    i <- value
  } else {
    i <- nauf.info(x)
    if (identical(which, 0) || identical(which, 0L)) {
      i[names(value)] <- value
    } else if (length(which) == 1) {
      i[[which]] <- value
    } else {
      i[which] <- value
    }
  }
  
  .nauf.info(x) <- i
  
  return(x)
}


nauf.info.default <- function(x) {
  stop("methods for classes c(", add_quotes(class(x), ", "), ") do not exist.")
}


`.nauf.info<-.default` <- function(x, value) {
  stop("methods for classes c(", add_quotes(class(x), ", "), ") do not exist.")
}


nauf.info.nauf.terms <- function(x) {
  return(attr(x, "nauf.info"))
}


`.nauf.info<-.nauf.terms` <- function(x, value) {
  attr(x, "nauf.info") <- value
  return(x)
}


nauf.info.nauf.frame <- function(x) {
  return(attr(attr(x, "terms"), "nauf.info"))
}


`.nauf.info<-.nauf.frame` <- function(x, value) {
  attr(attr(x, "terms"), "nauf.info") <- value
  return(x)
}


nauf.info.list <- function(x) {
  if (is.glmod(x)) {
    return(attr(attr(x$fr, "terms"), "nauf.info"))
  } else {
    stop("'x' is a list but not a lmod/glmod")
  }
}


`.nauf.info<-.list` <- function(x, value) {
  if (is.glmod(x)) {
    attr(attr(x$fr, "terms"), "nauf.info") <- value
  } else {
    stop("'x' is a list but not a lmod/glmod")
  }
  return(x)
}


nauf.info.nauf.glm <- function(x) {
  return(attr(x$terms, "nauf.info"))
}


`.nauf.info<-.nauf.glm` <- function(x, value) {
  attr(x$terms, "nauf.info") <- value
  attr(attr(x$model, "terms"), "nauf.info") <- value
  return(x)
}


nauf.info.nauf.lmerMod <- function(x) {
  return(attr(attr(x@frame, "terms"), "nauf.info"))
}


`.nauf.info<-.nauf.lmerMod` <- function(x, value) {
  attr(attr(x@frame, "terms"), "nauf.info") <- value
  return(x)
}


nauf.info.nauf.glmerMod <- function(x) {
  return(attr(attr(x@frame, "terms"), "nauf.info"))
}


`.nauf.info<-.nauf.glmerMod` <- function(x, value) {
  attr(attr(x@frame, "terms"), "nauf.info") <- value
  return(x)
}


nauf.info.nauf.stanreg <- function(x) {
  if (is.nauf.stanfer(x)) {
    return(attr(x$terms, "nauf.info"))
  }
  return(attr(attr(x$glmod$fr, "terms"), "nauf.info"))
}


`.nauf.info<-.nauf.stanreg` <- function(x, value) {
  if (is.nauf.stanfer(x)) {
    attr(x$terms, "nauf.info") <- value
    attr(attr(x$model, "terms"), "nauf.info") <- value
  } else {
    attr(attr(x$glmod$fr, "terms"), "nauf.info") <- value
  }
  return(x)
}


nauf.info.ref.grid <- function(x) {
  return(attr(x$ref.grid@model.info$terms, "nauf.info"))
}


`.nauf.info<-.nauf.ref.grid` <- function(x, value) {
  attr(x$ref.grid@model.info$terms, "nauf.info") <- value
  return(x)
}

