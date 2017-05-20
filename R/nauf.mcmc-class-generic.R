

# a is array[iteration, chain, coefficient]
# x is vector[feature] or matrix[obs, feature]
# returns array[iteration, chain, prediction]
`%*a%` <- function(x, a) {
  x <- t(x)
  
  d <- dim(a)
  d[3] <- ncol(x)
  
  ans <- array(dim = d)
  for (chain in 1:d[2]) {
    ans[, chain, ] <- as(a[, chain, ] %*% x, "matrix")
  }
  
  return(nauf.mcmc(ans))
}


#' @export
nauf.mcmc <- function(samples, chains = NULL, nms = NULL) {
  if (is.nauf.mcmc(samples)) return(samples)
  
  if (is.matrix(samples)) {
    if (is.null(chains)) {
      if (nrow(samples) %% chains) {
        stop("There are ", nrow(samples), " iterations so there cannot be ",
          chains, " chains.")
      }
      
      i <- nrow(samples) / chains
      a <- array(dim = c(i, chains, ncol(samples)))
      for (ch in 1:chains) {
        a[, ch, ] <- samples[(i * (ch - 1) + 1):(i * ch), , drop = FALSE]
      }
      
      if (!is.null(nms)) dimnames(a)[3] <- list(parameters = nms)
      first_class(a) <- "nauf.mcmc"
      return(a)
    }
    
  } else if (is.array(samples) && length(dim(samples)) == 3) {
    if (!is.null(nms)) dimnames(samples)[3] <- list(parameters = nms)
    first_class(samples) <- "nauf.mcmc"
    return(samples)
  }
  
  stop("'samples' must be an array with three dimensions ",
    "[iterations, chains, parameters] or a matrix [iterations, parameters] ",
    "with 'chains' specified and not permuted.")
}


#' @export
summary.nauf.mcmc <- function(object, probs = c(0.025, 0.5, 0.975), psign = TRUE,
                              ...) {
  summ <- rstan::monitor(object, warmup = 0, probs = probs, print = FALSE)
  
  nc <- ncol(summ)
  colnames(summ)[c(1:3, nc - 1)] <- c("Mean", "MCE", "SD", "Neff")
  summ <- summ[, c(1, 3:(nc - 1), 2, nc), drop = FALSE]
  
  if (length(w <- which(colnames(summ) == "50%"))) {
    colnames(summ)[w] <- "Median"
    summ <- summ[, c(1:2, w, setdiff(3:nc, w)), drop = FALSE]
  }
  
  summ <- as.data.frame(summ)
  rownames(summ) <- dimnames(object)[[3]]
  
  if (psign) {
    summ[["P(Sign(Mean))"]] <- apply(object, 3, psgn)
  }
  
  first_class(summ) <- "summ.nauf.mcmc"
  return(summ)
}


#' @export
as.array.nauf.mcmc <- function(x, ...) {
  drop_class(x) <- "nauf.mcmc"
  return(x)
}


#' @export
as.matrix.nauf.mcmc <- function(x, ...) {
  x <- as.array(x)
  d <- dim(x)
  m <- apply(x, 3, FUN = function(y) y)
  if (!is.matrix(m)) m <- t(m)
  dimnames(m) <- list(iterations = NULL, parameters = dimnames(x)[[3]])
  return(m)
}


#' @export
as.data.frame.nauf.mcmc <- function(x, row.names = NULL, optional = TRUE, ...) {
  return(as.data.frame(summary(x, ...)))
}


#' @export
as.data.frame.summ.nauf.mcmc <- function(x, row.names = NULL, optional = TRUE,
                                         ...) {
  class(x) <- "data.frame"
  if (!is.null(row.names)) rownames(x) <- row.names
  return(x)
}


#' @export
print.summ.nauf.mcmc <- function(x, mce = FALSE, rhat = NULL, ...) {
  class(x) <- "data.frame"
  if (!mce) x <- x[-which(colnames(x) == "MCE")]
  w <- any(x$Rhat > 1.1)
  if (w) {
    warning("Some Rhat's are greater than 1.1; chains have not converged.")
  }
  if (isFALSE(rhat) || (is.null(rhat) && !w)) {
    x <- x[-which(colnames(x) == "Rhat")]
  }
  print(x, ...)
}


#' @export
print.nauf.mcmc <- function(x, ...) {
  print(summary(x, ...), ...)
}

