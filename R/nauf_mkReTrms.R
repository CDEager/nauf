

#' @export
nauf_mkReTrms <- function(object) {
  if (!is.nauf.frame(object) || !length(object@nauf$vars$bars)) {
    stop("'object' must be a nauf.frame with random effects", call. = FALSE)
  }
  
  bars <- lme4::findbars(object@nauf$formula)
  names(bars) <- lme4:::barnames(bars)
  term.names <- vapply(bars, safeDeparse, "")
  blist <- list()
  for (b in 1:length(bars)) blist[[b]] <- nauf_mkBlist(object, b)
  names(blist) <- names(bars)
  
  nl <- vapply(blist, `[[`, 0L, "nl")
  if (any(diff(nl) > 0)) {
      ord <- rev(order(nl))
      blist <- blist[ord]
      nl <- nl[ord]
      term.names <- term.names[ord]
  }
  
  Ztlist <- lapply(blist, `[[`, "sm")
  Zt <- do.call(Matrix::rBind, Ztlist)
  names(Ztlist) <- term.names
  q <- nrow(Zt)
  cnms <- lapply(blist, `[[`, "cnms")
  nc <- lengths(cnms)
  nth <- as.integer((nc * (nc + 1))/2)
  nb <- nc * nl
  if (sum(nb) != q) {
    stop(sprintf("total number of RE (%d) not equal to nrow(Zt) (%d)", 
      sum(nb), q))
  }
  boff <- cumsum(c(0L, nb))
  thoff <- cumsum(c(0L, nth))
  
  Lambdat <- Matrix::t(do.call(Matrix::sparseMatrix, do.call(Matrix::rBind,
    lapply(seq_along(blist), function(i) {
      mm <- matrix(seq_len(nb[i]), ncol = nc[i], byrow = TRUE)
      dd <- diag(nc[i])
      ltri <- lower.tri(dd, diag = TRUE)
      ii <- row(dd)[ltri]
      jj <- col(dd)[ltri]
      data.frame(i = as.vector(mm[, ii]) + boff[i], j = as.vector(mm[, jj]) +
        boff[i], x = as.double(rep.int(seq_along(ii), rep.int(nl[i],
        length(ii))) + thoff[i]))
    }))))
    
  thet <- numeric(sum(nth))
  ll <- list(Zt = Matrix::drop0(Zt), theta = thet, Lind = as.integer(Lambdat@x),
    Gp = unname(c(0L, cumsum(nb))))
  ll$lower <- -Inf * (thet + 1)
  ll$lower[unique(Matrix::diag(Lambdat))] <- 0
  ll$theta[] <- is.finite(ll$lower)
  Lambdat@x[] <- ll$theta[ll$Lind]
  ll$Lambdat <- Lambdat
  fl <- lapply(blist, `[[`, "ff")
  fnms <- names(fl)
  if (length(fnms) > length(ufn <- unique(fnms))) {
    fl <- fl[match(ufn, fnms)]
    asgn <- match(fnms, ufn)
  } else {
    asgn <- seq_along(fl)
  }
  names(fl) <- ufn
  fl <- do.call(data.frame, c(fl, check.names = FALSE))
  attr(fl, "assign") <- asgn
  ll$flist <- fl
  ll$cnms <- cnms
  ll$Ztlist <- Ztlist
  
  return(ll)
}


nauf_mkBlist <- function(nmf, bar) {
  ccn <- bar + 1
  bar <- nmf@nauf$bars[[bar]]
  mm <- nauf_mm(nmf, ccn)
  if (!bar$effects$intercept) mm <- mm[, -1, drop = FALSE]
  regf <- fac2sparse(factor(interaction(nmf[, bar$group$factors, drop = FALSE]),
    levels = bar$group$levels), to = "d", drop.unused.levels = FALSE)
  sm <- Matrix::KhatriRao(regf, t(mm))
  dimnames(sm) <- list(rep(levels(regf), each = ncol(mm)), rownames(mm))
  return(list(ff = regf, sm = sm, nl = nlevels(regf), cnms = colnames(mm)))
}

