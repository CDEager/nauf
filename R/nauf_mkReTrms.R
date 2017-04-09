

nauf_mkReTrms <- function(fr, lvs = NULL) {
  # based on lme4::mkReTrms
  if (!inherits(fr, "nauf.frame")) {
    stop("'fr' was not created with nauf_model.frame")
  }
  bars <- lme4::findbars(attr(fr, "formula"))
  if (!length(bars)) {
    stop("No random effects terms specified in formula", call. = FALSE)
  }
  stopifnot(is.list(bars), vapply(bars, is.language, NA), inherits(fr, 
    "data.frame"))
  
  names(bars) <- lme4_barnames(bars)
  term.names <- vapply(bars, lme4_safeDeparse, "")
  blist <- nauf_mkBlist(bars, fr, lvs)
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
  nth <- as.integer((nc * (nc + 1)) / 2)
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


nauf_mkBlist <- function (bars, fr, lvs) {
  blist <- list()
  
  for (b in 1:length(bars)) {
    vars <- varnms(barform(bars[[b]], 3))
    ff <- interaction(fr[, vars, drop = FALSE])
    
    if (is.null(lvs)) {
      ff <- droplevels(ff)
      if (all(is.na(ff))) {
        stop("Invalid grouping factor specification, ", deparse(bars[[b]][[3]]),
          call. = FALSE)
      }
      
    } else {  # implies predict method with new data
      ff <- factor(ff, levels = lvs[[paste(vars, collapse = ":")]],
        ordered = FALSE)
    }
    
    mm <- nauf_mm(fr, b + 1)
    sm <- Matrix::KhatriRao(Matrix::fac2sparse(ff, to = "d",
      drop.unused.levels = FALSE), t(mm))
    dimnames(sm) <- list(rep(levels(ff), each = ncol(mm)), rownames(mm))
    
    blist[[b]] <- list(ff = ff, sm = sm, nl = nlevels(ff), cnms = colnames(mm))
  }
  
  names(blist) <- names(bars)
  
  return(blist)
}

