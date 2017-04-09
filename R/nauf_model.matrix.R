

#' @export
nauf_model.matrix <- function(object = NULL, data = NULL, ...) {
  mc <- match.call()

  if (inherits(object, "nauf.frame")) return(nauf_mm(object))
  if (inherits(data, "nauf.frame")) return(nauf_mm(data))
  if (inherits(object, "formula")) {
    names(mc)[2] <- "formula"
    mc[[1]] <- quote(nauf::nauf_model.frame)
    mf <- eval(mc, parent.frame())
    return(nauf_mm(mf))
  }
  
  stop("If 'object' is not a formula, then either 'object' or 'data' must be\n",
    "  a model frame created by nauf_model.frame")
}


nauf_mm <- function(mf, ccn = 1) {
  attr(mf, "na.action") <- "na.pass"
  formula <- attr(mf, "formula")
  mt <- attr(mf, "terms")
  info <- attr(mt, "nauf.info")
  ncs <- info$ncs_scale
  ufc <- info$uf
  cc <- info$cc[[ccn]]
  
  if (ccn == 1) {
    formula <- stats::delete.response(stats::terms(lme4::nobars(formula)))
  } else {
    formula <- stats::terms(barform(lme4::findbars(formula)[[ccn - 1]], 2))
    ccmain <- names(ufc)[names(ufc) %in% names(cc)]
    if (length(ccmain)) {
      torm <- which(names(cc) %in% ccmain)
      for (j in ccmain) {
        cj <- cc[[j]]
        lvs <- ufc[[j]][[cj]]
        mf[, j] <- factor(mf[, j], ordered = FALSE, levels = lvs)
        contr <- named_contr_sum(lvs, ncs)
        colnames(contr) <- paste0(".c", cj, ".", colnames(contr))
        contrasts(mf[, j]) <- contr
      }
      cc <- cc[-torm]
    }
  }

  mm <- stats::model.matrix(formula, mf)

  if (length(cc)) {
    mmlist <- list()
    cnms <- character()
    asgn <- attr(mm, "assign")
    asgn_cc <- sort(unique(unlist(lapply(cc, function(n) n$assign))))
    mmrm <- which(asgn %in% asgn_cc)
    asgn <- asgn[-mmrm]
    if (length(asgn)) {
      mmlist[[1]] <- mm[, -mmrm, drop = FALSE]
      cnms <- colnames(mmlist[[1]])
    }
    fmat <- attr(formula, "factors")

    for (i in cc) {
      uf <- names(i$factors)
      mfi <- mf
      for (u in uf) {
        cj <- i$factors[u]
        lvs <- varlist$uf[[u]][[cj]]
        mfi[, u] <- factor(mfi[, u], ordered = FALSE, levels = lvs)
        contr <- named_contr_sum(lvs, ncs)
        if (cj > 1) {
          colnames(contr) <- paste0(".c", cj, ".", colnames(contr))
        }
        contrasts(mfi[, u]) <- contr
      }

      imat <- fmat[, i$assign, drop = FALSE]
      main <- rownames(imat)[rowSums(imat) > 0]
      imat <- imat[main, , drop = FALSE]
      form <- paste("~", paste(main, collapse = "+"))
      for(j in 1:ncol(imat)) {
        form <- paste(form, paste(main[imat[, j] > 0], collapse = "*"),
          sep = "+")
      }
      ti <- stats::terms(stats::formula(form))
      asgn_ti <- which(colnames(attr(ti, "factors")) %in% colnames(imat))

      mmi <- model.matrix(ti, mfi)
      asgn_mmi <- attr(mmi, "assign")
      keep <- which(asgn_mmi %in% asgn_ti)
      mmi <- mmi[, keep, drop = FALSE]
      asgn_mmi <- asgn_mmi[keep]
      asgn_mmi <- as.numeric(factor(asgn_mmi))
      asgn <- c(asgn, i$assign[asgn_mmi])
      cnms <- c(cnms, colnames(mmi))

      mmlist[[length(mmlist) + 1]] <- mmi
    }

    mm <- do.call(cbind, args = mmlist)
    names(asgn) <- cnms
    asgn <- sort(asgn)
    mm <- mm[, names(asgn), drop = FALSE]
    names(asgn) <- NULL
    attr(mm, "assign") <- asgn
  }

  attr(mm, "contrasts") <- NULL
  mm[is.na(mm)] <- 0

  return(mm)
}

