

#' @export
nauf_model.frame <- function(formula, data = NULL, subset = NULL,
                             na.action = na.pass, drop.unused.levels = TRUE,
                             xlev = NULL, contrasts = NULL, ncs_scale = NULL,
                             ...) {
  # Ignore na.action, contrasts, drop.unused.levels, and xlev
  mc <- match.call()
  mc$na.action <- na.pass
  mc$drop.unused.levels <- TRUE
  mc$xlev <- NULL
  mc$contrasts <- NULL
  if (!isTRUE(all.equal(na.action, na.pass))) {
    warning("Ignoring 'na.action', must be na.pass")
  }
  if (!isTRUE(drop.unused.levels)) {
    warning("Ignoring 'drop.unused.levels', must be TRUE")
  }
  if (!is.null(xlev)) {
    warning("Ignoring 'xlev', must be NULL")
  }
  if (!is.null(contrasts)) {
    warning("Ignoring 'contrasts', must be NULL")
  }
  
  standardized_scale <- attr(formula, "standardized.scale")
  if (!is.null(standardized_scale)) {
    if (!is.null(ncs_scale)) {
      warning("Ignoring 'ncs_scale' because 'formula' comes from a",
      " standardized object")
    }
    ncs <- standardized_scale
  } else if (!is.null(ncs_scale)) {
    ncs <- ncs_scale
  } else {
    ncs <- 1
  }
  if (!is.scalar(ncs, 1)) {
    stop("The scale for sum contrasts must be a single positive number")
  }
  
  formula <- stats::formula(formula)
  class(formula) <- "formula"
  bars <- lme4::findbars(formula)
  sb_form <- lme4::subbars(formula)
  if (stats::is.empty.model(sb_form)) stop("There are no predictors in 'formula'")
  fe_form <- stats::terms(lme4::nobars(formula))
  if (!attr(fe_form, "response")) stop("'formula' must have a response")
  if (!attr(fe_form, "intercept")) {
    stop("There must be a fixed effects intercept")
  }
  fe_form <- stats::formula(stats::delete.response(fe_form))
  groups <- check_groups(formula)
  
  if (!is.data.frame(data)) stop("'data' must be a data.frame")
  
  mc$formula <- sb_form
  mc["ncs_scale"] <- NULL
  mc[[1]] <- quote(stats::model.frame)
  mf <- eval(mc, parent.frame())
  cnms <- colnames(mf)
  extras <- find_extras(mf)
  mt <- attr(mf, "terms")
  
  fmat <- attr(mt, "factors")
  fmat <- fmat[!(rownames(fmat) %in% groups), , drop = FALSE]
  if (any(fmat > 1)) {
    warning("'nauf' has not been tested for models that violate the ",
      "interaction hierarchy")
  }
  
  for (g in groups) mf[[g]] <- factor(mf[[g]], ordered = FALSE)
  for (j in which(!extras)[-1]) mf[[j]] <- charlogbin_to_uf(mf[[j]])
  
  uf <- sapply(mf, is.uf) & !extras
  hasna <- sapply(mf, anyNA)
  if (hasna[1] || any(hasna & !uf)) {
    stop("Only unordered factor predictors and random effects grouping factors",
      " can have NA values")
  }
  uf[1] <- FALSE
  uf[groups] <- FALSE
  of <- sapply(mf, is.ordered) & !extras
  mat <- sapply(mf, is.matrix) & !extras
  num <- sapply(mf, is.numeric) & !mat & !extras
  nauf <- uf & hasna
  rgrp <- cnms %in% groups
  names(rgrp) <- cnms
  nagrp <- rgrp & hasna
  num[1] <- mat[1] <- of[1] <- FALSE
  
  vars <- list()
  
  vars$groups <- list()
  for (j in groups) {
    vars$groups[[j]] <- levels(mf[[j]])
  }
  
  vars$uf <- list()
  for (j in cnms[uf]) {
    mf[[j]] <- standardize::named_contr_sum(mf[[j]], ncs, FALSE)
    vars$uf[[j]] <- list()
    vars$uf[[j]][[1]] <- levels(mf[[j]])
  }
  
  vars$of <- list()
  for (j in cnms[of]) {
    vars$of[[j]] <- list(levels = levels(mf[[j]]),
      contrasts = contrasts(mf[[j]]))
  }
  
  vars$num <- list()
  for (j in cnms[num]) {
    vars$num[[j]] <- mean(mf[[j]])
  }
  
  vars$mat <- list()
  for (j in cnms[mat]) {
    vars$mat[[j]] <- colMeans(mf[[j]])
  }
  
  vars$extras <- cnms[extras]
  
  temp <- contrast_changes(fe_form, mf, vars$uf)
  cc <- temp$cc
  vars$uf <- temp$uf
  for (b in bars) {
    temp <- contrast_changes(b, mf, vars$uf, cc)
    cc <- temp$cc
    vars$uf <- temp$uf
  }
  
  attr(mt, "dataClasses")[names(vars$uf)] <- "factor"
  attr(mt, "dataClasses")[names(vars$groups)] <- "factor"
  attr(mt, "nauf.info") <- list(
    resp = cnms[1],
    groups = vars$groups,
    uf = vars$uf,
    of = vars$of,
    num = vars$num,
    mat = vars$mat,
    extras = vars$extras,
    cc = cc,
    hasna = hasna,
    ncs_scale = ncs)
  class(mf) <- c("data.frame", "nauf.frame")
  class(mt) <- c("nauf.terms", "terms", "formula")
  class(formula) <- c("nauf.formula", "formula")
  attr(mf, "terms") <- mt
  attr(mf, "formula") <- formula
  
  return(mf)
}


contrast_changes <- function(form, mf, lvs, cc = NULL) {
  # TODO (CDEager): rewrite this so it requires less copying
  ccn <- length(cc) + 1
  if (re <- ccn > 1) {
    cc[[ccn]] <- list()
  } else {
    cc <- list(list())
  }
  
  mt <- attr(mf, "terms")
  allmat <- attr(mt, "factors")
  rn <- rownames(allmat)
  
  if (re) {
    group <- varnms(barform(form, 3))
    mf <- mf[!rowna(mf[, group, drop = FALSE]), , drop = FALSE]
    form <- barform(form, 2)
    if (!length(attr(stats::terms(form), "factors"))) {
      return(list(uf = lvs, cc = cc))
    }
  }
  
  fmat <- attr(stats::terms(form), "factors") > 0
  rn <- rownames(fmat)
  mf <- mf[, rn, drop = FALSE]
  uf <- sapply(mf, is.uf)
  hasna <- sapply(mf, anyNA)
  nauf <- uf & hasna
  ufmat <- fmat[uf, , drop = FALSE]
  naufmat <- fmat[nauf, , drop = FALSE]
  check_inter <- which(colSums(ufmat) > 1 & colSums(naufmat) > 0)
  
  if (re) {
    check_main <- rn[uf]
    check_main <- check_main[check_main %in% colnames(fmat)]
    for (j in check_main) {
      jlvs <- levels(droplevels(mf[[j]]))
      cj <- in_list(jlvs, lvs[[j]])
      if (!cj) {
        cj <- length(lvs[[j]]) + 1
        lvs[[j]][[cj]] <- jlvs
      }
      if (cj > 1) cc[[ccn]][[j]] <- cj
    }
  }
  
  if (length(check_inter)) {
    check_inter <- unique(lapply(check_inter, function(x) {
      sort(rownames(ufmat)[ufmat[, x] > 0])
    }))

    for (i in check_inter) {
      checki <- nauf_interaction(mf, i)
      if (checki$changed) {
        nn <- paste(i, collapse = ":")
        cc[[ccn]][[nn]] <- list(factors = numeric(),
          assign = which((colSums(ufmat) == length(i)) & apply(
            ufmat[i, , drop = FALSE], 2, all)))
        for (j in i) {
          jlvs <- checki$levels[[j]]
          cj <- in_list(jlvs, lvs[[j]])
          if (!cj) {
            cj <- length(lvs[[j]]) + 1
            lvs[[j]][[cj]] <- jlvs
          }
          cc[[ccn]][[nn]]$factors[j] <- cj
        }
      }
    }
  }

  return(list(uf = lvs, cc = cc))
}

