

# na.action, contrasts, drop.unused.levels, and xlev are ignored
#' @export
nauf_model.frame <- function(formula, data, family = gaussian, sumcoef = 1,
                             na.action = "na.pass", contrasts = NULL,
                             drop.unused.levels = TRUE, xlev = NULL, ...) {
  mc <- match.call()
  
  if (is.character(formula)) formula <- stats::formula(formula)
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a valid formula")
  }
  attributes(formula) <- NULL
  if (!attr(stats::terms(formula), "response")) {
    stop("'formula' has no response")
  }
  formula[[length(formula)]] <- lme4::expandDoubleVerts(
    formula[[length(formula)]])
  
  if (!is.data.frame(data) || is.nauf.frame(data)) {
    stop("'data' must be a data.frame")
  }
  attr(data, "terms") <- NULL
  
  family <- get_family(family)
  
  if (!is.scalar(sumcoef, 1)) {
    stop("'sumcoef' must be a single positive number")
  }

  regfunc <- get_regfunc(formula, family)

  mfvars <- get_frame_vars(formula, data, family, sumcoef, ...)

  mt <- attr(stats::model.frame(lme4::subbars(formula), data,
    na.action = "na.pass", ...), "terms")
  nauf_list <- list(call = mc, formula = formula, family = family,
    regfunc = regfunc, sumcoef = sumcoef, vars = mfvars$vars)
  mt <- nauf.terms(mt, nauf_list)  ## update so predvars is done in nauf.terms

  mf <- nauf.frame(mfmt$frame, nauf = nauf_list)
  attr(mf, "terms") <- mt
  attr(mf, "formula") <- formula
  attr(mf, "offset") <- mc$offset

  return(mf)
}


get_frame_vars <- function(formula, data, family, sumcoef, ...) {
  if (is.mixed(formula)) {
    sb <- lme4::subbars(formula)
    bars <- lme4::findbars(formula)
    barinfo <- lapply(bars, get_bar_info)
    names(barinfo) <- barnames(bars)
    groups <- group_factors(bars)
    
    mf <- nauf_mf(formula = sb, data = data, sumcoef = sumcoef,
      groups = groups, ...)
    mt <- attr(mf, "terms")
    hasna <- attr(mt, "hasna")
    vars <- attr(mt, "varlist")
    
    nrf <- stats::delete.response(mt)
    fef <- lme4::nobars(nrf)
    ref <- lme4::subbars(lme4:::reOnly(nrf))
    nrf <- lme4::subbars(nrf)
    
    predvars <- attr(attr(model.frame(nrf, data, na.action = "na.pass"),
      "terms"), "predvars")
    predvars.fixed <- attr(attr(model.frame(fef, data, na.action = "na.pass"),
      "terms"), "predvars")
    predvars.random <- attr(attr(model.frame(ref, data, na.action = "na.pass"),
      "terms"), "predvars")
    
    if (length(vars$nu)) {
      nuf <- stats::formula(paste("~", paste(names(vars$nu), collapse = "+")))
      predvars.nu <- attr(attr(model.frame(nuf, data, na.action = "na.pass"),
        "terms"), "predvars")
    } else {
      predvars.nu <- NULL
    }
    
    for (b in 1:length(barinfo)) {
      dn <- barinfo[[b]]$group$datanames
      fs <- barinfo[[b]]$effects$terms
      attributes(fs) <- NULL
      db <- data[!rowna(data[, dn, drop = FALSE]), , drop = FALSE]
      vars <- attr(attr(nauf_mf(fs, data = db, sumcoef = sumcoef,
        groups = groups, varlist = vars), "terms"), "varlist")
    }
    
    vars$hasna <- hasna
    vars$bars <- barinfo
    
  } else {
    mf <- nauf_mf(formula = formula, data = data, sumcoef = sumcoef, ...)
    mt <- attr(mf, "terms")
    vars <- attr(mt, "varlist")
    vars$hasna <- attr(mt, "hasna")
    vars$bars <- NULL
    
    nrf <- stats::delete.response(mt)
    predvars <- attr(attr(model.frame(nrf, data, na.action = "na.pass"),
      "terms"), "predvars")
    predvars.fixed <- predvars
    predvars.random <- NULL
    
    if (length(vars$nu)) {
      nuf <- stats::formula(paste("~", paste(names(vars$nu), collapse = "+")))
      predvars.nu <- attr(attr(model.frame(nuf, data, na.action = "na.pass"),
        "terms"), "predvars")
    } else {
      predvars.nu <- NULL
    }
  }
  
  attr(mf, "terms") <- NULL
  vars$predvars <- predvars
  vars$predvars.fixed <- predvars.fixed
  vars$predvars.random <- predvars.random
  vars$predvars.nu <- predvars.nu
  
  return(list(frame = mf, vars = vars))
}


nauf_mf <- function(formula, data, na.action = "na.pass",
                    drop.unused.levels = TRUE, xlev = NULL,
                    sumcoef = 1, varlist = NULL, groups = NULL, ...) {

  re <- !is.null(varlist)
  attributes(formula) <- NULL
  mf <- stats::model.frame(formula, data, na.action = "na.pass",
    drop.unused.levels = TRUE, xlev = NULL, ...)
  mt <- attr(mf, "terms")
  
  v <- predictors(mt)
  v <- v[!(v %in% groups)]
  mf[, v] <- charlogbin_to_uf(mf[, v, drop = FALSE])
  
  if (re) {
    resp <- varlist$resp
    offs <- varlist$off
    regfs <- varlist$regf
    nus <- varlist$nu
    ufc <- varlist$uf
    ofc <- varlist$of
    cc <- varlist$cc
    ccn <- length(cc) + 1
    cc[[ccn]] <- list()
    
    for (g in names(regfs)[names(regfs) %in% colnames(mf)]) {
      mf[, g] <- factor(mf[, g], ordered = FALSE, levels = regfs[[g]]$levels)
    }
    
    # these for loops will only change mf[, j] in very rare circumstances
    # where there are structural problems with the data, i.e. when
    # a random effect group has NA values, and in the subset where they are
    # not NA, an ordered factor loses levels or a numeric variable becomes
    # binary (causing to be converted to uf above)
    for (j in names(ofc)[names(ofc) %in% v]) {
      if (!is.ordered(mf[, j]) || !isTRUE(all.equal(ofc[[j]]$levels,
      levels(mf[, j])))) {
        warning("In applicable subset of the random effects, ", j, "\n",
          "   has fewer levels than in the fixed effects")
      }
        mf[, j] <- factor(mf[, j], ordered = TRUE, levels = ofc[[j]]$levels)
        contrasts(mf[, j]) <- ofc[[j]]$contrasts
    }
    for (j in names(nus)[names(nus) %in% v]) {
      if (!is.numeric(mf[, j])) {
        mf[, j] <- as.numeric(as.character(mf[, j]))
        warning("In fixed effects, ", j, "has more than two values (numeric)\n",
          "   but in the applicable subset of the random effects, it does not.")
      }
    }
    
  } else {
    regfs <- nus <- mats <- ufc <- ofc <- list()
    cc <- list(list())
    ccn <- 1
    for (g in groups) {
      mf[, g] <- droplevels(factor(mf[, g], ordered = FALSE))
      regfs[[g]] <- list(levels = levels(mf[, g]), hasna = anyNA(mf[, g]))
    }
    ynms <- non_term_names(mt)
    resp <- colnames(mf)[1]
    
    # offs may contain other things besides offsets (e.g. weights, mustart)
    # at this point nothing is being done with them besides recording their
    # names.
    offs <- ynms[!(ynms %in% resp)]
  }
  
  uf <- sapply(mf, is.uf)[v]
  of <- sapply(mf, is.ordered)[v]
  nu <- !uf & !of
  hasna <- colna(mf)[v]
  nauf <- hasna & uf
  attr(mt, "dataClasses")[v[uf]] <- "factor"
  fmat <- attr(mt, "factors")

  if (any(hasna & !uf)) {
    stop("the following have NA values but are not unordered factors:\n",
      paste(v[hasna & !uf], collapse = "\n"))
  }
  if (any(fmat > 1)) {
    warning("nauf has not been tested for models that violate ",
      "the interaction hierarchy.")
  }
  if (!attr(mt, "intercept")) {
    warning("nauf has not been tested for zero-intercept models.")
  }

  if (!re) {
    for (j in v[uf]) {
      ufc[[j]] <- list()
      mf[, j] <- named_contr_sum(mf[, j], sumcoef, FALSE)
      ufc[[j]][[1]] <- levels(mf[, j])
    }
    
    ofcfun <- as.character(getOption("contrasts")[2])
    for (j in v[of]) {
      if (empty_cells(mf[, j])) {
        mf[, j] <- droplevels(mf[, j])
        warning("Dropping unused levels from ordered factor ", j)
      }
      
      # if you don't do this then the attributes aren't equal
      # in testing for the ordered factors
      if (is.null(attr(mf[, j], "contrasts"))) {
        contrasts(mf[, j]) <- ofcfun
      }
      attr(mf[, j], "contrasts") <- contrasts(mf[, j])
      rownames(attr(mf[, j], "contrasts")) <- levels(mf[, j])
      
      ofc[[j]] <- list(levels = levels(mf[, j]),
        contrasts = contrasts(mf[, j]))
    }
    
    for (j in v[nu]) {
      nus[[j]] <- get_num_trm_means(j, data)
    }
    
  } else {
    for (j in v[uf]) {
      mf[, j] <- named_contr_sum(mf[, j], sumcoef, FALSE)
      cj <- in_list(levels(mf[, j]), ufc[[j]])
      if (!cj) {
        ufc[[j]][[length(ufc[[j]]) + 1]] <- levels(mf[, j])
        cj <- length(ufc[[j]])
      }
      if (cj > 1) {
        cc[[ccn]][[j]] <- cj
      }
    }
  }

  if (any(nauf)) {
    vmat <- fmat[v, , drop = FALSE]
    ufmat <- vmat[v[uf], , drop = FALSE]
    naufmat <- ufmat[v[nauf], , drop = FALSE]
    has2uf <- colSums(ufmat) > 1
    hasnauf <- colSums(naufmat) > 0
    check_naufi <- which(has2uf & hasnauf)

    if (length(check_naufi)) {
      check_cc <- list()
      for (j in 1:length(check_naufi)) {
        # loop since using apply() could return matrix instead of list
        check_cc[[j]] <- sort(rownames(ufmat)[ufmat[, check_naufi[j]] > 0])
      }
      check_cc <- unique(check_cc)

      for (i in check_cc) {
        checki <- nauf_interaction(mf, i)
        if (checki$changed) {
          nn <- paste(i, collapse = ":")
          cc[[ccn]][[nn]] <- list(factors = list(),
            assign = which((colSums(ufmat) == length(i)) & apply(
              ufmat[i, , drop = FALSE], 2, all)))
          for (j in i) {
            lvs <- checki$levels[[j]]
            cj <- in_list(lvs, ufc[[j]])
            if (!cj) {
              ufc[[j]][[length(ufc[[j]]) + 1]] <- lvs
              cj <- length(ufc[[j]])
            }
            cc[[ccn]][[nn]]$factors[[j]] <- cj
          }
        }
      }
    }
  }
    
  if (!re) {
    attr(mt, "hasna") <- colna(mf)
  }

  varlist <- list(
    resp = resp,
    off = offs,
    regf = regfs,
    nu = nus,
    uf = ufc,
    of = ofc,
    cc = cc)
  
  attr(mt, "varlist") <- varlist
  attr(mt, "sumcoef") <- sumcoef
  attr(mf, "terms") <- mt
  
  return(mf)
}

