

nauf.terms <- function(object, nauf_list = NULL) {
  if (!isS4(object) && inherits(object, "formula") &&
  inherits(object, "terms")) {
    cl <- class(object)
    class(object) <- c("nauf.terms", cl[cl != "nauf.terms"])
    currlist <- attr(object, "nauf")
    if (is.null(currlist) && is.null(nauf_list)) {
      stop("resulting nauf.terms object must have a valid 'nauf' attribute")
    }
    if (!is.null(nauf_list)) {
      attr(attr(object, "terms"), "nauf") <- nauf_list
    }
    attr(attr(object, "terms"), "predvars.fixed") <- nauf_list$predvars.fixed
    attr(attr(object, "terms"), "predvars.random") <- nauf_list$predvars.random
    return(object)
  }
  stop("'object' isn't S3 or doesn't inherit from formula and terms")
}


nauf.frame <- setClass("nauf.frame",
  contains = "data.frame",
  slots = list(nauf = "list")
)


nauf.glm <- function(object) {
  if (!isS4(object) && inherits(object, "lm")) {
    cl <- class(object)
    class(object) <- c("nauf.glm", cl[cl != "nauf.glm"])
    return(object)
  }
  stop("'object' isn't S3 or doesn't inherit from 'lm'")
}


#' @importClassesFrom lme4 lmerMod
nauf.lmerMod <- setClass("nauf.lmerMod", contains = "lmerMod")


#' @importClassesFrom lme4 glmerMod
nauf.glmerMod <- setClass("nauf.glmerMod", contains = "glmerMod")


nauf.grid <- function(object) {
  ### TODO: fix for when beta contains NAs; also fix in model.matrix generics
  if (is.nauf.glmerMod(object)) {
    dffun <- function(k, dfargs) NA
    dfargs <- list()
  } else if (!is.nauf.glm(object)) {
    dffun <- function(k, dfargs) dfargs$df
    dfargs <- list(df = mod$df)
  } else if (is.nauf.lmerMod(object)) {
    dffun <- function(k, dfargs) pbkrtest::Lb_ddf(k, dfargs$unadjV, dfargs$adjV)
    dfargs <- list(unadjV = lme4::vcov(mod), adjV = pbkrtest::vcovAdj(mod))
  } else {
    stop("must supply a nauf model")
  }
  
  mf <- model.frame(object)
  vars <- mf@nauf$vars
  lvs <- list()
  for (v in names(vars$uf)) {
    lvs[[v]] <- vars$uf[[v]][[1]]
    if (vars$hasna[v]) lvs[[v]] <- c(lvs[[v]], NA)
  }
  for (v in names(vars$of)) {
    lvs[[v]] <- vars$of[[v]]$levels
  }
  
  gmf <- expand.grid(lvs)
  
  nus <- list()
  for (v in names(vars$nu)) {
    for (j in names(vars$nu[[v]])) {
      lvs[[j]] <- vars$nu[[v]][[j]]
      nus[[j]] <- vars$nu[[v]][[j]]
    }
  }
  
  if (length(nus)) {
    p <- vars$predvars.nu
    nus <- eval(p, envir = nus)
    for (j in 1:length(vars$nu)) {
      if (is.matrix(nus[[j]])) {
        gmf[[names(vars$nu)[j]]] <- matrix(as.vector(nus[[j]]), nrow(gmf),
          ncol(nus[[j]]), byrow = TRUE)
      } else {
        gmf[[names(vars$nu)[j]]] <- nus[[j]]
      }
    }
  }
  
  attr(gmf, "terms") <- NULL
  mf@frame <- gmf
  gmm <- nauf_model.matrix(mf)
  asgn <- attr(model.matrix(object), "assign")
  summ <- summary(object)
  
  rgargs <- list(
    model.info = list(
      call = summ$call,
      terms = terms(object),
      xlev = lvs[sapply(lvs, length) > 1]),
    roles = list(
      predictors = names(lvs),
      responses = character(),
      multresp = character()),
    grid = gmf,
    levels = lvs,
    matlevs = list(),
    linfct = gmm,
    bhat = summ$coefficients[, 1],
    nbasis = matrix(),
    V = vcov(object),
    dffun = dffun,
    dfargs = dfargs,
    misc = list(
      estName = "prediction",
      estType = "prediction",
      infer = c(FALSE, FALSE),
      level = 0.95,
      adjust = "none",
      famSize = ncol(gmf) - 1,
      avgd.over = character(),
      assign = asgn)
    post.beta = matrix()
  )

  d <- data.frame(y = 1:5, x = c(2, 5, 3, 6, 1))
  rg <- lm(y ~ x, d)
  rg <- lsmeans::ref.grid(rg, data = d)
  for (i in names(rgargs)) {
    slot(rg, i) <- rgargs[[i]]
  }
  
  nrg <- list(ref.grid = rg)
  class(nrg) <- c("nauf.grid", "list")
  
  return(nrg)
}

