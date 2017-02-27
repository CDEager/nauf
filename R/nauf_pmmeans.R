
#' Create a reference grid with \code{NA} unordered factors.
#'
#' \code{nauf_grid} creates a reference grid for a model fit with
#' \code{\link{nauf_reg}} which can be used to obtain predicted marginal
#' means and pairwise comparisons for factors with \code{\link{nauf_pmmeans}}.
#' See \code{\link[lsmeans]{ref.grid}} for details of reference grid properties.
#'
#' A reference grid creates a data frame which contains all possible
#' combinations of the factors in a regression model, holding all covariates
#' at their mean values.  There are many options for
#' \code{\link[lsmeans]{ref.grid}} which are not currently supported for
#' \code{nauf} models.  The main functionality which is not currently supported
#' is that the reference grid cannot be created specifying certain levels for
#' factors (i.e. the \code{at} argument), and the behavior for polynomial
#' covariates has not been tested.  A direct call to
#' \code{\link[lsmeans]{ref.grid}} will result in warnings (or possibly errors),
#' and inference made with the resulting object will be misleading and/or
#' incorrect.  Only \code{nauf_grid} should be used.
#'
#' The reference grid returned by \code{nauf_grid} contains combinations of
#' factors which are not actually possible given the data set.  For example,
#' if factor \code{f1} has levels \code{A} and \code{B}, and factor \code{f2}
#' is \code{NA} when \code{f1 = A}, and takes values \code{C} and \code{D}
#' when \code{f1 = B}, the reference grid will still contain the combinations
#' \code{f1 = A, f2 = C}; \code{f1 = A, f2 = D}; and \code{f1 = B, f2 = NA},
#' even though these combinations are not possible.  This is because it is
#' impossible to know without the researcher's knowledge which combinations
#' make sense.  In many cases, this is inconsequential for the computation
#' of predicted marginal means, since the coding of unordered factors in
#' \code{nauf} regressions will average over the effects.  In cases where
#' these rows in the reference grid will cause invalid estimates and pairwise
#' comparisons, subsetting arguments can be used in the call to
#' \code{\link{nauf_pmmeans}} to ensure only the correct subsets are considered.
#'
#' @param mod A regression model fit with \code{\link{nauf_reg}}.
#'
#' @return A \code{\link[lsmeans]{ref.grid-class}} object which includes
#'   \code{NA} values for relevant unordered factors in the grid expansion.
#'
#' @examples
#' \dontrun{
#' mod <- nauf_reg(y ~ f1 + f2 + x)
#' rg <- nauf_grid(mod)
#' }
#'
#' @seealso \code{\link{nauf_pmmeans}}
#'
#' @importFrom utils combn
#' @export
nauf_grid <- function(mod) {
  if (!is.nauf(mod) || !inherits(mod, "lm")) {
    stop("must supply a nauf model")
  }

  asgn <- attr(nauf_model_matrix(mod), "assign")
  mtnr <- stats::delete.response(mod$terms)
  dc <- attr(mtnr, "dataClasses")
  cn <- attr(mod$terms, "factors")
  cn <- rownames(cn)[rowSums(cn) > 0]
  dc <- dc[names(dc) %in% cn]
  lvs <- xlvs <- mod$xlevels
  hasna <- attr(mod$terms, "hasna")
  dat <- mod$model
  for (n in names(hasna)[hasna]) {
    lvs[[n]] <- c(lvs[[n]], NA)
    dat[, n] <- addNA(dat[, n])
  }
  for (n in names(dc)[dc == "numeric"]) {
    lvs[[n]] <- mean(dat[, n])
  }
  g <- expand.grid(lvs)
  counts <- data.frame(xtabs(~ ., dat[, names(lvs)[sapply(lvs, length) > 1]]))
  g[, ncol(g) + 1] <- counts$Freq
  colnames(g)[ncol(g)] <- ".wgt."

  mf <- nauf_on(g)
  attr(mf, "terms") <- mtnr
  mefc <- attr(mtnr, "mefc")
  for (f in names(mefc)) {
    contrasts(mf[, f]) <- mefc[[f]]$contrasts
  }
  mm <- nauf_model_matrix(mf)
  summ <- summary(mod)

  rgargs <- list(
    model.info = list(
      call = mod$call,
      terms = mod$terms,
      xlev = mod$xlevels),
    roles = list(
      predictors = all.vars(mtnr),
      responses = character(),
      multresp = character()),
    grid = g,
    levels = lvs,
    matlevs = list(),
    linfct = mm,
    bhat = summ$coefficients[, 1],
    nbasis = matrix(),
    V = summ$cov,
    dffun = function(k, dfargs) dfargs$df,
    dfargs = list(df = mod$df),
    misc = list(
      estName = "prediction",
      estType = "prediction",
      infer = c(FALSE, FALSE),
      level = 0.95,
      adjust = "none",
      famSize = ncol(g) - 1,
      avgd.over = character(),
      assign = asgn),
    post.beta = matrix())

  d <- data.frame(y = 1:5, x = c(2,5,3,6,1))
  rg <- lm(y ~ x, d)
  rg <- lsmeans::ref.grid(rg, data = d)
  for (i in names(rgargs)) {
    slot(rg, i) <- rgargs[[i]]
  }

  return(rg)
}


#' Predicted marginal means with \code{NA} unordered factor values.
#'
#' \code{nauf_pmmeans} computes predicted marginal means (also called
#' least-squares means, a term which is less generally applicable) for a factor
#' or an interaction involving factors in a model fit with \code{\link{nauf_reg}},
#' and optionally pairwise comparisons of these group means.  There are several
#' subsetting arguments which allow for the predictions to be generated over
#' only the relevant subset of the reference grid created by
#' \code{\link{nauf_grid}} (see 'Details').
#'
#' First a description for the estimates generated when all elements of
#' \code{vars} are factors is given, and then the covariate estimates are
#' described.
#'
#' Because, by design, certain factors in a regression fit with
#' \code{\link{nauf_reg}} are applicable only within subsets of the data,
#' using \code{\link[lsmeans]{lsmeans}} directly will generate and average
#' over estimates which are invalid.  \code{nauf_pmmeans} allows the user
#' to generate estimates for a factor (or factors) for only the subsets of the
#' data where the estimate makes sense.  The first argument, \code{object}, can
#' be either a model object fit with \code{\link{nauf_reg}} or a reference grid
#' made with \code{\link{nauf_grid}}.  If a model object is provided, the
#' funciton first calls \code{\link{nauf_grid}} to create a reference grid.
#' Depending on the type and complexity of the model, creating the reference
#' grid can be time-consuming, and so it is recommended that a reference
#' grid be created first, and used as the \code{object} argument, especially
#' when multiple calls to \code{nauf_pmmeans} will be made for the same model.
#'
#' There are four ways to specify a subset of the reference grid to use in the
#' the computation of predicted marginal means, all of which default to
#' \code{NULL}, which indicates that the entire reference grid should be used.
#' If any of the factors specified in \code{vars} have \code{NA} values,
#' this will lead to incorrect estimates.
#'
#' The \code{keep_level} argument, if specified, should be a named list
#' of character vectors, where the names of the vectors are the names of factors
#' in the regression, and the elements of the vectors are levels which should
#' be kept in the reference grid.  For each element in \code{keep_level},
#' any row which has a value for that factor (including \code{NA}) which is
#' not in the chacter vector is dropped from the reference grid.  The
#' \code{drop_level} argument works in the opposite way: any row which contains
#' the specified levels is dropped from the reference grid.
#'
#' The \code{keep_group} argument, if specified, should be a list where each
#' element is a named list specifying the factor levels which define a group
#' of factor combinations which should be kept.  For each group entry (that is,
#' each element in the \code{keep_group} list), any row in the reference grid
#' which matches any combination of the factor levels in the group entry is
#' kept.  When there are multiple group entries in \code{keep_group}, any
#' reference grid row which matches at least one group is kept.  The
#' \code{drop_group} argument works in the opposte way: any row in the
#' reference grid which matches at least one group is dropped.
#'
#' As an example, consider a situation where factor \code{f1} takes levels
#' \code{A} and \code{B}, and factor \code{f2} takes levels \code{C} and
#' \code{D} when \code{f1 = A}, but \code{f2 = NA} whenver \code{f1 = B}.  That
#' is, \code{f2} is only contrastive for the group \code{f1 = A}.  In this
#' case, calling \code{nauf_pmmeans} with \code{vars = "f2"}, but not
#' specifying any of the subsetting arguments, would generate mean estimates
#' for \code{f2} averaging over the groups \code{f1 = B, f2 = C} and
#' \code{f1 = B, f2 = D}, which are not real groups.  To avoid this, there
#' are several possible ways of subsetting.  You could specify
#' \code{keep_level = list(f1 = "A")}, or \code{drop_level = list(f1 = "B")},
#' or \code{keep_group = list(list(f1 = "A"))}, or
#' \code{drop_group = list(list(f1 = "B"))}.  Any of these specifications
#' would result in two estimates, one for \code{f2 = C} which corresponds
#' to \code{f1 = A, f2 = C}, and one for \code{f2 = D} which corresponds
#' to \code{f1 = A, f2 = D}.
#'
#' Continuing with this example, the correct estimates for \code{f1} will
#' depend on the interpretation of \code{f2 = NA}.  If \code{f2 = NA} is
#' similar in interpretation to \code{f2 = C}, but is coded as \code{NA}
#' because \code{f1 = B, f2 = D} is not possible, then it would make sense
#' to compare the groups \code{f1 = A, f2 = C} and \code{f1 = B, f2 = NA}.
#' This is most easily accomplished by using the \code{drop_group} argument:
#' \code{drop_group = list(list(f1 = "B", f2 = c("C", "D")))}.  This doesn't
#' technically drop the non-existent groups \code{f1 = A, f2 = NA}, but this
#' won't affect the end result in this case because of the way \code{NA} values
#' are coded.  To explicitly drop this group (or similar groups when you aren't
#' sure whether leaving them in affects the outcome), you can do
#' \code{drop_group = list(list(f1 = "B", f2 = c("C", "D")), list(f1 = "A",
#' f2 = NA))}.  If \code{f2 = NA} cannot be interpreted as being similar to
#' to either \code{f2 = C} or \code{f2 = D}, then it makes sense to compare
#' \code{f1 = A}, averaging over \code{f1 = A, f2 = C} and \code{f1 = A, f2 = D},
#' to \code{f1 = B, f2 = NA}.  In this case, no subsetting arguments are
#' required (but may be used).
#'
#' The recommended use of the subsetting arguments is to create a list of
#' of explicit group definitions for \code{keep_group} as its own object, and
#' then use the relevant elements of that list in calls to \code{nauf_pmmeans}.
#'
#' If \code{vars} contains a combination of factors and covariates, then the
#' estimates generated are the slope for the full interaction of all of the
#' covariates in \code{vars} for each group in the full interaction of all of
#' the factors in \code{vars}.  As such, an error is thrown if \code{vars}
#' contains only covariates but no factors, or if \code{vars} includes
#' covariates but the full interaction of the variables listed in \code{vars}
#' is not in the model.  So, for example, for factor \code{f1} and covariates
#' \code{x1} and \code{x2}, with the full interaction \code{f1 * x1 * x2} in
#' the model, \code{vars = c("f1", "x1")} would return estimates (and optionally
#' pairwise comparisons) for the effect of an increase of 1 unit in \code{x1}
#' for each level of \code{f1} (excluding \code{NA} if \code{f1} has \code{NA}
#' values).  Calling \code{nauf_pmmeans} with  \code{vars = c("f1", "x1", "x2")}
#' would return the effect of the \emph{interaction term only} for each level
#' of \code{f1}; that is, the effect estimates would not include the main
#' effects of \code{x1} and \code{x2}.
#'
#' @param object A regression model fit with \code{\link{nauf_reg}} or a
#'   reference grid created with \code{\link{nauf_grid}}.
#' @param vars A character vector specifying the names of variables included
#'   in the model for which estimates are to be generated for the full
#'   interaction term.  At least one must be a factor.  See 'Details'.
#' @param keep_level A named list of character vectors specifying the levels
#'   of factors which should be considered in the reference grid. See 'Details'.
#' @param drop_level A named list of character vectors specifying the levels
#'   of factors which should not be considered in the reference grid.  See
#'   'Details.'
#' @param keep_group A list of groups which should be considered in the
#'   reference grid.  Any row in the grid which does not match at least one
#'   group exactly is not considered.  See 'Details'.
#' @param drop_group A list of groups which should not be considered in the
#'   reference grid.  Any row in the grid which matches at least one group
#'   exactly is not considered.  See 'Details'.
#' @param pairwise A logical (default \code{TRUE}) indicating whether or not
#'   pairwise comparisons for the group means should be computed.  If the
#'   there is only one group mean, \code{pairwise} is ignored.
#'
#' @return An object of class \code{lsm.list} which inherits from \code{nauf}
#'   and contains one or more objects of class \code{lsmobj}.  See for
#'   \code{\link[lsmeans]{ref.grid-class}} and \code{\link[lsmeans]{summary}}
#'   for arguments to alter printing options, including the confidence level
#'   for means (default is 0.95), the p-value adjustment method for pairwise
#'   comparisons (default is "tukey"), and whether to give the results
#'   on the scale of the dependent variable, or whether to use the inverse
#'   link function (e.g. for logistic regression, whether results should be
#'   given in log-odds or probability). The object additionally has an attribute
#'   \code{specs} which is a list of the \code{vars}, \code{keep_level},
#'   \code{drop_level}, \code{keep_group}, and \code{drop_group} arguments
#'   passed in the call to \code{nauf_pmmeans}.
#'
#' @examples
#' \dontrun{
#' # see 'Details' section for a description of the sample problem
#' # the recommended usage is demonstrated here
#' mod <- nauf_reg(y ~ f1 + f2, data = mydata)
#' rg <- nauf_grid(mod)
#' groups <- list()
#'
#' groups$f2 <- list(list(f1 = "A", f2 = c("C", "D")))
#' f2_pmm <- nauf_pmmeans(rg, vars = "f2", keep_group = groups$f2)
#'
#' # for the case where f2 = NA is similar in interpretation to f2 = C
#' groups$f1 <- list(
#'   A = list(f1 = "A", f2 = "C"),
#'   B = list(f1 = "B", f2 = NA))
#' f1_pmm <- nauf_pmmeans(rg, vars = "f1", keep_group = groups$f1)
#'
#' # for the case where f2 = NA is not similar to f2 = C or f2 = D
#' groups$f1 <- list(
#'   A = list(f1 = "A", f2 = c("C", "D")),
#'   B = list(f1 = "B", f2 = NA))
#' f1_pmm <- nauf_pmmeans(rg, vars = "f1", keep_group = groups$f1)
#'
#' # to compare all three possible groups
#' groups$f1_f2 <- list(
#'   A = list(f1 = "A", f2 = c("C", "D")),
#'   B = list(f1 = "B", f2 = NA)
#' f1_f2_pmm <- nauf_pmmeans(rg, vars = c("f1", "f2"),
#'   keep_group = groups$f1_f2)
#' }
#'
#' @seealso \code{\link{nauf_grid}}
#'
#' @export
nauf_pmmeans <- function(object, vars, keep_level = NULL, drop_level = NULL,
                         keep_group = NULL, drop_group = NULL,
                         pairwise = TRUE) {

  if (inherits(object, "ref.grid")) {
    if (is.nauf(object@model.info$terms)) {
      rg <- object
    } else {
      stop("must supply a ref.grid or model object fit with nauf functions")
    }
  } else if (is.nauf(object) && inherits(object, "lm")) {
    rg <- nauf_grid(object)
  } else {
    stop("must supply a ref.grid or model object fit with nauf functions")
  }
  if (!length(vars) || !is.character(vars) ||
  !all(vars %in% names(rg@levels))) {
    stop("'vars' must be a character vector specifying factor(s) in the model")
  }

  facs <- unlist(lapply(rg@levels, function(x) length(x) > 1))
  covs <- names(facs)[!facs]
  facs <- names(facs)[facs]
  vars_covs <- vars[vars %in% covs]
  vars_facs <- vars[vars %in% facs]

  if (!length(vars_facs)) {
    stop("at least one element of 'vars' must be a factor in the model")
  }

  if (!is.null(keep_level)) {
    if (!is.list(keep_level) || !all(names(keep_level) %in% facs)) {
      stop("'keep_level' must be a named list of factor levels to keep_level")
    }
    k <- rep(TRUE, nrow(rg@grid))
    for (n in names(keep_level)) {
      if (!all(keep_level[[n]] %in% rg@levels[[n]])) {
        stop("'keep_level' contains invalid factor levels")
      }
      k <- k & rg@grid[, n] %in% keep_level[[n]]
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }

  if (!is.null(drop_level)) {
    if (!is.list(drop_level) || !all(names(drop_level) %in% facs)) {
      stop("'drop_level' must be a named list of factor levels to drop_level")
    }
    k <- rep(TRUE, nrow(rg@grid))
    for (n in names(drop_level)) {
      if (!all(drop_level[[n]] %in% rg@levels[[n]])) {
        stop("'drop_level' contains invalid factor levels")
      }
      k <- k & !(rg@grid[, n] %in% drop_level[[n]])
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }

  if (!is.null(keep_group)) {
    if (!is.list(keep_group) || is.data.frame(keep_group) ||
    !all(unlist(lapply(keep_group, function(x) is.list(x) &
    !is.data.frame(x))))) {
      stop("'keep_group' must be a list of groups expressed as ",
        "named lists of factor levels")
    }
    k <- rep(FALSE, nrow(rg@grid))
    for (g in keep_group) {
      if (!all(names(g) %in% facs) || !all(sapply(names(g),
      function(x) all(g[[x]] %in% rg@levels[[x]])))) {
        stop("'keep_group' must be a list of groups expressed as ",
          "named lists of factor levels")
      }
      k <- k | apply(sapply(names(g),
        function(x) rg@grid[, x] %in% g[[x]]), 1, all)
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }

  if (!is.null(drop_group)) {
    if (!is.list(drop_group) || is.data.frame(drop_group) ||
    !all(unlist(lapply(drop_group, function(x) is.list(x) &
    !is.data.frame(x))))) {
      stop("'drop_group' must be a list of groups expressed as ",
        "named lists of factor levels")
    }
    k <- rep(FALSE, nrow(rg@grid))
    for (g in drop_group) {
      if (!all(names(g) %in% facs) || !all(sapply(names(g),
      function(x) all(g[[x]] %in% rg@levels[[x]])))) {
        stop("'drop_group' must be a list of groups expressed as ",
          "named lists of factor levels")
      }
      k <- k | apply(sapply(names(g),
        function(x) !(rg@grid[, x] %in% g[[x]])), 1, any)
    }
    rg@grid <- rg@grid[k, , drop = FALSE]
    rg@linfct <- rg@linfct[k, , drop = FALSE]
  }

  if (length(vars_covs)) {
    cn <- paste(paste(vars_covs, collapse = ":"), ":", paste(vars_facs,
      collapse = ","), sep = "")

    mt <- rg@model.info$terms
    ni <- attr(mt, "factors")
    ni <- which(apply(ni, 2, function(x) all(x[vars] > 0) & !any(
      x[!(names(x) %in% vars)] > 0)))
    if (!length(ni)) {
      stop("when 'vars' includes covariates, the full interaction term ",
        "must be in the model")
    }
    ni <- which(rg@misc$assign != ni)

    mf <- nauf_on(rg@grid)
    attr(mf, "terms") <- stats::delete.response(mt)
    mf[, vars_covs] <- 1
    mf[, covs[!(covs %in% vars_covs)]] <- 0
    mefc <- attr(attr(mf, "terms"), "mefc")
    for (f in names(mefc)) {
      contrasts(mf[, f]) <- mefc[[f]]$contrasts
    }

    mm <- nauf_model_matrix(mf)
    mm[, ni] <- 0
    rg@linfct <- mm

  } else {
    cn <- paste(vars_facs, collapse = ",")
  }

  estgrid <- expand.grid(rg@levels[vars_facs])
  if (length(vars_covs)) {
    estgrid <- estgrid[apply(!is.na(estgrid), 1, all), , drop = FALSE]
  }
  k <- rep(TRUE, nrow(estgrid))
  est <- list()
  for (i in 1:nrow(estgrid)) {
    est[[i]] <- rep(TRUE, nrow(rg@grid))
    for (j in vars_facs) {
      est[[i]] <- est[[i]] & (rg@grid[, j] %in% estgrid[i, j])
    }
    if ((n <- sum(est[[i]]))) {
      est[[i]] <- est[[i]] / n
    } else {
      k[i] <- FALSE
    }
  }
  if (sum(k) < 1) {
    stop("invalid specification; must be at least one group")
  }
  est <- est[k]
  estgrid <- estgrid[k, , drop = FALSE]
  names(est) <- apply(estgrid, 1, function(x) paste(x, collapse = ","))
  est.lsm <- lsmeans::contrast(rg, est)
  est.lsm@roles$predictors <- vars
  colnames(est.lsm@grid) <- cn
  est.lsm@misc$estName <- "pmmean"
  est.lsm@misc$infer <- c(TRUE, FALSE)
  est.lsm@misc$famSize <- nrow(est.lsm@grid)

  res <- list(pmmeans = est.lsm)
  class(res) <- c("lsm.list", "list")

  if (pairwise & sum(k) > 1) {
    combs <- utils::combn(length(est), 2)
    contr <- as.list(as.data.frame(apply(combs, 2,
      function(x) est[[x[1]]] - est[[x[2]]])))
    names(contr) <- apply(combs, 2,
      function(x) paste(names(est)[x], collapse = " - "))
    contr <- lsmeans::contrast(rg, contr)
    contr@misc$estType <- "pairs"
    contr@misc$adjust <- "tukey"
    contr@misc$methDesc <- "pairwise differences"
    contr@misc$famSize <- est.lsm@misc$famSize
    res$contrasts <- contr
  }

  attr(res, "specs") <- list(
    vars = vars,
    keep_level = keep_level,
    keep_group = keep_group,
    drop_level = drop_level,
    drop_group = drop_group)

  return(nauf_on(res))
}

