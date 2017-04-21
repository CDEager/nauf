
#' Create a model frame using \code{nauf} contrasts.
#'
#' \code{nauf_model.frame} creates a model frame which employs
#' \code{\link{nauf_contrasts}} for unordered factors.
#'
#' First, the default method for \code{\link[stats]{model.frame}} is called.
#' Then any variable in the resulting model frame that is not an unordered
#' factor but has only two unique non-\code{NA} values is coerced to an
#' unordered factor.  Unordered factors are then assigned contrasts with
#' \code{\link[standardize]{named_contr_sum}}, passing \code{ncs_scale} as the
#' function's \code{scale} argument.  Then, necessary contrast changes in
#' interaction terms and random effects slopes are determined as described in
#' \code{\link{nauf_contrasts}}.
#'
#' The recommended usage is to first \code{\link[standardize]{standardize}} the
#' regression variables, and then use the \code{formula} and \code{data}
#' elements in the resulting \code{standardized} object as arguments to
#' \code{nauf_model.frame}.  When this is done, \code{ncs_scale} is obtained
#' from the \code{standardized.scale} attribute of the \code{formula}, unless
#' \code{ncs_scale} is specified as a value which does not match the
#' \code{standardized} scale, in which case the explicitly specified
#' \code{ncs_scale} argument is used with a warning.  If
#' \code{\link[standardize]{standardize}} is not used prior to calling
#' \code{nauf_model.frame}, then \code{ncs_scale} defaults to \code{1} unless
#' explicitly specified in the function call, in which case the specified value
#' is used.
#'
#' Changes from the following default values are ignored with a warning:
#' \describe{
#'   \item{na.action = na.pass}{This default value is required in order for
#'     \code{NA} values to be treated as defined in
#'     \code{\link{nauf_contrasts}}.}
#'   \item{drop.unused.levels = TRUE}{This default value is set because
#'     \code{nauf_model.frame} assumes that \code{data} is not new data.  To
#'     create a \code{\linkS4class{nauf.frame}} with new data, the \code{terms}
#'     attribute of an already existing \code{\linkS4class{nauf.frame}} (which
#'     has class \code{\linkS4class{nauf.terms}}) can be used as the
#'     \code{formula} argument to \code{\link[stats]{model.frame}}.}
#'   \item{xlev = NULL}{This default is necessary for the same reasons as the
#'     default value for \code{drop.unused.levels}.}
#'   \item{contrasts = NULL}{For unordered factors, contrasts are automatically
#'     created with \code{\link[standardize]{named_contr_sum}}, as sum contrasts
#'     are necessary to implement \code{\link{nauf_contrasts}}.  To specify
#'     custom contrasts for ordered factors, the custom contrasts should be
#'     explicitly assigned to the ordered factor in \code{data} (this is
#'     automatically done if \code{\link[standardize]{standardize}} is used
#'     first as recommended).}
#' }
#'
#' @param formula,data,subset,... See \code{\link[stats]{model.frame}}.
#' @param na.action,drop.unused.levels,xlev,contrasts Changes from default
#'   values for these arguments are ignored with a warning.
#' @param ncs_scale A positive number passed as the \code{scale} argument to
#'   \code{\link[standardize]{named_contr_sum}} for all unordered factor
#'   contrasts.  The default is to first check whether \code{formula} comes from
#'   a \code{standardized} object returned by
#'   \code{\link[standardize]{standardize}}. If it is, then the \code{scale}
#'   argument from the \code{\link[standardize]{standardize}} call is used.  If
#'   it is not, then \code{ncs_scale} is set to \code{1}.  The value for
#'   \code{ncs_scale} can also be set explicitly.  If it is set explicitly and
#'   \code{formula} is from a \code{standardized} object with a different scale
#'   than the explicitly set value, then the explicitly set value
#'   is used and a warning is issued.
#'
#' @return A \code{\linkS4class{nauf.frame}}.
#'
#' @examples
#' dat <- plosives
#' dat$spont[dat$dialect == "Valladolid"] <- NA
#' form <- intdiff ~ voicing * dialect * spont +
#'   (1 + voicing * spont | speaker) + (1 + dialect | item)
#'
#' ## default behavior when standardize is not used
#' # defaults to ncs_scale = 1
#' mf <- nauf_model.frame(form, dat)
#'
#' # uses specified ncs_scale = 0.5
#' mf_0.5 <- nauf_model.frame(form, dat, ncs_scale = 0.5)
#'
#' ## standardize first (recommended use)
#' sdat <- standardize(form, dat)
#' sdat_0.5 <- standardize(form, dat, scale = 0.5)
#'
#' # uses ncs_scale = 1 from attr(sdat$formula, "standardized.scale")
#' mf_sdat <- nauf_model.frame(sdat$formula, sdat$data)
#'
#' # uses ncs_scale = 0.5 from attr(sdat_0.5$formula, "standardized.scale")
#' mf_sdat_0.5 <- nauf_model.frame(sdat_0.5$formula, sdat_0.5$data)
#'
#' \dontrun{
#' ## not recommended
#' # uses specified ncs_scale = 0.5 and issues a warning since
#' # attr(sdat$formula, "standardized.scale") = 1
#' mf_warning <- nauf_model.frame(sdat$formula, sdat$data, ncs_scale = 0.5)
#' }
#'
#' @seealso \code{\link{nauf_contrasts}} for a description of the contrasts
#'   applied to unordered factors, \code{\link{nauf_model.matrix}} for obtaining
#'   a fixed effects model matrix, and \code{\link{nauf_glFormula}} for
#'   obtaining both fixed effects and random effects model matrices.
#'
#' @export
nauf_model.frame <- function(formula, data = NULL, subset = NULL,
                             na.action = na.pass, drop.unused.levels = TRUE,
                             xlev = NULL, contrasts = NULL,
                             ncs_scale = attr(formula, "standardized.scale"),
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
  if (is.null(ncs <- ncs_scale)) ncs <- 1
  if (!is.scalar(ncs, 1)) {
    stop("The scale for sum contrasts must be a single positive number")
  }
  if (!is.null(standardized_scale) && standardized_scale != ncs) {
    warning("'formula' is from a standardized object with scale ",
      standardized_scale, " but ncs_scale was specified as ", ncs)
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

