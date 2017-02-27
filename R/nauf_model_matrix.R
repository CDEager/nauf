
#' Get the model matrix for a model with unordered factor \code{NA} values.
#'
#' \code{nauf_model_matrix} extracts the model matrix from a regression fit
#' with \code{\link{nauf_reg}}, constructs the model matrix specified in
#' a model frame built with \code{\link{nauf_model_frame}}, or constructs
#' a model matrix after calling \code{\link{nauf_model_frame}}, depending
#' on the value of \code{object} and \code{data}.  See 'Details'.
#'
#' Exactly what happens depends on the value of the arguments \code{object}
#' and \code{data}.  The following are organized in order of evaluation:
#'
#' If \code{object} inherits from \code{\link[=nauf-class]{nauf}} and \code{lm}
#' (i.e. it is a regression fit with \code{\link{nauf_reg}}), then the model
#' matrix stored in the regression fit is returned (it is the same as simply
#' calling \code{object$x}).
#'
#' If \code{object} inherits from \code{\link[=nauf-class]{nauf}} and
#' \code{data.frame} (i.e. it is a model frame created by
#' \code{\link{nauf_model_frame}}), then the function begins by calling the
#' default method for \code{\link[stats]{model.matrix}}, and then checking
#' whether the \code{ccna} attribute of the model frame's terms object indicates
#' that any unordered factor interactions require a change in contrasts (see
#' \code{\link{nauf_interaction}}).  For each interaction which requires changes,
#' a copy of the model frame is made and in the copy the factor is re-leveled,
#' dropping levels which are not applicable in the interaction, and assigning
#' contrasts to the new factor with \code{\link{named_contr_sum}}.  Then
#' a model matrix for all terms which contain the relevant unordered factor
#' interaction is generated, and the columns in the original model matrix
#' are replaced with the revised contrasts, appending "ccna_" at the beginning
#' of the column names so that they are easily spotted in the regression output.
#' Finally, all \code{NAs} in the model matrix are set to zero and the matrix
#' is returned.
#'
#' If \code{object} does not meet either of the above specifications, then
#' \code{data} is checked.  If it is a model frame created by
#' \code{\link{nauf_model_frame}}, then a model matrix is created in the way
#' just described.
#'
#' If none of the previous three conditions are met, and \code{object} is a
#' \code{\link[stats]{formula}} and \code{data} is a \code{data.frame}, then
#' \code{\link{nauf_model_frame}} is called, passing along any additional
#' arguments in \code{...}, and a model matrix is built with the resulting
#' model frame.
#'
#' @param object A model frame fit by \code{\link{nauf_model_frame}}, a
#'   regression model fit by \code{\link{nauf_reg}}, a
#'   \code{link[stats]{formula}}, or \code{NULL}.  See 'Details'.
#' @param data A data.frame or \code{NULL}.  See 'Details'.
#' @param ... Additional arguments to pass to \code{\link{nauf_model_frame}}
#'
#' @return A model matrix with unordered factor \code{NAs} set to zero.  It
#'   has an \code{assign} attribute, but no \code{contrasts} attribute.
#'
#' @examples
#' \dontrun{
#' fr <- nauf_model_frame(y ~ f1 * f2 * (f3 + x1), mydata)
#' mod <- nauf_reg(y ~ f1 * f2 * (f3 + x1), mydata)
#'
#' # all of the following return the same
#' x <- nauf_model_matrix(y ~ f1 * f2 * (f3 + x1), mydata)
#' x <- nauf_model_matrix(fr)
#' x <- nauf_model_matrix(data = fr)
#' x <- nauf_model_matrix(mod)
#' }
#'
#' @export
nauf_model_matrix <- function(object = NULL, data = NULL, ...) {
  if (is.nauf(object) & (is.data.frame(object) | inherits(object, "lm"))) {
    return(model.matrix(object))
  } else if (is.nauf(data) & is.data.frame(data)) {
    return(model.matrix(data))
  } else if ("formula" %in% class(object) & is.data.frame(data)) {
    return(model.matrix(nauf_model_frame(formula = object, data = data, ...)))
  }
  stop("either 'object' must be a formula and 'data' a ",
    "data.frame, or one of the two must be a model frame created ",
    "by nauf_model_frame.")
}


#' @describeIn nauf_model_matrix S3 method for class 'nauf'
#'
#' @export
model.matrix.nauf <- function(object = NULL, data = NULL, ...) {
  if (is.nauf(object) & inherits(object, "lm")) {
    return(object$x)
  } else if (is.data.frame(object) & is.nauf(object)) {
    mf <- object
  } else if (is.data.frame(data) & is.nauf(data)) {
    mf <- data
  } else {
    stop("must supply a data.frame or regression model which ",
      "inherits from nauf")
  }

  mf <- nauf_off(mf)
  if (is.null((mt <- attr(mf, "terms")))) {
    stop("'mf' must be a model frame created by nauf_model_frame")
  }
  mt <- nauf_off(mt)
  mta <- attributes(mt)
  if (is.null((ccna <- mta$ccna))) {
    stop("'mf' must be a model frame created by nauf_model_frame")
  }

  mm <- stats::model.matrix(mt, mf)

  if (length(ccna)) {
    # TODO (CDEager): this entire section could be made more efficient
    mmlist <- list()
    cnms <- character()
    asgn <- attr(mm, "assign")
    asgn_ccna <- sort(unique(unlist(lapply(ccna, function(n) n$assign))))
    mmrm <- which(asgn %in% asgn_ccna)
    asgn <- asgn[-mmrm]
    if (length(asgn)) {
      mmlist[[1]] <- mm[, -mmrm, drop = FALSE]
      cnms <- colnames(mmlist[[1]])
    }
    fmat <- mta$factors

    for (i in ccna) {
      uf <- names(i$factors)
      mfi <- mf
      for (u in uf) {
        mfi[, u] <- factor(mfi[, u], ordered = FALSE,
          levels = i$factors[[u]]$levels)
        mfi[, u] <- named_contr_sum(mfi[, u], FALSE)
      }

      imat <- fmat[, i$assign, drop = FALSE]
      main <- rownames(imat)[rowSums(imat) > 0]
      imat <- imat[main, , drop = FALSE]
      form <- paste("~", paste(main, collapse = "+"))
      for(j in 1:ncol(imat)) {
        form <- paste(form, paste(main[imat[, j] > 0], collapse = "*"),
          sep = "+")
      }
      ti <- stats::terms(stats::as.formula(form))
      asgn_ti <- which(colnames(attr(ti, "factors")) %in% colnames(imat))

      mmi <- model.matrix(ti, mfi)
      asgn_mmi <- attr(mmi, "assign")
      keep <- which(asgn_mmi %in% asgn_ti)
      mmi <- mmi[, keep, drop = FALSE]
      asgn_mmi <- asgn_mmi[keep]
      asgn_mmi <- as.numeric(factor(asgn_mmi))
      asgn <- c(asgn, i$assign[asgn_mmi])

      cnmsi <- colnames(mmi)
      sep <- "_"
      new_cnmsi <- paste("ccna", cnmsi, sep = sep)
      # this while-loop shouldn't really ever be needed
      while (any(new_cnmsi %in% cnms)) {
        sep <- paste(sep, "_", sep = "")
        new_cnmsi <- paste("ccna", cnmsi, sep = sep)
      }
      colnames(mmi) <- new_cnmsi
      cnms <- c(cnms, new_cnmsi)

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

