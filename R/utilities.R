

# reorganized mapply with default for simplify set to FALSE, lower case
# arguments (just don't use with functions where partial matching could
# occur), and MoreArgs -> same.  Just a more intuitive way for me to code.
mlapply <- function(..., same = NULL, simplify = FALSE, use.names = TRUE,
                    fun = NULL) {
  return(mapply(FUN = fun, ..., MoreArgs = same, SIMPLIFY = simplify,
    USE.NAMES = use.names))
}


# args is a list of lists, each of which is passed to do.call(fun)
doapply <- function(fun, args, simplify = FALSE, use.names = TRUE) {
  return(mapply(FUN = do.call, args = args, MoreArgs = list(what = fun),
    SIMPLIFY = simplify, USE.NAMES = use.names))
}


attrapply <- function(x, ..., same = NULL, simplify = FALSE, use.names = TRUE) {
  return(mapply(FUN = structure, .Data = x, ..., MoreArgs = same,
    SIMPLIFY = simplify, USE.NAMES = use.names))
}


# split a matrix (based on split.data.frame but with col option)
msplit <- function(x, f, byrow = TRUE, drop = FALSE, ...) {
  if (byrow) {
    return(lapply(split(x = seq_len(nrow(x)), f = f, drop = drop, ...),
      function(ind) x[ind, , drop = FALSE]))
  }
  return(lapply(split(x = seq_len(ncol(x)), f = f, drop = drop, ...),
    function(ind) x[, ind, drop = FALSE]))
}


# call matrix allowing rownames and colnames to be specified separately
# (useful for mapply), and then optionally call as.data.frame on the result
matframe <- function(data = NA, nrow = 1, ncol = 1, byrow = FALSE,
                     dimnames = NULL, rownames = NULL, colnames = NULL,
                     as_data_frame = TRUE, ...) {
  if (!is.null(rownames) || !is.null(colnames)) {
    if (!is.null(dimnames)) {
      stop("If 'dimnames' is specified, can't specify rownames or colnames.")
    }
    dimnames <- list(rownames, colnames)
  }
  
  ans <- matrix(data = data, nrow = nrow, ncol = ncol, byrow = byrow,
    dimnames = dimnames)
    
  if (as_data_frame) return(as.data.frame(ans, ...))
  
  return(ans)
}


nsplit <- function(x, ord = NULL) {
  if (is.null(ord)) {
    ans <- split(x, names(x))
  } else {
    ans <- split(x[ord], names(x)[ord])
  }
  if (is.list(x) && !is.data.frame(x)) ans <- lapply(ans, unname)
  return(ans)
}


levs_and_contr <- function(fac) {
  return(list(levels = levels(fac), contrasts = contrasts(fac)))
}


# divide non-neg vec by self to make simplex
as_simplex <- function(x) {
  return(x / sum(x))
}


# convert matrix to list with element for each column
list_mat_cols <- function(x) {
  return(split(x, c(col(x))))
}


# find b0 from either terms or formula
has_resp <- function(x) {
  if (!inherits(x, "terms")) x <- stats::terms(x)
  return(attr(x, "response") > 0)
}


# for warnings and errors
add_quotes <- function(x, collapse = " ") {
  return(paste0("'", x, "'", collapse = collapse))
}


# convert effects or group from ranef bar to formula
barform <- function(x, n) {
  return(eval(substitute(~ foo, list(foo = x[[n]]))))
}


# names of variables
varnms <- function(x) {
  if (!inherits(x, "terms")) x <- stats::terms(x)
  x <- attr(x, "factors")
  if (!is.matrix(x)) return(NULL)
  return(rownames(x))
}


# does formula have intercept
has_b0 <- function(x) {
  if (!inherits(x, "terms")) x <- stats::terms(x)
  return(attr(x, "intercept") > 0)
}


# find model frame elements that are neither response nor predictor/group
find_extras <- function(x) {
  cnms <- colnames(x)
  op <- substr(cnms, 1, 1) == "("
  off <- sapply(cnms, function(x) nchar(x) >  7 && substr(x, 1, 7) == "offset(")
  fmat <- attr(attr(x, "terms"), "factors")
  fmat <- rownames(fmat)[rowSums(fmat) == 0]
  fmat <- cnms %in% fmat[-1]
  extras <- op | off | fmat
  extras[1] <- FALSE
  names(extras) <- cnms
  return(extras)
}


# are any cells in the full contingency table of x empty
empty_cells <- function(x) {
  return(any(!xtabs(~ ., x)))
}


# what is the first element in lst which matches x exactly
in_list <- function(x, lst) {
  if (is.character(x)) x <- list(x)
  m <- match(x, lst, nomatch = 0L)
  names(m) <- names(x)
  return(m)
}


# reorder false/true interpretation types
reorder_ft <- function(x) {
  if (in_list(tolower(x), list(c("0", "1"), c("false", "true"), c("f", "t"),
  c("no", "yes"), c("n", "y")))) {
    x <- x[2:1]
  }
  return(x)
}


# is x a single number whose sign is in sgn
is.scalar <- function(x, sgn = c(-1, 0, 1)) {
  return(is.numeric(x) && is.vector(x) && length(x) == 1 && sign(x) %in% sgn)
}


# return family information
#' @importFrom MASS negative.binomial
get_family <- function(object) {
  if (inherits(object, "family")) {
    if (!is.null(object$family) && !is.null(object$link)) {
      return(object)
    }
    stop("'family' not recognized")
  }
  
  if (is.nauf.stanreg(object)) {
    return(object$family)
  }

  if (inherits(object, c("lm", "merMod"))) {
    if (inherits(object, "lm")) {
      if (!inherits(object, "glm")) return(gaussian())
      return(get_family(object$family))
    }
    if (is.nauf.lmerMod(object)) return(gaussian())
    return(get_family(object@resp$family))
  }

  if (is.name(object)) object <- as.character(object)

  if (isTRUE(all.equal(object, MASS::negative.binomial)) ||
  (is.character(object) && object %in% c("negbin", "nb", "negative.binomial",
  "negative binomial"))) {
    return("negbin")
  }

#  if (is.character(object) && object %in% c("ordinal", "ordered")) {
#    return("ordinal")
#  }

  if (is.character(object)) {
    tryfunc <- tryCatch(object <- get(object, mode = "function",
      envir = parent.frame()), error = function(e) e)
    if (inherits(tryfunc, "error")) stop("'family' not recognized")
  }

  if (is.function(object)) {
    object <- object()
    if (!is.null(object$family) && !is.null(object$link)) {
      return(object)
    }
  }

  stop("'family' not recognized")
}


# is x character, logical, or anything with only two unique
# non-NA values but that isn't already an unordered factor
is.charlogbin <- function(x) {
  return(NCOL(x) == 1 && (is.character(x) || is.logical(x) ||
    (!is.uf(x) && nval(x) == 2)))
}


# convert things that aren't factors but should be (according to charlogbin)
# to unordered factors.  x can be a data.frame or a vector
charlogbin_to_uf <- function(x) {
  if (is.data.frame(x)) {
    uf <- which(sapply(x, is.charlogbin))
    for (j in uf) x[[j]] <- factor(x[[j]], ordered = FALSE)
  } else if (is.charlogbin(x)) {
    x <- factor(x, ordered = FALSE)
  }
  return(x)
}


# how many unique values are in x (x a factor, vector, matrix, or data.frame)
nval <- function(x, rm.na = TRUE) {
  if ((L <- length(dim(x))) > 2) stop("'x' must have at most two dims")
  x <- unique(x)
  if (!L) {
    if (rm.na) return(sum(!is.na(x)))
    return(length(x))
  }
  if (rm.na) return(sum(!rowna(x)))
  return(nrow(x))
}


# is x an unordered factor
is.uf <- function(x) {
  return(is.factor(x) && !is.ordered(x))
}


# are any/all of the observations in each row of 'x' NA?
rowna <- function(x, f = any) {
  if (length(dim(x)) != 2) stop("'x' must have two dims")
  return(apply(is.na(x), 1, f))
}


# are any/all of the observations in each column of 'x' NA?
colna <- function(x, f = any) {
  if (length(dim(x)) != 2) stop("'x' must have two dims")
  if (is.data.frame(x)) {
    return(sapply(x, function(n) f(is.na(n))))
  }
  return(apply(is.na(x), 2, f))
}


# get char vec of var names involved in random effects groups
# and make sure all groups have intercepts
check_groups <- function(formula) {
  bars <- lme4::findbars(formula)
  if (!length(bars)) return(NULL)
  g <- lapply(bars, function(x) sort(varnms(barform(x, 3))))
  gu <- unique(g)
  g <- g[sapply(bars, function(x) has_b0(barform(x, 2)))]
  nob0 <- which(!sapply(gu, in_list, g))
  if (length(nob0)) {
    gu <- paste("  ", sapply(gu[nob0], function(x) paste(x, collapse = ":")))
    stop("The following random effects groups do not have random intercepts:\n",
      paste(gu, collapse = "\n"))
  }
  return(sort(unique(unlist(gu))))
}


# glmer vs lmer, etc
is.linear <- function(object){
  return(isTRUE(all.equal(get_family(object), gaussian())))
}


psgn <- function(x, f = mean) {
  return(mean(sign(x) == sign(f(x))))
}


isFALSE <- function(x) {
  return(identical(FALSE, x))
}


condition_on_re <- function(re.form, ReForm, REForm, REform) {
  re.form <- lme4_reFormHack(re.form, ReForm, REForm, REform)
  if (is.null(re.form)) return(TRUE)
  if (is.na(re.form) || isTRUE(all.equal(~ 0, re.form))) return(FALSE)
  stop("Currently, only NULL, NA, and ~0 are supported re.form values.")
}


levasgn <- function(x, add_NEW = FALSE) {
  flist <- get_reTrms(x)$flist
  lvs <- lapply(flist, levels)[attr(flist, "assign")]
  if (add_NEW) lvs <- lapply(lvs, c, "_NEW_")
  return(lvs)
}


last_names <- function(object, nms) {
  d <- length(dim(object))
  if (missing(nms)) {
    if (d) return(dimnames(object)[[d]])
    return(names(object))
  }
  if (d) {
    if (is.list(nms)) {
      dimnames(object)[d] <- nms
    } else {
      dimnames(object)[[d]] <- nms
    }
  } else {
    if (is.list(nms)) nms <- nms[[1]]
    names(object) <- nms
  }
  return(object)
}






###### manipulating classes ######

`drop_class<-` <- function(x, which, value) {
  if (missing(which)) {
    if (isS4(x)) stop("'x' is S4; cannot alter class")
    class(x) <- setdiff(class(x), value)
  } else {
    class(attr(x, which)) <- setdiff(class(attr(x, which)), value)
  }
  return(x)
}


`first_class<-` <- function(x, which, value) {
  if (missing(which)) {
    if (isS4(x)) stop("'x' is S4; cannot alter class")
    class(x) <- union(value, class(x))
  } else {
    class(attr(x, which)) <- union(value, class(attr(x, which)))
  }
  return(x)
}


`last_class<-` <- function(x, which, value) {
  if (missing(which)) {
    if (isS4(x)) stop("'x' is S4; cannot alter class")
    class(x) <- c(setdiff(class(x), value), value)
  } else {
    class(attr(x, which)) <- c(setdiff(class(attr(x, which)), value), value)
  }
  return(x)
}


`add_class<-` <- function(x, which, value) {
  if (!is.list(value) || length(intersect(value[[1]], value[[2]]))) {
    stop("should use first/last_class if 'value' is not a list; ",
      "elements should have no overlap")
  }
  if (missing(which)) {
    if (isS4(x)) stop("'x' is S4; cannot alter class")
    class(x) <- c(value[[1]], setdiff(class(x), do.call(c, value)), value[[2]])
  } else {
    class(attr(x, which)) <- c(value[[1]], setdiff(class(attr(x, which)),
      do.call(c, value)), value[[2]])
  }
  return(x)
}



###### is ######


is.mixed <- function(object) {
  return(length(lme4::findbars(formula(object))) > 0)
}


is.glmod <- function(object) {
  return(is.list(object) && all(c("fr", "X", "reTrms") %in% names(object)))
}


is.nauf.formula <- function(object) {
  return(inherits(object, "nauf.formula"))
}


is.nauf.terms <- function(object) {
  return(inherits(object, "nauf.terms"))
}


is.nauf.frame <- function(object) {
  return(inherits(object, "nauf.frame"))
}


is.nauf.glm <- function(object) {
  return(inherits(object, "nauf.glm"))
}


is.nauf.lmerMod <- function(object) {
  return(inherits(object, "nauf.lmerMod"))
}


is.nauf.glmerMod <- function(object) {
  return(inherits(object, "nauf.glmerMod"))
}


is.nauf.merMod <- function(object) {
  return(inherits(object, c("nauf.lmerMod", "nauf.glmerMod")))
}


is.nauf.model <- function(object) {
  return(inherits(object, c("nauf.glm", "nauf.lmerMod", "nauf.glmerMod",
    "nauf.stanreg")))
}


is.nauf.ref.grid <- function(object) {
  return(inherits(object, "nauf.ref.grid"))
}


is.nauf.mer.anova <- function(object) {
  return(inherits(object, "nauf.mer.anova"))
}


is.nauf.pmm <- function(object) {
  return(inherits(object, "nauf.pmm"))
}


is.summ.nauf.pmm <- function(object) {
  return(inherits(object, "summ.nauf.pmm"))
}


is.nauf.stanreg <- function(object) {
  return(inherits(object, "nauf.stanreg"))
}


is.nauf.stanmer <- function(object) {
  return(inherits(object, "nauf.stanreg") && inherits(object, "lmerMod"))
}


is.nauf.stanfer <- function(object) {
  return(inherits(object, "nauf.stanreg") && !inherits(object, "lmerMod"))
}


is.nauf.mcmc <- function(object) {
  return(inherits(object, "nauf.mcmc"))
}


is.summ.nauf.mcmc <- function(object) {
  return(inherits(object, "summ.nauf.mcmc"))
}


is.nauf.postmm <- function(object) {
  return(inherits(object, "nauf.postmm"))
}


is.summ.nauf.postmm <- function(object) {
  return(inherits(object, "summ.nauf.postmm"))
}

