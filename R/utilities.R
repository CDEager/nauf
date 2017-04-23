

as_simplex <- function(x) {
  return(x / sum(x))
}


list_mat_cols <- function(x) {
  return(split(x, c(col(x))))
}


add_quotes <- function(x, collapse = " ") {
  return(paste0("'", x, "'", collapse = collapse))
}


# d is a data.frame
# lvs is a named list of ok levels
# nms specifies the elements in d and lvs to compare
# f is the function to apply to the rows of the logical matrix
match_row <- function(d, lvs, nms = names(lvs), f = all) {
  return(apply(sapply(nms, function(v) {
    d[[v]] %in% lvs[[v]]
  }), 1, f))
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
  m <- which(sapply(lst, function(n) isTRUE(all.equal(n, x))))
  if (length(m)) return(m[1])
  return(0)
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
  if (is.nauf.glm(object)) {
    if (!inherits(object, "glm")) return(gaussian())
    return(get_family(object$family))
  }
  if (is.nauf.lmerMod(object)) return(gaussian())
  if (is.nauf.glmerMod(object)) return(get_family(object@resp$family))
  
  if (inherits(object, "family")) {
    if (!is.null(object$family) && !is.null(object$link)) {
      return(object)
    }
    stop("'family' not recognized")
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
  if (!L) x <- matrix(x)
  if (rm.na) x <- x[!rowna(x), , drop = FALSE]
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

