

# means of the variables which go into a numeric term
# for use in pmmeans grid
get_num_trm_means <- function(term, data) {
  formula <- stats::formula(paste("~", term))
  vars <- all.vars(formula)
  vmiss <- vars[!(vars %in% colnames(data))]
  if (length(vmiss)) {
    vmiss <- paste(paste(" ", vmiss), collapse = "\n")
    stop("Unable to evalutate the term ", term,
      " because the follwing are not in 'data':\n", vmiss)
  }
  vmat <- sapply(data[, vars, drop = FALSE], function(x) length(dim(x)) > 1)
  vmat <- names(vmat)[vmat]
  if (length(vmat)) {
    vmat <- paste(paste(" ", vmat), collapse = "\n")
    stop("Unable to evalutate the term ", term,
      " because the follwing are matrices in 'data':\n", vmiss)
  }
  return(lapply(data[, vars, drop = FALSE], mean))
}


# which regression function should be called
get_regfunc <- function(formula, family) {
  if (is.mixed(formula)) {
    if (is.character(family)) {
      return(quote(nauf::nauf_glmer.nb))
    } else if (family$family == "gaussian" && family$link = "identity") {
      return(quote(nauf::nauf_lmer))
    } else {
      return(quote(nauf::nauf_glmer))
    }
  } else {
    if (is.character(family)) {
      return(quote(nauf::nauf_glm.nb))
    } else if (family$family == "gaussian" && family$link = "identity") {
      return(quote(nauf::nauf_lm))
    } else {
      return(quote(nauf::nauf_glm))
    }
  }
}


# what is the first element in lst which matches x exactly
in_list <- function(x, lst) {
  m <- which(sapply(lst, function(n) isTRUE(all.equal(n, x))))
  if (length(m)) return(m[1])
  return(0)
}


# reorder false/true interpretation types
reorder_ft <- function(x) {
  ft <- list(
    c("0", "1"),
    c("false", "true"),
    c("f", "t"),
    c("no", "yes"),
    c("n", "y"))
  
  if (in_list(tolower(x), ft)) {
    x <- x[2:1]
  }
  
  return(x)
}


# is x a single number whose sign is in sgn
is.scalar <- function(x, sgn = c(-1, 0, 1)) {
  return(is.numeric(x) && is.vector(x) && length(x) == 1 && sign(x) %in% sgn)
}


# is x character, logical, or anything with only two unique
# non-NA values but that isn't already an unordered factor
charlogbin <- function(x) {
  return((is.null(dim(x)) && any(c("character", "logical") %in% class(x))) ||
    (!is.uf(x) && length(sort(unique(x))) == 2))
}


# convert things that aren't factors but should be (according to charlogbin)
# to unordered factors.  x can be a data.frame or a vector
charlogbin_to_uf <- function(x) {
  if (is.data.frame(x)) {
    for (j in 1:ncol(x)) {
      x[[j]] <- charlogbin_to_uf(x[[j]])
    }
  } else if (charlogbin(x)) {
    x <- factor(x, ordered = FALSE)
  }
  return(x)
}


# is x an unordered factor
is.uf <- function(x) {
  return(is.factor(x) && !is.ordered(x))
}


# are any cells in the full contingency table of x empty
empty_cells <- function(x) {
  return(any(!xtabs(~ ., x)))
}


# get the names of variables that are not involved in any terms
non_term_names <- function(x) {
  if (inherits(x, "terms")) {
    mt <- x
  } else if (is.data.frame(x)) {
    mt <- attr(x, "terms")
  } else {  # try to get terms; may fail
    mt <- terms(x)
  }
  fmat <- attr(mt, "factors")
  if (!is.matrix(fmat)) return(NULL)
  return(rownames(fmat)[rowSums(fmat) == 0])
}


# get the names of the non-response and non-offset vars from a terms obj
predictors <- function(x) {
  if (inherits(x, "terms")) {
    mt <- x
  } else if (is.data.frame(x)) {
    mt <- attr(x, "terms")
  } else {  # try to get terms; may fail
    mt <- terms(x)
  }
  fmat <- attr(mt, "factors")
  if (!is.matrix(fmat)) return(NULL)
  return(rownames(fmat)[rowSums(fmat) > 0])
}


# are any/all of the observations in each row of 'x' NA?
rowna <- function(x, f = any) {
  if (is.null(dim(x)) || length(dim(x)) != 2) stop("'x' must have two dims")
  return(apply(is.na(x), 1, f))
}


# are any/all of the observations in each column of 'x' NA?
colna <- function(x, f = any) {
  if (is.null(dim(x)) || length(dim(x)) != 2) stop("'x' must have two dims")
  if (is.data.frame(x)) {
    return(sapply(x, function(n) f(is.na(n))))
  }
  return(apply(is.na(x), 2, f))
}


# bar is a single element from lme4::findbars
get_bar_info <- function(mf, bar) {
  fg <- stats::terms(stats::as.formula(substitute(~foo, list(foo = bar[[3]]))))
  d <- all.vars(fg)
  g.fac <- rownames(attr(attr(fg, "terms"), "factors"))
  
  fs <- stats::terms(stats::as.formula(substitute(~foo, list(foo = bar[[2]]))))
  b0 <- attr(attr(fs, "terms"), "itercept")
  s.fac <- attr(attr(fs, "terms"), "factors")
  if (!is.matrix(s.fac)) s.fac <- NULL  # intercept only
  
  lvs <- levels(interaction(mf[, g.fac, drop = FALSE], drop = TRUE))
  
  return(
    list(group = list(terms = fg, datanames = d, factors = g.fac, levels = lvs),
    effects = list(terms = fs, intercept = b0, factors = s.fac)))
}


# object is result of lme4::findbars or is a formula with re terms
group_factors <- function(object) {
  if (inherits(object, "formula")) object <- lme4::findbars(object)
  g <- lapply(object, get_bar_info)
  g <- sort(unique(unlist(lapply(g, function(x) x$factors))))
  return(unname(g))
}


# return family information
#' @importFrom stringr str_locate_all
get_family <- function(family) {
  if (inherits(family, "family")) {
    if (!is.null(family$family)) {
      return(family)
    }
    stop("'family' not recognized")
  }
  
  if (isTRUE(all.equal(family, MASS::negative.binomial)) ||
  (is.character(family) && family %in% c("negbin", "nb", "negative.binomial",
  "negative binomial"))) {
    return("negbin")
  }
  
  if (is.character(family)) {
    tryfunc <- tryCatch(family <- get(family, mode = "function",
      envir = parent.frame()), error = function(e) e)
    if (inherits(tryfunc, "error")) stop("'family' not recognized")
  }
  
  if (is.function(family)) {
    family <- family()
    if (is.null(family$family)) {
      stop("'family' not recognized")
    }
    return(family)
  }
  
  stop("'family' not recognized")
}

