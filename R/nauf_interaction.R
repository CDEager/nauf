

nauf_interaction <- function(x, cols = colnames(x)) {
  if (!is.data.frame(x)) stop("'x' must be a data.frame'")
  if (!is.character(cols)) {
    if (is.numeric(cols) || is.logical(cols)) {
      cols <- colnames(x)[cols]
    } else {
      stop("'cols' must be a character, numeric, or logical vector")
    }
  }
  cols <- unique(cols)
  nm <- paste(cols, collapse = ":")

  # coerce char/logical/binary to unordered factor and record levels
  x <- droplevels(charlogbin_to_uf(x[, cols, drop = FALSE]))
  uf <- sapply(x, is.uf)
  if (sum(uf) < 2) {
    stop("interaction ", nm, " does not involve at least two unordered factors")
  }
  x <- x[, uf]
  mlvs <- lapply(x, function(n) reorder_ft(sort(levels(n))))

  # remove any unique combination which involves NAs
  x <- unique(x)
  x <- droplevels(x[!rowna(x), , drop = FALSE])
  if (!nrow(x)) stop("No unique applicable combinations in ", nm)
  
  # remove levels which are redundant
  # e.g. f1 in [a, b, c], f2 in [d, e, f, g]
  #    when f1 = a, f2 in [d, e, f]
  #    when f1 = b, f2 in    [e, f, g]
  #    when f1 = c, f2 = NA
  #    then at this point [c] has been dropped form f1,
  #    but we still need to drop [d, g] from f2
  if (empty_cells(x)) {
    torm <- vector("list", ncol(x))
    
    for (j in 1:length(torm)) {
      for (i in levels(x[, j])) {
        check <- droplevels(x[x[, j] == i, -j, drop = FALSE])
        if (any(sapply(check, nlevels) == 1)) {
          torm[[j]] <- c(torm[[j]], i)
        }
      }
    }
    
    for (j in 1:length(torm)) {
      x[x[, j] %in% torm[[j]], j] <- NA
    }
    
    x <- unique(droplevels(x[!rowna(x), , drop = FALSE]))
    if (!nrow(x)) stop("No unique applicable combinations in ", nm)
  }
  
  ilvs <- lapply(x, function(n) reorder_ft(sort(levels(n))))
  changed <- !isTRUE(all.equal(mlvs, ilvs))
  if (any(unlist(lapply(x, nlevels)) < 2)) {
    stop("At least one factor in ", nm,
      " has only one level when NAs are removed")
  }
  
  # there can still be collinearity, but in this case it isn't structural
  # so warning rather than error (i.e. a column will be dropped from
  # the model matrix if the interaction is included in a regression)
  if (empty_cells(x)) warning("Collinearity in the interaction ", nm)

  return(list(levels = ilvs, changed = changed))
}

