

## this file contains internal functions taken from afex v0.17-8
## which are altered as little as possible and preceded by 'afex_'


afex_get_mixed_warnings <- function(x) {
  full_model_name <- names(x)[[2]]
  ntry <- function(x) tryCatch(x, error = function(e) NULL)
  if (is.list(x$full)) {
    warnings1 <- c(full = lapply(x[[2]], function(y) y@optinfo$warnings),
      lapply(x[[3]], function(y) y@optinfo$warnings))
    warnings2 <- c(full = lapply(x[[2]], function(y) ntry(y@optinfo$conv$lme4$messages)),
      lapply(x[[3]], function(y) ntry(y@optinfo$conv$lme4$messages)))
  } else {
    warnings1 <- c(full = list(x[[full_model_name]]@optinfo$warnings),
      lapply(x[[3]], function(y) y@optinfo$warnings))
    warnings2 <- c(full = list(ntry(x[[full_model_name]]@optinfo$conv$lme4$messages)),
      lapply(x[[3]], function(y) ntry(y@optinfo$conv$lme4$messages)))
  }
  warnings <- mapply(function(x, y) c(unlist(x), y), warnings1, warnings2,
    SIMPLIFY = FALSE)
  warn <- vapply(warnings, function(y) !length(y) == 0, NA)
  for (i in names(warn)[warn]) {
    warning("lme4 reported (at least) the following warnings for '",
      i, "':\n  * ", paste(warnings[[i]], collapse = "\n  * "),
      call. = FALSE)
  }
}

