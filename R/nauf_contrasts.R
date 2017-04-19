

#' @export
nauf_contrasts <- function(object) {
  if (is.nauf.model(object)) {
    object <- model.frame(object)
  }
  if (is.nauf.frame(object)) {
    object <- attr(object, "terms")
  }
  if (is.nauf.terms(object)) {
    info <- attr(object, "nauf.info")
  } else {
    stop("Must supply a nauf.terms, nauf.frame, or nauf model")
  }
  
  hasna <- info$hasna
  uf <- info$uf
  of <- info$of
  ncs <- info$ncs_scale
  contr <- list()
  
  for (v in names(uf)) {
    v.contr <- lapply(uf[[v]], standardize::named_contr_sum, scale = ncs)
    if (length(uf[[v]]) > 1) {
      for (cj in 2:length(uf[[v]])) {
        colnames(v.contr[[cj]]) <- paste0(".c", cj, ".",
          colnames(v.contr[[cj]]))
        z <- matrix(0, nrow(v.contr[[1]]), ncol(v.contr[[cj]]),
          dimnames = list(rownames(v.contr[[1]]), colnames(v.contr[[cj]])))
        z[rownames(v.contr[[cj]]), ] <- v.contr[[cj]]
        v.contr[[cj]] <- z
      }
    }
    v.contr <- do.call(cbind, args = v.contr)
    if (hasna[v]) {
      v.contr <- rbind(v.contr, 0)
      rownames(v.contr)[nrow(v.contr)] <- NA
    }
    contr[[v]] <- v.contr
  }
  
  if (length(of)) contr <- c(contr, lapply(of, function(x) x$contrasts))
  
  if (!length(contr)) return(NULL)
  
  return(contr)
}

