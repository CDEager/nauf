
#' Create named sum contrasts for a factor.
#'
#' \code{named_contr_sum} creates sum contrasts for a factor which are named
#' with the levels of the factor rather than with number (e.g. if a factor
#' \code{f1} has levels \code{A}, \code{B}, and \code{C}, then rather than
#' creating contrast columns \code{f11} and \code{f12}, it creates columns
#' \code{f1A} and \code{f1B}.
#'
#' If \code{x} is a factor, then its levels are sorted (which eliminates
#' \code{NA} as a level if \code{\link[base]{addNA}} has been used, but
#' does not remove any \code{NA} values), and the length of the resulting
#' vector is taken as the number of levels to form sum contrasts with (and
#' the levels are used to name the contrast columns). If \code{x} is a
#' character vector, its unique values are sorted and used as the column names,
#' and its length is used as number of levels.  If \code{x} is a numeric vector
#' with more than one element, then it is sorted and converted to character
#' and treated as if \code{x} had been passed as a character vector.  If
#' \code{x} is a numeric vector with only one element, then its value is used
#' as the number of factor levels, and the columns of the contrast matrix
#' are not given names, resulting in essentially a call to
#' \code{\link[stats]{contr.sum}}.
#'
#' At this point, if there are fewer than two levels, an error is thrown.
#' Otherwise, \code{\link[stats]{contr.sum}} is called to obtain the contrast
#' matrix, and the rows and columns are given the names of the levels; that is,
#' for a factor with \code{K} levels, there are \code{K - 1} columns whose names
#' are the first \code{K - 1} levels, and \code{K} rows whose names are the
#' \code{K} levels.  The names are in alphabetical order excpet for if there are
#' two levels and they are \code{0} and \code{1} or \code{FALSE} and \code{TRUE},
#' in which case they are assigned in the reverse order so that the dummy
#' variable coefficient represents a change toward \code{1} rather than the
#' less intuitive default of representing a change toward \code{0}.
#'
#' Lastly, if \code{x} is a factor and \code{return_factor = TRUE}, then the
#' resulting contrasts are applied to \code{x} and \code{x} is returned rather
#' than the contrast matrix.
#'
#' Using sum contrasts is important for handling not-applicable \code{NA} values
#' in a factor since in the model matrix the \code{NAs} can be set to zero for
#' all of the dummy variables, averaging over the effect of the factor rather
#' than causing the \code{NAs} to be coded as the base level as would occur
#' with treatment contrasts (in sociolinguistics this is often referred to as
#' factor slashing).  In the \code{nauf} package, all unordered factors
#' are automatically coded with named sum contrasts (as are any variables not
#' coded as unordered factors, but which have only two unique values; e.g.
#' logicals, variables coded as integer \code{0 / 1}, etc.).  See
#' \code{\link{nauf_model_frame}} and \code{\link{nauf_interaction}} for details.
#'
#' @param x Either a factor or a numeric or character vector.
#' @param return_factor A logical indicating whether a factor (\code{TRUE}) or
#'   contrast matrix (\code{FALSE}, the default) should be returned.
#'   If \code{x} is not a factor, then \code{return_factor} is ignored.
#'
#' @return If \code{x} is a factor and \code{return_factor = TRUE}, then
#'   an unordered factor with named sum contrasts is returned.  Otherwise,
#'   a matrix of named sum contrasts is returned.
#'
#' @seealso \code{\link{nauf_interaction}} for interactions involving unordered
#'   factors where at least one has \code{NA} values,
#'   \code{\link{nauf_model_frame}} for automatic assignment of unordered factor
#'   contrasts, and \code{\link{nauf_model_matrix}} for the treatment of
#'   \code{NAs} in regressions.
#'
#' @examples
#' f <- factor(rep(c("a", "b", "c", NA), 2))
#' f <- addNA(f)
#' levels(f)  # NA listed as factor
#' contrasts(f)  # NA included in contrast matrix
#' named_contr_sum(f)  # contrast matrix with named columns and NA excluded
#' named_contr_sum(f, TRUE)  # return factor with same contrasts
#' 
#' named_contr_sum(letters[1:5])  # character method
#' named_contr_sum(5)  # same as contr.sum
#' named_contr_sum(c(1, 3, 6))  # converted to character
#' named_contr_sum(c(0, 1))  # treated differently (1 before 0)
#'
#' \dontrun{
#' named_contr_sum(1)  # must be >= 2
#' named_contr_sum(rep(c("a", NA), 3))  # only one unique non-NA
#' }
#'
#' @importMethodsFrom stats contrasts<- contrasts
#'
#' @export
named_contr_sum <- function(x, return_factor = FALSE) {
  if (!is.factor(x) & !is.numeric(x) & !is.character(x)) {
    stop("'x' must be a factor or a numeric or character vector")
  }
  
  if (is.factor(x)) {
    lvs <- sort(levels(x))
    n <- length(lvs)
  } else if (is.character(x)) {
    lvs <- sort(unique(x))
    n <- length(lvs)
  } else if (is.numeric(x)) {
    if (length(x) > 1) {
      lvs <- sort(as.character(unique(x)))
      n <- length(lvs)
    } else {
      n <- x
      lvs <- NULL
    }
  } else {
    stop("'x' must be a factor or a numeric or character vector")
  }
  
  if (n < 2) {
    stop("contrasts cannot be applied to factors with fewer than two levels")
  }
  
  if (n == 2 & !is.null(lvs)) {
    if (all(lvs == c("FALSE", "TRUE")) | all(lvs == c("0", "1"))) {
      lvs <- lvs[2:1]
    }
  }
  
  contr <- contr.sum(n)
  if (!is.null(lvs)) {
    rownames(contr) <- lvs
    colnames(contr) <- lvs[-n]
    if (is.factor(x) & return_factor) {
      x <- factor(x, ordered = FALSE, levels = lvs)
      contrasts(x) <- contr
      return(x)
    }
  }
  
  return(contr)
}



#' Unordered factor interactions involving \code{NAs}.
#'
#' The function determines whether unordered factor levels need to be dropped
#' in the computation of an interaction term due to \code{NA} values, treating
#' characters and logicals as unordered factors.
#'
#' If there are not at least two unordered factors in the specified columns,
#' an error is thrown. Otherwise, the unordered factors involved in the
#' interaction are subsetted.  If none of these has \code{NA} values, then
#' function stops and returns a list of the levels for each factor.  If any
#' of the unordered factors has \code{NA} values, then a table of all their
#' unique combinations in \code{x} is made, excluding any combination which
#' contains \code{NAs}.  Then the unique combinations are examined two columns
#' at a time, removing redundant rows, and then dropping levels which
#' only occur once within-column.  Then the two columns are collapsed as a new
#' column which is compared with the next column, continuing the process until
#' there is only one column.  If there is only one unique combination at any
#' point, then an error is thrown.  Otherwise, the levels for each factor
#' which should be used in the interaction are returned.  Note that only the
#' highest order interaction of the unordered factors is considered.  If
#' unordered factors \code{A}, \code{B}, and \code{C} are involved in an
#' interaction, it may be that \code{A} needs levels dropped in the third-order
#' interaction, but not in the interaction \code{A:B}, for example.  The
#' function would need to be called separately with \code{cols = c("A", "B")}
#' to determine the levels for the second-order interaction.
#'
#' This function is implemented in \code{\link{nauf_model_frame}}, and is
#' necessary in cases where \code{NA} values are collinear with non-\code{NA}
#' values in other unordered factors.  For example, consider a hypothetical
#' study where participants speak one of three different dialects of Spanish,
#' which we will denote \code{A}, \code{B}, and \code{C}.  Suppose that in
#' dialects \code{A} and \code{B}, some speakers are bilingual in the (same)
#' second language while others are monolingual, but in dialect \code{C}, all
#' speakers are monolingual.  In this case, we can code two factors.  The
#' factor \code{dialect} has levels \code{A}, \code{B}, and \code{C}.  The
#' factor \code{bilingual} has levels \code{TRUE} and \code{FALSE}, with all
#' observations pertaining to \code{dialect C} coded as \code{NA}, since the
#' the factor is not contrastive within the group (i.e. being monolingual in
#' \code{dialect C} is different from being monolingual in the other two).
#' Using \code{\link{named_contr_sum}}, we assign the following contrasts (in
#' the \code{bilingual} table, the \code{NA} row is shown here to be explicit;
#' it is not a part of the contrast matrix returned by the function because
#' the \code{\link[stats]{contrasts}} requires the number of rows to be one
#' more than the number of columns; the \code{NAs} are set to zero by
#' \code{\link{nauf_model_matrix}}):
#'
#' Contrasts for main effect of \code{dialect}
#' \tabular{rr}{
#'     \tab A \tab B\cr
#'   A \tab  1 \tab  0\cr
#'   B \tab  0 \tab  1\cr
#'   C \tab -1 \tab -1
#' }
#' 
#' Contrasts for main effect of \code{bilingual}
#' \tabular{r}{
#'         \tab TRUE\cr
#'   TRUE  \tab  1\cr
#'   FALSE \tab -1\cr
#'   NA    \tab  0
#' }
#' 
#' This setup allows the regression coefficient \code{bilingualTRUE} to
#' only apply when either \code{dialectA} or \code{dialectB} is \code{1},
#' and to always be multiplied by \code{0} when \code{dialectA} and
#' \code{dialectB} are \code{-1} (i.e. for \code{dialect C}).  In the regression
#' output, the predicted grand mean for each dialect is obtained by applying
#' the contrast coding to the \code{dialect} coefficients.  That is, we have:
#'
#' \code{
#'   mean(A) = (Intercept) + dialectA  
#'   mean(B) = (Intercept) + dialectB 
#'   mean(C) = (Intercept) - dialectA - dialectB
#' }
#'
#' For dialects \code{A} and \code{B}, this estimate averages over the effect
#' of \code{bilingual}, and the estimates for the sub-groups can be obtained by:
#'
#' \code{
#'   mean(A:TRUE)  = mean(A) + bilingualTRUE
#'   mean(A:FALSE) = mean(A) - bilingualTRUE
#'   mean(B:TRUE)  = mean(B) + bilingualTRUE
#'   mean(B:FALSE) = mean(B) - bilingualTRUE
#' }
#'
#' Provided that any covariates have been put on unit scale (see 
#' \code{\link[base]{scale}}) and orthogonal polynomial contrasts have been set
#' for any ordered factors (see \code{\link[stats]{contr.poly}}), both of which
#' are good ideas for most regressions anyway, then these types of contrasts
#' also result in the intercept being the corrected mean, and the output is
#' easier to interpret.  An additional (and more important) advantage to
#' implementing \code{NA} values, however, comes when factors such as the ones
#' in this example interact.  Assume now that we are interested in the
#' interaction \code{dialect * bilingual}, rather than just the main effects.
#' Coding the observations from dialect \code{C} as \code{bilingual = FALSE}
#' yields a rank-deficient matrix (i.e. there are collinear columns).  If we
#' simply multiply the contrasts for the main effects outlined above in the
#' interaction (which is what the default for \code{link[stats]{model.matrix}}
#' does), the result is not collinear; hoever, it is still undesirable:
#'
#' Undesirable contrasts for \code{dialect:bilingual} interaction
#' \tabular{llrr}{
#'    \tab  \tab dialectA:bilingualTRUE \tab dialectB:bilingualTRUE\cr
#'   A \tab TRUE  \tab  1 \tab  0\cr
#'   A \tab FALSE \tab -1 \tab  0\cr
#'   B \tab TRUE  \tab  0 \tab  1\cr
#'   B \tab FALSE \tab  0 \tab -1\cr
#'   C \tab NA    \tab  0 \tab  0
#' }
#' 
#'
#' These contrasts are undesirable because the same interaction term could be
#' expressed with one column with some releveling, and in cases more complicated
#' than this one, there could also be collinearity.  So, instead, we should use
#' the following contrasts:
#'
#' \tabular{llr}{
#'    \tab  \tab dialectA:bilingualTRUE\cr
#'   A \tab TRUE  \tab  1\cr
#'   A \tab FALSE \tab -1\cr
#'   B \tab TRUE  \tab -1\cr
#'   B \tab FALSE \tab  1\cr
#'   C \tab NA    \tab  0
#' }
#' 
#' Determining these contrasts manually can be tedious (especially as the number
#' of factors and levels grows), hence \code{nauf_interaction}.
#' For this example \code{nauf_interaction} determines that \code{C} is
#' redundant in the interaction term, and so it drops it from the levels
#' of \code{dialect}, while keeping \code{bilingual} the same.  In this case,
#' \code{nauf_interaction} would return a list with an element \code{levels}
#' which has sub-elements \code{dialect} and \code{bilingual}, each of which is
#' a character vector (\code{c("A", "B")} and \code{c("TRUE", "FALSE")},
#' respectively), along with a second element \code{changed = TRUE} to indicate
#' that at least one level was dropped from at least one factor.  See
#' \code{\link{nauf_model_matrix}} for implementation of interaction contrasts
#' which don't match main effect contrasts, and how this effects interpretation
#' of the regression output.  See the examples section for code
#' implementing this example.
#'
#' @param x A data.frame.
#' @param cols A vector specifying columns in \code{x} involved in an
#'   interaction.  At least two must be unordered factors.  Defaults to
#'   all columns in \code{x}.
#'
#' @return A list with two elements.  The first, \code{levels}, is a named list
#'   with one entry for each unordered factor in \code{cols}.  Each entry is a
#'   character vector of the factor levels which should be coded for in the full
#'   interaction term (with sum contrasts).  The second element, \code{changed},
#'   is a logical indicating whether or not \code{levels} is the same as the
#'   levels in the main effects.
#'
#' @seealso \code{\link{named_sum_contr}} for the contrasts that \code{nauf}
#'   assigns to unordered factors,
#'   \code{\link{nauf_model_frame}} for automatic assignment of unordered factor
#'   contrasts, and \code{\link{nauf_model_matrix}} for the treatment of
#'   \code{NAs} in regressions.
#'
#' @examples
#' dat <- data.frame(
#'   dialect = c("A", "A", "B", "B", "C"),
#'   bilingual = c(TRUE, FALSE, TRUE, FALSE, NA))
#' nauf_interaction(dat)  # drops dialect C
#' 
#' dat <- as.data.frame(lapply(dat, function(n) rep(n, 10)))
#' dat$x <- rnorm(nrow(dat))
#' dat$f <- rep(c("U", "V"), 25)
#' nauf_interaction(dat, c("dialect", "bilingual")) # same as initial example
#' nauf_interaction(dat)  # ignores x, checks dialect:bilingual:f
#' 
#' nauf_interaction(dat, c("dialect", "f"))
#' # sees that neither has NAs and returns the levels of the factors
#' 
#' \dontrun{
#' nauf_interaction(dat, "dialect")  # error; only one column
#' nauf_interaction(dat, c("dialect", "x"))  # error; only one unordered factor
#' }
#'
#' @export
nauf_interaction <- function(x, cols = colnames(x)) {
  if (!is.data.frame(x)) stop("'x' must be a data.frame'")
  if (length(cols) < 2) stop("must supply more than one factor")
  if (!is.character(cols)) {
    if (is.numeric(cols)) {
      cols <- colnames(x)[cols]
    } else {
      stop("'cols' must be a character, numeric, or logical vector")
    }
  }
  
  x <- x[, cols]
  makefac <- unlist(lapply(x, function(n) any(c("character", "logical") %in%
    class(n))))
  for (j in cols[makefac]) {
    x[, j] <- factor(x[, j])
  }
  uf <- unlist(lapply(x, function(n) is.factor(n) & !is.ordered(n)))
  if (sum(uf) < 2) {
    stop("interaction must involve at least two unordered factors")
  }
  
  x <- x[, uf]
  cols <- colnames(x)
  mainlvs <- lapply(x, function(n) sort(levels(n)))
  if (!any(unlist(lapply(x, anyNA)))) {
    return(list(levels = mainlvs, changed = FALSE))
  }
  
  # unique character x with no NAs
  ux <- as.matrix(unique(x))
  ux <- ux[apply(!is.na(ux), 1, all), , drop = FALSE]
  rownames(ux) <- 1:nrow(ux)
  
  # reduce ux
  rx <- ux
  while (ncol(rx) > 1) {
    tab <- xtabs(~ rx[, 1] + rx[, 2])
    d <- dimnames(tab)
    d1 <- d[[1]][rowSums(tab) > 1]
    d2 <- d[[2]][colSums(tab) > 1]
    if (!length(d1) || !length(d2)) {
      stop("collinearity in the interaction term ",
        paste(colnames(ux), collapse = ":"))
    }
    rows <- which(rx[, 1] %in% d1 & rx[, 2] %in% d2)
    rx[, 2] <- paste(rx[, 1], rx[, 2], sep = ".")
    rx <- rx[rows, -1, drop = FALSE]
  }
  
  ux <- data.frame(ux[rownames(rx), , drop = FALSE], stringsAsFactors = FALSE)
  lvs <- lapply(ux, function(n) sort(unique(n)))
  changed <- !is.logical(all.equal(mainlvs, lvs))
  
  return(list(levels = lvs, changed = changed))
}


# 'main effect factor contrasts'
# get factor contrasts applied to main effects terms from a model fit with nauf
get_mefc <- function(object) {
  if (inherits(object, "lm")) {
    a <- attributes(object$terms)
    if (length(a$mefc)) {
      uf <- a$mefc[unlist(lapply(a$mefc, function(n) !n$ordered))]
      of <- a$mefc[unlist(lapply(a$mefc, function(n) n$ordered))]
      uf <- if (length(uf)) lapply(uf, function(n) n$contrasts)
      of <- if (length(of)) lapply(of, function(n) n$contrasts)
      for (u in names(uf)) {
        if (a$hasna[u]) {
          uf[[u]] <- rbind(uf[[u]], 0)
          rownames(uf[[u]])[nrow(uf[[u]])] <- NA
        }
      }
      return(list(unordered = uf, ordered = of))
    }
  }
  return(NULL)
}


# 'contrast changes due to NA values'
# get factor contrast changes due to NA values from a model fit with nauf
get_ccna <- function(object) {
  if (inherits(object, "lm")) {
    a <- attributes(object$terms)
    if (length(a$ccna)) {
      ccna <- lapply(a$ccna, function(n) lapply(n$levels, named_contr_sum))
      uf <- get_mefc(object)$unordered
      for (j in 1:length(ccna)) {
        for (u in names(ccna[[j]])) {
          lvs <- rownames(uf[[u]])
          rn <- rownames(ccna[[j]][[u]])
          lvs <- lvs[!(lvs %in% rn)]
          if (length(lvs)) {
            z <- matrix(0, length(lvs), ncol(ccna[[j]][[u]]),
              dimnames = list(lvs, colnames(ccna[[j]][[u]])))
            ccna[[j]][[u]] <- rbind(ccna[[j]][[u]], z)
          }
        }
      }
      return(ccna)
    }
  }
  return(NULL)
}

#' Contrasts for a model fit with \code{\link{nauf_reg}}.
#'
#' \code{nauf_contrasts} lists all factor contrasts in a model fit with
#' \code{\link{nauf_reg}}, including changes in contrasts for interactions
#' involving \code{NAs}.
#'
#' @param object A model fit using \code{\link{nauf_reg}}.
#'
#' @return A list with two elements: \code{mefc} (main effect factor contrasts)
#'   and \code{ccna} (contrast changes due to \code{NAs}).
#' 
#'   \code{mefc} is a list with two elements \code{unordered} and \code{ordered},
#'   each of which is a named list with an entry for each unordered or ordered
#'   factor in the regression which contains the contrasts applied to the main
#'   effect.
#'
#'   \code{ccna} is a named list with an element for each combination of
#'   unordered factors which required a change in factor contrasts from the
#'   main effects, with an entry for each of the factors containing the
#'   contrasts applied to any term involving only those unordered factors
#'   (but possibly also ordered factors and covariates).
#'
#' @examples
#' \dontrun{
#' m <- nauf_reg(y ~ f1 * f2 * (f3 + x1 + x2), data = mydata)
#' nauf_contrasts(m)
#' }
#' @export
nauf_contrasts <- function(object) {
  return(list(mefc = mefc(object), ccna = ccna(object)))
}
