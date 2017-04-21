context("nauf_model.frame")

set.seed(1)
dat <- rbind(
  expand.grid(f1 = c("a", "b"), f2 = c(TRUE, FALSE),
    f3 = NA, f4 = c(2, 3), stringsAsFactors = FALSE),
  expand.grid(f1 = c("c"), f2 = c(TRUE, FALSE),
    f3 = c("u", "v", "w"), f4 = c(2, 3), stringsAsFactors = FALSE),
  expand.grid(f1 = c("d"), f2 = NA,
    f3 = c("u", "v", "w"), f4 = c(2, 3), stringsAsFactors = FALSE))
dat <- dat[order(dat$f1, dat$f2, dat$f3, dat$f4), ]
rownames(dat) <- NULL
dat[25:26, 4] <- NA
dat <- as.data.frame(lapply(dat, function(n) rep(n, 10)))
dat$o <- factor(rep(1:3, length.out = nrow(dat)), ordered = TRUE)
contrasts(dat$o) <- contr.poly(3)
dat$x <- stats::rnorm(nrow(dat))
dat$z <- stats::rnorm(nrow(dat))
dat$y <- stats::rnorm(nrow(dat))
dat$w <- w <- stats::runif(nrow(dat))
dat$os <- os <- stats::rnorm(nrow(dat))

form1 <- y ~ f2 * f4 * (f1 + f3) + f1 * f3 * (x + o + poly(z, 3))

mf1 <- stats::model.frame(form1, dat, subset = dat$x < 2, offset = os,
  weights = w, na.action = na.pass)
  
for (j in c("f1", "f2", "f3", "f4")) {
  mf1[[j]] <- standardize::named_contr_sum(mf1[[j]], 1, FALSE)
}
  
ni1 <- list(
  resp = "y",
  groups = list(),
  uf = list(
    f2 = list(c("TRUE", "FALSE")),
    f4 = list(c("2", "3")),
    f1 = list(c("a", "b", "c", "d"),
      c("a", "b", "c"),
      c("c", "d")),
    f3 = list(c("u", "v", "w"))),
  of = list(o = list(
    levels = c("1", "2", "3"),
    contrasts = contr.poly(3))),
  num = list(x = mean(subset(dat, x < 2)$x)),
  mat = list("poly(z, 3)" = colMeans(poly(dat$z, 3)[dat$x < 2, ])),
  extras = c("(offset)", "(weights)"),
  cc = list(list(
    "f1:f2" = list(factors = c(2, 1), assign = 9),
    "f1:f3" = list(factors = c(3, 1), assign = c(13, 22, 23, 24)),
    "f1:f2:f4" = list(factors = c(2, 1, 1), assign = 20))),
  hasna = c(F, T, T, F, T, F, F, F, F, F),
  ncs_scale = 1
)
names(ni1$cc[[1]][[1]][[1]]) <- c("f1", "f2")
names(ni1$cc[[1]][[2]][[1]]) <- c("f1", "f3")
names(ni1$cc[[1]][[3]][[1]]) <- c("f1", "f2", "f4")
names(ni1$cc[[1]][[1]][[2]]) <- "f2:f1"
names(ni1$cc[[1]][[2]][[2]]) <- c("f1:f3", "f1:f3:x", "f1:f3:o",
  "f1:f3:poly(z, 3)")
names(ni1$cc[[1]][[3]][[2]]) <- "f2:f4:f1"
names(ni1$hasna) <- c("y", "f2", "f4", "f1", "f3", "x", "o",
  "poly(z, 3)", "(offset)", "(weights)")

class(mf1) <- c("data.frame", "nauf.frame")
mt <- attr(mf1, "terms")
class(mt) <- c("nauf.terms", "terms", "formula")
attr(mt, "nauf.info") <- ni1
attr(mt, "dataClasses")[names(ni1$uf)] <- "factor"
attr(mf1, "terms") <- mt
attr(mf1, "formula") <- form1
class(attr(mf1, "formula")) <- c("nauf.formula", "formula")

allcontr <- list(
  f2 = rbind(contr.sum(2), 0),
  f4 = rbind(contr.sum(2), 0),
  f1 = cbind(contr.sum(4), rbind(contr.sum(3), 0), rbind(0, 0, contr.sum(2))),
  f3 = rbind(contr.sum(3), 0),
  o = contr.poly(3))
dimnames(allcontr[[1]]) <- list(c("TRUE", "FALSE", NA), "TRUE")
dimnames(allcontr[[2]]) <- list(c("2", "3", NA), "2")
dimnames(allcontr[[3]]) <- list(c("a", "b", "c", "d"),
  c("a", "b", "c", ".c2.a", ".c2.b", ".c3.c"))
dimnames(allcontr[[4]]) <- list(c("u", "v", "w", NA), c("u", "v"))
rownames(allcontr[[5]]) <- c("1", "2", "3")

test_that("basic functionality works", {
  expect_equal((nmf1 <- nauf_model.frame(form1, dat, subset = dat$x < 2, offset = os,
    weights = w)), mf1)
  expect_equal(nauf_contrasts(nmf1, TRUE), allcontr)
  expect_equal(nauf_contrasts(nmf1), allcontr[-5])
})

mf2 <- mf1
attr(attr(mf1, "terms"), "nauf.info")$ncs_scale <- 0.5
for (j in c("f1", "f2", "f3", "f4")) {
  mf1[[j]] <- standardize::named_contr_sum(mf1[[j]], 0.5, FALSE)
}

test_that("scale works", {
  expect_equal(nauf_model.frame(form1, dat, subset = dat$x < 2, offset = os,
    weights = w, ncs_scale = 0.5), mf1)
})

attr(form1, "standardized.scale") <- 0.5
test_that("standardized scale works", {
  expect_equal(nauf_model.frame(form1, dat, subset = dat$x < 2, offset = os,
    weights = w), mf1)
  expect_warning(expect_equal(nauf_model.frame(form1, dat, subset = dat$x < 2,
    offset = os, weights = w, ncs_scale = 1), mf2))
})

test_that("ignored argument warnings work", {
  expect_warning(nauf_model.frame(form1, dat, contrasts = list()))
  expect_warning(nauf_model.frame(form1, dat, xlev = list()))
  expect_warning(nauf_model.frame(form1, dat, drop.unused.levels = FALSE))
  expect_warning(nauf_model.frame(form1, dat, na.action = na.omit))
})

test_that("zero intercept errors work", {
  expect_error(nauf_model.frame(y ~ 0, dat))
  expect_error(nauf:::check_groups(y ~ x + (0 + x | subj)))
  expect_equal(inherits(nauf:::check_groups(y ~ x + (1 + x || subj)), "error"),
    FALSE)
  expect_equal(inherits(nauf:::check_groups(y ~ x + (1 | subj)), "error"),
    FALSE)
})

d <- expand.grid(f1 = c("a", "b", "c"), subj = 1:10)
d$subj[d$subj > 5] <- NA
d$f1[!is.na(d$subj) & d$f1 == "c"] <- NA
d$y <- stats::rnorm(nrow(d))

ni2 <- list(
  resp = "y",
  groups = list(subj = paste(1:5)),
  uf = list(f1 = list(c("a", "b", "c"), c("a", "b"))),
  of = list(),
  num = list(),
  mat = list(),
  extras = character(),
  cc = list(list(), list(f1 = 2)),
  hasna = c(F, T, T),
  ncs_scale = 1
)
names(ni2$hasna) <- c("y", "f1", "subj")

test_that("ranef group main eff works", {
  expect_equal(attr(attr(nauf_model.frame(y ~ f1 + (1 + f1 | subj), d),
    "terms"), "nauf.info"), ni2)
})

rm(list = ls())
