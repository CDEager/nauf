context("nauf_model_frame")

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
dat$x <- stats::rnorm(nrow(dat))
dat$z <- stats::rnorm(nrow(dat))
dat$y <- stats::rnorm(nrow(dat))
dat$w <- w <- stats::runif(nrow(dat))
dat$os <- os <- stats::rnorm(nrow(dat))

form1 <- y ~ f2 * f4 * (f1 + f3) + f1 * f3 * (x + o + poly(z, 3))
form2 <- y ~ f2 * f4 * (f1 + f3) + f1 * f3 * (x + o + poly(z, 3)) + offset(os)
form3 <- y ~ f1 * f3 * f4 + f1 * f3 * (x + o + poly(z, 3))

mf1 <- stats::model.frame(form1, dat, subset = dat$x < 2, offset = os,
  weights = w, na.action = "na.pass")
mf2 <- stats::model.frame(form2, dat, subset = dat$x < 2,
  weights = w, na.action = "na.pass")

for (j in c("f1", "f2", "f3", "f4")) {
  mf1[, j] <- named_contr_sum(mf1[, j], FALSE)
  mf2[, j] <- named_contr_sum(mf2[, j], FALSE)
  attr(attr(mf1, "terms"), "dataClasses")[j] <- "factor"
}
contrasts(mf2$o) <- contrasts(mf1$o) <- contr.poly(3)

hasna <- c(T, T, F, T, F, F, F)
names(hasna) <- c("f2", "f4", "f1", "f3", "x", "o", "poly(z, 3)")

mefc <- list()
mefc$f2 <- list(ordered = FALSE, levels = c("TRUE", "FALSE"),
  contrasts = named_contr_sum(c("TRUE", "FALSE")))
mefc$f4 <- list(ordered = FALSE, levels = c("2", "3"),
  contrasts = named_contr_sum(c("2", "3")))
mefc$f1 <- list(ordered = FALSE, levels = c("a", "b", "c", "d"),
  contrasts = named_contr_sum(c("a", "b", "c", "d")))
mefc$f3 <- list(ordered = FALSE, levels = c("u", "v", "w"),
  contrasts = named_contr_sum(c("u", "v", "w")))
mefc$o <- list(ordered = TRUE, levels = c("1", "2", "3"),
  contrasts = contr.poly(3))

ccna <- list()
ccna[["f1:f2"]] <- list(factors = list(f1 = list(), f2 = list()))
ccna[["f1:f3"]] <- list(factors = list(f1 = list(), f3 = list()))
ccna[["f1:f2:f4"]] <- list(factors = list(f1 = list(), f2 = list(), f4 = list()))
for (i in 1:length(ccna)) {
  for (j in 1:length(ccna[[i]]$factors)) {
    ccna[[i]]$factors[[j]]$ordered <- FALSE
  }
}
ccna[["f1:f2"]]$factors$f1$levels <- c("a", "b", "c")
ccna[["f1:f2"]]$factors$f2$levels <- c("FALSE", "TRUE")
ccna[["f1:f3"]]$factors$f1$levels <- c("c", "d")
ccna[["f1:f3"]]$factors$f3$levels <- c("u", "v", "w")
ccna[["f1:f2:f4"]]$factors$f1$levels <- c("a", "b", "c")
ccna[["f1:f2:f4"]]$factors$f2$levels <- c("FALSE", "TRUE")
ccna[["f1:f2:f4"]]$factors$f4$levels <- c("2", "3")
for (i in 1:length(ccna)) {
  for (j in 1:length(ccna[[i]]$factors)) {
    ccna[[i]]$factors[[j]]$contrasts <- named_contr_sum(
      ccna[[i]]$factors[[j]]$levels)
  }
}

attr(attr(mf1, "terms"), "hasna") <- hasna
attr(attr(mf1, "terms"), "mefc") <- mefc
attr(attr(mf1, "terms"), "ccna") <- ccna
attr(mf1, "terms") <- nauf_on(attr(mf1, "terms"))
attr(attr(mf2, "terms"), "hasna") <- hasna
attr(attr(mf2, "terms"), "mefc") <- mefc
attr(attr(mf2, "terms"), "ccna") <- ccna
attr(mf2, "terms") <- nauf_on(attr(mf2, "terms"))
mf1 <- nauf_on(mf1)
mf2 <- nauf_on(mf2)

mt1 <- attr(mf1, "terms")
mt2 <- attr(mf2, "terms")

mfa <- stats::model.frame(~ x, dat, na.action = "na.pass")
hasna <- FALSE
names(hasna) <- "x"
attr(attr(mfa, "terms"), "hasna") <- hasna
attr(attr(mfa, "terms"), "mefc") <- list()
attr(attr(mfa, "terms"), "ccna") <- list()
attr(mfa, "terms") <- nauf_on(attr(mfa, "terms"))
mfa <- nauf_on(mfa)

mfb <- stats::model.frame(~ f2, dat, na.action = "na.pass")
hasna <- TRUE
names(hasna) <- "f2"
attr(attr(mfb, "terms"), "hasna") <- hasna
attr(attr(mfb, "terms"), "mefc") <- mefc["f2"]
attr(attr(mfb, "terms"), "ccna") <- list()
attr(mfb, "terms") <- nauf_on(attr(mfb, "terms"))
mfb <- nauf_on(mfb)
mfb$f2 <- named_contr_sum(mfb$f2, FALSE)

test_that("default behavior works", {
  expect_equal(nauf_model_frame(form1, dat,
    subset = dat$x < 2, offset = os, weights = w), mf1)
  expect_equal(nauf_model_frame(form2, dat,
    subset = dat$x < 2, weights = w), mf2)
  expect_equal(nauf_model_frame(mt1, dat,
    subset = dat$x < 2, offset = os, weights = w), mf1)
  expect_equal(nauf_model_frame(mt2, dat,
    subset = dat$x < 2, weights = w), mf2)
  expect_equal(nauf_model_frame(nauf_off(mt1), dat,
    subset = dat$x < 2, offset = os, weights = w), mf1)
  expect_equal(nauf_model_frame(nauf_off(mt2), dat,
    subset = dat$x < 2, weights = w), mf2)
  expect_equal(nauf_model_frame(~ x, dat), mfa)
  expect_equal(nauf_model_frame(~ f2, dat), mfb)
  expect_equal(is.nauf(nauf_model_frame(~ f2, dat)), TRUE)
})

test_that("warnings work", {
  expect_warning(expect_equal(nauf_model_frame(form1, dat,
    subset = dat$x < 2, offset = os, weights = w, na.action = "na.omit"), mf1))
  expect_warning(expect_equal(nauf_model_frame(form1, dat,
    subset = dat$x < 2, offset = os, weights = w, xlev = list(f1 = "b")), mf1))
  expect_warning(expect_equal(nauf_model_frame(form1, dat,
    subset = dat$x < 2, offset = os, weights = w,
    drop.unused.levels = FALSE), mf1))
  expect_warning(nauf_model_frame(y ~ f1:x, dat))
})

dat$o[1:10] <- NA

test_that("errors work", {
  expect_error(nauf_model_frame(~ 0, dat))
  expect_error(nauf_model_frame(form1, dat))
})

rm(list = ls())
