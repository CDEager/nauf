context("nauf_class and nauf_contrasts")

set.seed(1)
dat <- rbind(
  expand.grid(f1 = c("a", "b"), f2 = c(TRUE, FALSE), stringsAsFactors = FALSE),
  expand.grid(f1 = c("c", "d"), f2 = NA), stringsAsFactors = FALSE)
rownames(dat) <- NULL
dat <- as.data.frame(lapply(dat, function(n) rep(n, 50)))
dat$o <- factor(sample(rep(1:3, length.out = nrow(dat))), ordered = TRUE)
dat$x <- stats::rnorm(nrow(dat))
dat$z <- stats::rnorm(nrow(dat))
b <- stats::rnorm(36)
e <- stats::rnorm(nrow(dat))

form <- ~ f1 * f2 * (o + x + poly(z, 2))
mm <- nauf_model_matrix(form, dat)
dat$mod <- as.vector(mm %*% b)
dat$gau <- dat$mod + e
m <- nauf_reg(gau ~ f1 * f2 * (o + x + poly(z, 2)), dat)

mefc <- attr(m$terms, "mefc")
ccna <- attr(m$terms, "ccna")
mefc <- list(
  unordered = list(
    f1 = mefc$f1$contrasts,
    f2 = mefc$f2$contrasts),
  ordered = list(
    o = mefc$o$contrasts))
mefc$unordered$f2 <- rbind(mefc$unordered$f2, 0)
rownames(mefc$unordered$f2)[3] <- NA
ccna[[1]] <- ccna[[1]]$factors
ccna[[1]]$f1 <- rbind(ccna[[1]]$f1$contrasts, 0, 0)
rownames(ccna[[1]]$f1)[3:4] <- c("c", "d")
ccna[[1]]$f2 <- mefc$unordered$f2
contr <- list(mefc = mefc, ccna = ccna)

m2 <- nauf_reg(gau ~ x, dat)
contr2 <- list(mefc = list(unordered = list(), ordered = list()), ccna = list())
m3 <- nauf_reg(gau ~ f2, dat)
contr3 <- list(mefc = list(unordered = list(f2 = contr[[1]][[1]][[2]]),
  ordered = list()), ccna = list())

test_that("is.nauf works", {
  expect_equal(is.nauf(m), TRUE)
  expect_equal(is.nauf(m$terms), TRUE)
  expect_equal(is.nauf(m$x), FALSE)
  expect_equal(is.nauf(m$model), TRUE)
  expect_equal(is.nauf(terms(m$model)), TRUE)
  expect_equal(is.nauf(dat), FALSE)
})

test_that("nauf_contrasts works", {
  expect_equal(nauf_contrasts(m), contr)
  expect_equal(nauf_contrasts(m$terms), contr)
  expect_equal(nauf_contrasts(m$model), contr)
  expect_error(nauf_contrasts(m$x))
  expect_equal(nauf_contrasts(m2), contr2)
  expect_equal(nauf_contrasts(m3), contr3)
})

u <- 5

test_that("nauf_class works", {
  expect_is(nauf_on(dat), "nauf")
  expect_equal(nauf_off(nauf_on(dat)), dat)
  expect_warning(nauf_on(u))
})

rm(list = ls())
