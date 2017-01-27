context("nauf_reg")

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
dat$bin <- stats::rbinom(nrow(dat), 1, stats::plogis(dat$mod / 2))
dat$poi <- stats::rpois(nrow(dat), exp((dat$mod + 8) / 5))
dat$ngb <- stats::rnbinom(nrow(dat), size = 15,
  prob = stats::plogis(dat$mod / 2))

m_gau <- nauf_reg(gau ~ f1 * f2 * (o + x + poly(z, 2)), dat)
m_bin <- nauf_reg(bin ~ f1 * f2 * (o + x + poly(z, 2)), dat, family = binomial)
m_poi <- nauf_reg(poi ~ f1 * f2 * (o + x + poly(z, 2)), dat, family = poisson)
m_ngb <- nauf_reg(ngb ~ f1 * f2 * (o + x + poly(z, 2)), dat, family = "negbin")

test_that("linear works", {
  expect_silent(nauf_reg(gau ~ f1 * f2 * (o + x + poly(z, 2)), dat))
  expect_equal(names(coef(m_gau)), colnames(mm))
  expect_equal(anova(m_gau)[c(1, 6), 1], c(3, 1))
  expect_equal(m_gau$x, mm)
  expect_equal(predict(m_gau, newdata = dat), m_gau$fitted)
  expect_equal(m_gau$model, nauf_model_frame(
    gau ~ f1 * f2 * (o + x + poly(z, 2)), dat))
  expect_equal(model.frame(m_gau), m_gau$model)
  expect_equal(model.matrix(m_gau), m_gau$x)
})

test_that("binomial works", {
  expect_silent(nauf_reg(bin ~ f1 * f2 * (o + x + poly(z, 2)),
    dat, family = binomial))
  expect_equal(names(coef(m_bin)), colnames(mm))
  expect_equal(anova(m_bin)[c(2, 7), 1], c(3, 1))
  expect_equal(m_bin$x, mm)
  expect_equal(predict(m_bin, newdata = dat), qlogis(m_bin$fitted))
  expect_equal(m_bin$model, nauf_model_frame(
    bin ~ f1 * f2 * (o + x + poly(z, 2)), dat))
  expect_equal(model.frame(m_bin), m_bin$model)
  expect_equal(model.matrix(m_bin), m_bin$x)
})

test_that("poisson works", {
  expect_silent(nauf_reg(poi ~ f1 * f2 * (o + x + poly(z, 2)),
    dat, family = poisson))
  expect_equal(names(coef(m_poi)), colnames(mm))
  expect_equal(anova(m_poi)[c(2, 7), 1], c(3, 1))
  expect_equal(m_poi$x, mm)
  expect_equal(predict(m_poi, newdata = dat), log(m_poi$fitted))
  expect_equal(m_poi$model, nauf_model_frame(
    poi ~ f1 * f2 * (o + x + poly(z, 2)), dat))
  expect_equal(model.frame(m_poi), m_poi$model)
  expect_equal(model.matrix(m_poi), m_poi$x)
})

test_that("negbin works", {
  expect_silent(nauf_reg(ngb ~ f1 * f2 * (o + x + poly(z, 2)),
    dat, family = "negbin"))
  expect_equal(names(coef(m_ngb)), colnames(mm))
  expect_warning(expect_equal(anova(m_ngb)[c(2, 7), 1], c(3, 1)),
    "tests made without re-estimating 'theta'")
  expect_equal(m_ngb$x, mm)
  expect_equal(predict(m_ngb, newdata = dat), log(m_ngb$fitted))
  expect_equal(m_ngb$model, nauf_model_frame(
    ngb ~ f1 * f2 * (o + x + poly(z, 2)), dat))
  expect_equal(model.frame(m_ngb), m_ngb$model)
  expect_equal(model.matrix(m_ngb), m_ngb$x)
})

rm(list = ls())
