context("nauf_pmmeans")

set.seed(1)

d <- expand.grid(
  f1 = c("a", "b", "c"),
  f2 = c("d", "e", "f"),
  f3 = c("g", "h", "i"))
d$f2[d$f1 == "b" & d$f2 == "d"] <- NA
d$f2[d$f1 == "a"] <- NA
d <- rbind(d, d)
d$x <- stats::rnorm(nrow(d), 3, 1)
mm <- nauf_model_matrix(~ f1 + f2 + f3 + x, data = d)
b <- rnorm(ncol(mm))
d$y <- as.vector(mm %*% b) + rnorm(nrow(d))

mod <- nauf_reg(y ~ f1 + f2 + f3 + x, data = d)
b <- summary(mod)$coefficients[, 1]
b0 <- b[1] + mean(d$x) * b["x"]

rg <- nauf_grid(mod)

f1 <- as.numeric(c(b0 + b["f1a"],
  b0 + b["f1b"] - 0.5 * b["f2d"],
  b0 - b["f1a"] - b["f1b"]))
f1_keep_group2 <- f1_keep_group <- list(
  a = list(f1 = "a", f2 = NA),
  b = list(f1 = "b", f2 = c("e", "f")),
  c = list(f1 = "c"))
f1_drop_group <- list(
  b = list(f1 = "b", f2 = c(NA, "d"))
)
f1_keep_group2$a <- f1_keep_group2$a[-2]
f1_pmm_keep_group <- nauf_pmmeans(rg, "f1", keep_group = f1_keep_group)
f1_pmm_drop_group <- nauf_pmmeans(rg, "f1", drop_group = f1_drop_group)
f1_pmm_keep_group2 <- nauf_pmmeans(rg, "f1", keep_group = f1_keep_group2)
f1_pmm_npw <- nauf_pmmeans(rg, "f1", keep_group = f1_keep_group,
  pairwise = FALSE)

test_that("keep_group and drop_group work", {
  expect_equal(summary(f1_pmm_keep_group[[1]])[, 2], f1)
  expect_equal(summary(f1_pmm_keep_group2[[1]])[, 2], f1)
  expect_equal(summary(f1_pmm_drop_group[[1]])[, 2], f1)
})

test_that("pairwise works", {
  expect_equal(f1_pmm_keep_group[[1]], f1_pmm_npw[[1]])
  expect_equal(length(f1_pmm_keep_group), 2)
  expect_equal(length(f1_pmm_npw), 1)
})

f1 <- f1[-1]
f1_pmm_keep <- nauf_pmmeans(mod, "f1", keep_level = list(f1 = c("b", "c")),
  keep_group = f1_keep_group)
f1_pmm_drop <- nauf_pmmeans(mod, "f1", drop_level = list(f1 = "a"),
  keep_group = f1_keep_group)

test_that("mod method, keep_level and drop_level work", {
  expect_equal(summary(f1_pmm_keep[[1]])[, 2], f1)
  expect_equal(summary(f1_pmm_drop[[1]])[, 2], f1)
})

d <- expand.grid(f1 = c("a", "b", "c"), f2 = c("d", "e"))
d$f2[d$f1 == "b"] <- NA
d <- rbind(d, d, d, d)
d$x1 <- stats::rnorm(nrow(d))
d$x2 <- stats::rnorm(nrow(d))
mm <- nauf_model_matrix(~ (f1 + f2) * x1 * x2, data = d)
b <- stats::rnorm(ncol(mm))
d$y <- as.vector(mm %*% b) + stats::rnorm(nrow(d))
mod <- nauf_reg(y ~ (f1 + f2) * x1 * x2, data = d)
b <- summary(mod)$coefficients[, 1]
f1x1 <- as.numeric(c(b["f1a:x1"], b["f1b:x1"],
  -1 * (b["f1a:x1"] + b["f1b:x1"])))
f2x1x2 <- as.numeric(c(b["f2d:x1:x2"], -1 * b["f2d:x1:x2"]))
f1x1_pmm <- nauf_pmmeans(mod, c("f1", "x1"))
f2x1x2_pmm <- nauf_pmmeans(mod, c("f2", "x1", "x2"))

test_that("covariate method works", {
  expect_equal(nrow(summary(f1x1_pmm)[[1]]), 3)
  expect_equal(nrow(summary(f2x1x2_pmm)[[1]]), 2)
  expect_equal(summary(f1x1_pmm)[[1]][, 2], f1x1)
  expect_equal(summary(f2x1x2_pmm)[[1]][, 2], f2x1x2)
})

d <- expand.grid(f1 = c("a", "b", "c"), f2 = c("d", "e"))
d$f2[d$f1 == "a"] <- NA
d <- rbind(d, d, d)
mm <- nauf_model_matrix(~ f1 * f2, data = d)
b <- stats::rnorm(ncol(mm))
d$y <- as.vector(mm %*% b) + stats::rnorm(nrow(d))
mod <- nauf_reg(y ~ f1 * f2, data = d)
b <- summary(mod)$coefficients[, 1]
f1f2_pmm <- nauf_pmmeans(mod, c("f1", "f2"), keep_group = list(
  list(f1 = "a", f2 = NA), list(f1 = c("b", "c"), f2 = c("d", "e"))))
f1f2 <- as.numeric(c(
  b[1] + b[3] + b[4] + b[5],
  b[1] - b[2] - b[3] + b[4] - b[5],
  b[1] + b[3] - b[4] - b[5],
  b[1] - b[2] - b[3] - b[4] + b[5],
  b[1] + b[2]))

test_that("factor interactions work", {
  expect_equal(summary(f1f2_pmm)[[1]][, 2], f1f2)
})


rm(list = ls())
