context("nauf_pmmeans")

xlev <- list(f1 = c("a", "b", NA), f2 = c("d", "e"), f3 = c("g", "h"))
d <- expand.grid(xlev)

lvs1 <- list(f1 = c("a", NA), f3 = "h")
lvs2 <- list(f3 = "h")
lvs3 <- list(f1 = NA)

mr1 <- nauf:::match_row(d, lvs1)
mr2 <- nauf:::match_row(d, lvs2)
mr3 <- nauf:::match_row(d, lvs3)
mr4 <- nauf:::match_row(d, lvs1, "f3")
mr5 <- nauf:::match_row(d, lvs1, f = any)

s1 <- c(7, 9, 10, 12)
s2 <- 7:12
s3 <- c(3, 6, 9, 12)
s4 <- 7:12
s5 <- c(1, 3, 4, 6, 7:12)

keep_group <- list(lvs2)
kg2 <- apply(sapply(keep_group, function(g) {
  nauf:::match_row(d, g)
}), 1, any)
drop_group <- list(lvs3)
dg3 <- !apply(sapply(drop_group, function(g) {
  nauf:::match_row(d, g)
}), 1, any)
drop_group <- list(lvs1)
dg4 <- !apply(sapply(drop_group, function(g) {
  nauf:::match_row(d, g)
}), 1, any)

eg <- data.frame(f1 = NA)
ec1 <- nauf:::estimate_contrasts(d, eg, "f1")
e1 <- list("1" = rep(c(0, 0, 1/4), 4))
eg <- data.frame(f1 = "a", f3 = "g")
ec2 <- nauf:::estimate_contrasts(d, eg, "f1")
e2 <- list("1" = rep(c(1/4, 0, 0), 4))
ec3 <- nauf:::estimate_contrasts(d, eg, "f3")
e3 <- list("1" = rep(c(1/6, 0), each = 6))
ec4 <- nauf:::estimate_contrasts(d, eg, NULL)
e4 <- list("1" = rep(1/12, 12))
eg <- data.frame(f1 = c("a", NA), f3 = c("g", "h"))
ec5 <- nauf:::estimate_contrasts(d, eg, c("f1", "f3"))
e5 <- list("1" = c(1/2, 0, 0, 1/2, rep(0, 8)),
  "2" = c(rep(0, 8), 1/2, 0, 0, 1/2))

keep1 <- nauf:::sub_est_grid(data.frame(f1 = "a"), d, "f1")
keep2 <- nauf:::sub_est_grid(data.frame(f1 = "a", f3 = "g"), d, c("f1", "f3"))
keep3 <- nauf:::sub_est_grid(data.frame(f1 = c("a", NA)), d, "f1")
keep4 <- nauf:::sub_est_grid(data.frame(f1 = c("a", NA), f3 = c("g", "h")), d,
  c("f1", "f3"))

test_that("subsetting works", {
  expect_equal(which(mr1), s1)
  expect_equal(which(mr2), s2)
  expect_equal(which(mr3), s3)
  expect_equal(which(mr4), s4)
  expect_equal(which(mr5), s5)
  expect_equal(kg2, mr2)
  expect_equal(dg3, !is.na(d$f1))
  expect_equal(dg4, !mr1)
  expect_equal(e1, ec1)
  expect_equal(e2, ec2)
  expect_equal(e3, ec3)
  expect_equal(e4, ec4)
  expect_equal(e5, ec5)
  expect_equal(which(keep1), c(1, 4, 7, 10))
  expect_equal(which(keep2), c(1, 4))
  expect_equal(which(keep3), c(1, 3, 4, 6, 7, 9, 10, 12))
  expect_equal(which(keep4), c(1, 4, 9, 12))
})


test_that("check_subset works", {
  expect_equal(nauf:::check_subset(xlev, list(f1 = "a"),
    list(f1 = NA, f2 = c("d", "e")), list(list(f1 = "a")), list(list(f1 = "a"),
    list(f1 = NA, f2 = c("d", "e")))), NULL)
  expect_equal(nauf:::check_subset(xlev, list(), list(), list(), list()), NULL)
  expect_equal(nauf:::check_subset(NULL, list(), list(), list(), list()), NULL)
  expect_error(nauf:::check_subset(NULL, list(f1 = "a"), list(), list(), list()))
  expect_error(nauf:::check_subset(xlev, list(f4 = "a"), list(), list(), list()))
  expect_error(nauf:::check_subset(xlev, "a", list(), list(), list()))
  expect_error(nauf:::check_subset(xlev, list(f1 = "d"), list(), list(), list()))
  expect_error(nauf:::check_subset(xlev, list(), list(), list(), list(
    f1 = "a")))
  expect_error(nauf:::check_subset(xlev, list(), list(), list(), list(
    list())))
  expect_error(nauf:::check_subset(xlev, list(), list(), list(), list(
    list(f1 = "d"))))
  expect_error(nauf:::check_subset(xlev, list(), list(), list(), list(
    list(f1 = NA, f2 = NA))))
  expect_error(nauf:::check_subset(xlev, list(), list(), list(), list(
    list(f1 = NA), list(f4 = "g"))))
})


set.seed(1)

d <- expand.grid(
  f1 = c("a", "b", "c"),
  f2 = c("d", "e", "f"),
  f3 = c("g", "h", "i"))
d$f2[d$f1 == "b" & d$f2 == "d"] <- NA
d$f2[d$f1 == "a"] <- NA
d <- rbind(d, d)
d$x <- stats::rnorm(nrow(d), 3, 1)
d$y <- stats::rnorm(nrow(d))
mod <- nauf_lm(y ~ f1 + f2 + f3 + x, data = d)
b <- summary(mod)$coefficients[, 1]
b0 <- b[1] + mean(d$x) * b["x"]

rg <- nauf_ref.grid(mod)

f1 <- as.numeric(c(b0 + b["f1a"],
  b0 + b["f1b"] - 0.5 * b["f2d"],
  b0 - b["f1a"] - b["f1b"]))
f1_keep_group2 <- f1_keep_group <- list(
  a = list(f1 = "a", f2 = NA),
  b = list(f1 = "b", f2 = c("e", "f")),
  c = list(f1 = "c"))
f1_drop_group <- list(
  b = list(f1 = "b", f2 = c(NA, "d")))
  
f1_keep_group2$a <- f1_keep_group2$a[-2]
f1_pmm_keep_group <- nauf_pmmeans(rg, pairwise ~ f1, keep_group = f1_keep_group)
f1_pmm_drop_group <- nauf_pmmeans(rg, pairwise ~ f1, drop_group = f1_drop_group)
f1_pmm_keep_group2 <- nauf_pmmeans(rg, pairwise ~ f1, keep_group = f1_keep_group2)
f1_pmm_npw <- nauf_pmmeans(rg, ~ f1, keep_group = f1_keep_group)  # failed

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
f1_pmm_keep <- nauf_pmmeans(mod, pairwise ~ f1,
  keep_level = list(f1 = c("b", "c")), keep_group = f1_keep_group)
f1_pmm_drop <- nauf_pmmeans(mod, pairwise ~ f1, drop_level = list(f1 = "a"),
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
d$y <- stats::rnorm(nrow(d))
mod <- nauf_lm(y ~ (f1 + f2) * x1 * x2, data = d)
b <- summary(mod)$coefficients[, 1]

f1x1 <- as.numeric(
  b["x1"] +
  c(b["f1a:x1"], b["f1b:x1"], -1 * (b["f1a:x1"] + b["f1b:x1"])) +
  mean(d$x2) * b["x1:x2"] + mean(d$x2) *
  c(b["f1a:x1:x2"], b["f1b:x1:x2"], -1 * (b["f1a:x1:x2"] + b["f1b:x1:x2"])))
  
f2 <- c(1, -1, 0)
x1x2 <- (mean(d$x1) + 1) * (mean(d$x2) + 1) - mean(d$x1) * mean(d$x2)
f2x1x2 <- as.numeric(
  b["x1"] + b["x2"] + x1x2 * b["x1:x2"] +
  f2 * (b["f2d:x1"] + b["f2d:x2"] + x1x2 * b["f2d:x1:x2"])
)

test_that("factor-covariate and covariate method works", {
  expect_warning(expect_equal(summary(
    nauf_pmmeans(mod, pairwise ~ f1 * x1))[[1]][, 3], f1x1))
  expect_equal(summary(nauf_pmmeans(mod, pairwise ~ f2 : x1 : x2))[[1]][, 4],
    f2x1x2)
  expect_warning(expect_equal(summary(nauf_pmmeans(mod, ~ x1))[[1]][, 2],
    as.numeric(b["x1"] + mean(d$x2) * b["x1:x2"])))
  expect_warning(expect_equal(summary(nauf_pmmeans(mod, ~ x1 + x2))[[1]][, 3],
    as.numeric(b["x1"] + b["x2"] + (mean(d$x1) + mean(d$x2) + 1) * b["x1:x2"])))
})

d <- expand.grid(f1 = c("a", "b", "c"), f2 = c("d", "e"))
d$f2[d$f1 == "a"] <- NA
d <- rbind(d, d, d)
d$y <- stats::rnorm(nrow(d))
mod <- nauf_lm(y ~ f1 * f2, data = d)
b <- summary(mod)$coefficients[, 1]
f1f2_pmm <- nauf_pmmeans(mod, pairwise ~ f1 * f2, keep_group = list(
  list(f1 = "a", f2 = NA), list(f1 = c("b", "c"), f2 = c("d", "e"))))
f1f2 <- as.numeric(c(
  b[1] + b[3] + b[4] + b[5],
  b[1] - b[2] - b[3] + b[4] - b[5],
  b[1] + b[3] - b[4] - b[5],
  b[1] - b[2] - b[3] - b[4] + b[5],
  b[1] + b[2]))

test_that("factor interactions work", {
  expect_equal(summary(f1f2_pmm)[[1]][, 3], f1f2)
})


rm(list = ls())
