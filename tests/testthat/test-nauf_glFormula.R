context("nauf_glFormula")

d <- expand.grid(
  f1 = c("a", "b", "c", "d"),
  f2 = c("e", "f"),
  f3 = c("g", "h", "i"),
  rg = c(1:10, NA)
)
d$f2[d$f1 == "d"] <- NA
d$f1[d$f1 == "c" & !is.na(d$rg)] <- NA
d$f3[d$f3 == "i" & !is.na(d$rg)] <- NA
d$x <- stats::rnorm(nrow(d))
d$y <- stats::rnorm(nrow(d))
d$count <- stats::rpois(nrow(d), 5)

lmod <- nauf_lFormula(y ~ f1 * f2 * x + (1 + f1 * x + f3 | rg), d)
glmod <- nauf_glFormula(count ~ f1 * f2 * x + (1 + f1 * x + f3 | rg), d,
  family = poisson)

xcols <- c("(Intercept)", "f1a", "f1b", "f1c", "f2e", "x",
  "f1.c3.a:f2e", "f1.c3.b:f2e", "f1a:x", "f1b:x", "f1c:x", "f2e:x",
  "f1.c3.a:f2e:x", "f1.c3.b:f2e:x")
cnms <- list(rg = c("(Intercept)", "f1.c2.a", "f1.c2.b", "x", "f3.c2.g",
  "f1.c2.a:x", "f1.c2.b:x"))
flist <- data.frame(rg = d$rg)
flist$rg <- factor(flist$rg)
attr(flist, "assign") <- 1

test_that("nauf_lFormula and nauf_glFormula work", {
  expect_equal(colnames(lmod$X), xcols)
  expect_equal(colnames(glmod$X), xcols)
  expect_equal(lmod$reTrms$cnms, cnms)
  expect_equal(glmod$reTrms$cnms, cnms)
  expect_equal(lmod$reTrms$flist, flist)
  expect_equal(glmod$reTrms$flist, flist)
  expect_equal(nrow(lmod$X), nrow(d))
  expect_equal(nrow(glmod$X), nrow(d))
  expect_equal(ncol(lmod$reTrms$Zt), nrow(d))
  expect_equal(ncol(glmod$reTrms$Zt), nrow(d))
  expect_equal(nrow(lmod$reTrms$Zt), 70)
  expect_equal(nrow(glmod$reTrms$Zt), 70)
  expect_equal(anyNA(lmod$X), FALSE)
  expect_equal(anyNA(glmod$X), FALSE)
  expect_equal(anyNA(lmod$reTrms$Zt), FALSE)
  expect_equal(anyNA(glmod$reTrms$Zt), FALSE)
  expect_equal(all(as.matrix(lmod$reTrms$Zt)[, is.na(d$rg)] == 0), TRUE)
  expect_equal(all(as.matrix(glmod$reTrms$Zt)[, is.na(d$rg)] == 0), TRUE)
})

rm(list = ls())
