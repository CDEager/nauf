context("nauf_model.matrix")

set.seed(1)
dat <- rbind(
  expand.grid(f1 = c("a", "b"), f2 = c(TRUE, FALSE), stringsAsFactors = FALSE),
  expand.grid(f1 = c("c", "d"), f2 = NA), stringsAsFactors = FALSE)
rownames(dat) <- NULL
dat <- as.data.frame(lapply(dat, function(n) rep(n, 10)))
dat$o <- factor(rep(1:3, length.out = nrow(dat)), ordered = TRUE)
dat$x <- stats::rnorm(nrow(dat))
dat$z <- stats::rnorm(nrow(dat))
dat$y <- stats::rnorm(nrow(dat))

form <- y ~ f1 * f2 * (o + x + poly(z, 2))
mf <- nauf_model.frame(form, dat)
mm <- stats::model.matrix(form, mf)
mf2 <- mf
mf2$f1 <- factor(mf2$f1, levels = c("a", "b"))
mf2$f1 <- standardize::named_contr_sum(mf2$f1, 1, FALSE)
colnames(contrasts(mf2$f1)) <- paste0(".c2.", colnames(contrasts(mf2$f1)))
mm2 <- stats::model.matrix(form, mf2)

asgn <- attr(mm, "assign")
d <- which(asgn %in% attr(attr(mf, "terms"), "nauf.info")$cc[[1]][[1]]$assign)
mm <- mm[, -d]
asgn <- asgn[-d]
mm <- cbind(mm, mm2[, c(9, 20:24)])
asgn <- c(asgn, 6, 13, 13, 14, 15, 15)
names(asgn) <- colnames(mm)
asgn <- sort(asgn)
mm <- mm[, names(asgn)]
names(asgn) <- NULL
mm[is.na(mm)] <- 0
attr(mm, "contrasts") <- NULL
attr(mm, "assign") <- asgn

test_that("all three methods work", {
  expect_equal(nauf_model.matrix(form, dat), mm)
  expect_equal(nauf_model.matrix(mf), mm)
  expect_equal(nauf_model.matrix(data = mf), mm)
  expect_equal(nrow(mm), nrow(dat))
})

test_that("errors work", {
  expect_error(nauf_model.matrix(stats::model.frame(form, dat)))
})

rm(list = ls())
