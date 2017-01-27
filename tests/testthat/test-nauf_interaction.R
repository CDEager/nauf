library(nauf)
context("Unordered factor interactions involving NAs")

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
dat <- rbind(dat, dat)
dat$o <- factor(rep(1:3, length.out = nrow(dat)), ordered = TRUE)
dat$x <- rnorm(nrow(dat))

f12 <- list(levels = list(
  f1 = c("a", "b", "c"),
  f2 = c("FALSE", "TRUE")),
changed = TRUE)

f13 <- list(levels = list(
  f1 = c("c", "d"),
  f3 = c("u", "v", "w")),
changed = TRUE)

f14 <- list(levels = list(
  f1 = c("a", "b", "c", "d"),
  f4 = c("2", "3")),
changed = FALSE)

f23 <- list(levels = list(
  f2 = c("FALSE", "TRUE"),
  f3 = c("u", "v", "w")),
changed = FALSE)

f24 <- list(levels = list(
  f2 = c("FALSE", "TRUE"),
  f4 = c("2", "3")),
changed = FALSE)

f34 <- list(levels = list(
  f3 = c("u", "v", "w"),
  f4 = c("2", "3")),
changed = FALSE)

# warning: only one level
f123 <- list(levels = list(
  f1 = "c",
  f2 = c("FALSE", "TRUE"),
  f3 = c("u", "v", "w")),
changed = TRUE)

f124 <- list(levels = list(
  f1 = c("a", "b", "c"),
  f2 = c("FALSE", "TRUE"),
  f4 = c("2", "3")),
changed = TRUE)

# warning: empty cells
f134 <- list(levels = list(
  f1 = c("c", "d"),
  f3 = c("u", "v", "w"),
  f4 = c("2", "3")),
changed = TRUE)

f234 <- list(levels = list(
  f2 = c("FALSE", "TRUE"),
  f3 = c("u", "v", "w"),
  f4 = c("2", "3")),
changed = FALSE)

# warning: only one level
f1234 <- list(levels = list(
  f1 = "c",
  f2 = c("FALSE", "TRUE"),
  f3 = c("u", "v", "w"),
  f4 = c("2", "3")),
changed = TRUE)


test_that("simple unordered factor cases work", {
  expect_equal(nauf_interaction(dat, c(1, 2)), f12)
  expect_equal(nauf_interaction(dat, c(1, 3)), f13)
  expect_equal(nauf_interaction(dat, c(1, 4)), f14)
  expect_equal(nauf_interaction(dat, c(2, 3)), f23)
  expect_equal(nauf_interaction(dat, c(2, 4)), f24)
  expect_equal(nauf_interaction(dat, c(3, 4)), f34)
})

test_that("higher order cases work", {
  expect_equal(nauf_interaction(dat, c(1, 2, 3)), f123)
  expect_equal(nauf_interaction(dat, c(1, 2, 4)), f124)
  expect_equal(nauf_interaction(dat, c(1, 3, 4)), f134)
  expect_equal(nauf_interaction(dat, c(2, 3, 4)), f234)
  expect_equal(nauf_interaction(dat), f1234)
})

test_that("warnings work", {
  expect_warning(nauf_interaction(dat, c(1, 2, 3)))
  expect_warning(nauf_interaction(dat, c(1, 3, 4)))
  expect_warning(nauf_interaction(dat))
})

test_that("errors work", {
  expect_error(nauf_interaction(dat, "d"))
  expect_error(nauf_interaction(dat, c("f1", "x")))
})

rm(dat, f12, f13, f14, f23, f24, f34, f123, f124, f134, f234, f1234)
