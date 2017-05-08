context("nauf_pmmeans")

dat <- fricatives
dat$uvoi[!(dat$lang == "Catalan" & dat$wordpos == "Medial")] <- NA
s.pvoi <- standardize(pvoi ~ lang * wordpos + uvoi + lang * dur, dat)
m.pvoi <- nauf_lm(s.pvoi$formula, s.pvoi$data)
rg.pvoi <- nauf_ref.grid(m.pvoi)
b <- coef(m.pvoi)

pmm <- list()

pmm[[1]] <- nauf_pmmeans(rg.pvoi, "lang", pairwise = TRUE,
  subset = list(
    list(lang = "Catalan", wordpos = "Initial", uvoi = NA),
    list(lang = "Catalan", wordpos = "Medial", uvoi = "Voiceless"),
    list(lang = "Spanish", wordpos = c("Initial", "Medial"), uvoi = NA)
  )
)

pmm[[2]] <- nauf_pmmeans(rg.pvoi, "uvoi", pairwise = TRUE, na_as_level = "uvoi")

pmm[[3]] <- nauf_pmmeans(rg.pvoi, "dur")

pmm[[4]] <- nauf_pmmeans(rg.pvoi, c("lang", "dur"), pairwise = TRUE)


pmm <- lapply(lapply(pmm, summary), `[[`, 1)


test_that("subset works", {
  expect_equal(pmm[[1]]$pmmean, as.numeric(c(
    b[1] + b[3] - (b[4] + b[2] + b[7]) / 2,
    b[1] - b[3] - (b[4] - b[7]) / 2)))
})


test_that("na_as_level works", {
  expect_equal(pmm[[2]]$pmmean, as.numeric(b[1] + c(b[2], -b[2], 0)))
})


test_that("covariate works", {
  expect_equal(pmm[[3]]$pmmean, as.numeric(b[6]))
})


test_that("factor covariate works", {
  expect_equal(pmm[[4]]$pmmean, as.numeric(b[6] + c(b[9], -b[9])))
})


rm(list = ls())
