context("nauf_glmer")

d <- plosives
d$vnpc1 <- prcomp(d[, 1:5], scale = TRUE)$x[, 1]
d <- droplevels(subset(d, voicing == "Voiceless"))
nc <- which(d$dialect != "Cuzco")
d$age[nc] <- NA
d$ling[nc] <- NA
d$ed[nc] <- NA
d$wordpos[!d$spont] <- NA
d$spont[d$dialect == "Valladolid"] <- NA

f <- vnpc1 ~ place + speechrate + dialect * (spont + sex) +
  age + (1 + spont | speaker) +
  (1 + dialect | item)

sdat <- standardize::standardize(f, d)
mod <- nauf_lmer(sdat$formula, sdat$data)
nd <- predict(sdat, subset(d, speaker == "s02"))
preds <- predict(mod, nd)
pfit <- predict(mod)
X <- lme4::getME(mod, "X")
X2 <- X[d$speaker == "s02", ]
fopreds <- predict(mod, nd, re.form = NA)
fofit <- predict(mod, re.form = NA)

test_that("nauf_lmer works", {
  expect_equal(is.nauf.lmerMod(mod), TRUE)
  expect_equal(preds, fitted(mod)[d$speaker == "s02"])
  expect_equal(pfit, fitted(mod))
  expect_equal(fopreds, drop(X2 %*% lme4::fixef(mod)))
  expect_equal(fofit, drop(X %*% lme4::fixef(mod)))
})


f <- I(vdur > 50) ~ dialect + spont + sex + (1 | speaker) + (1 | item)

sdat <- standardize::standardize(f, d, family = binomial)
mod <- nauf_glmer(sdat$formula, sdat$data, family = binomial)
nd <- predict(sdat, subset(d, speaker == "s02"))
preds <- predict(mod, nd)
pfit <- predict(mod)
X <- lme4::getME(mod, "X")
X2 <- X[d$speaker == "s02", ]
fopreds <- predict(mod, nd, re.form = NA)
fofit <- predict(mod, re.form = NA)

test_that("nauf_glmer works", {
  expect_equal(is.nauf.glmerMod(mod), TRUE)
  expect_equal(preds, qlogis(fitted(mod)[d$speaker == "s02"]))
  expect_equal(pfit, qlogis(fitted(mod)))
  expect_equal(fopreds, drop(X2 %*% lme4::fixef(mod)))
  expect_equal(fofit, drop(X %*% lme4::fixef(mod)))
  expect_equal(predict(mod, type = "response"), fitted(mod))
})


f <- vdur ~ dialect + (1 | item)

sdat <- standardize::standardize(f, d, family = "negbin")
mod <- suppressWarnings(nauf_glmer.nb(sdat$formula, sdat$data))
nd <- predict(sdat, subset(d, speaker == "s02"))
preds <- predict(mod, nd)
pfit <- predict(mod)
X <- lme4::getME(mod, "X")
X2 <- X[d$speaker == "s02", ]
fopreds <- predict(mod, nd, re.form = NA)
fofit <- predict(mod, re.form = NA)

test_that("nauf_glmer.nb works", {
  expect_equal(is.nauf.glmerMod(mod), TRUE)
  expect_equal(preds, log(fitted(mod)[d$speaker == "s02"]))
  expect_equal(pfit, log(fitted(mod)))
  expect_equal(fopreds, drop(X2 %*% lme4::fixef(mod)))
  expect_equal(fofit, drop(X %*% lme4::fixef(mod)))
  expect_equal(predict(mod, type = "response"), fitted(mod))
})

rm(list = ls())
