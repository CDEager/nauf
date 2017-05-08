context("nauf_glmer")

d <- fricatives
d$uvoi[!(d$lang == "Catalan" & d$wordpos == "Medial")] <- NA
d$c_speaker <- d$s_speaker <- d$speaker
d$c_speaker[d$lang != "Catalan"] <- NA
d$s_speaker[d$lang != "Spanish"] <- NA
d$fv <- d$pvoi == 1

sdat <- standardize(pvoi ~ lang * wordpos + uvoi +
  (1 + wordpos + uvoi | c_speaker) + (1 + wordpos | s_speaker),
  d)

mod <- nauf_lmer(sdat$formula, sdat$data)

preds <- predict(mod, sdat$data)
fitt <- fitted(mod)
pfit <- predict(mod)
X <- lme4::getME(mod, "X")
fopreds <- predict(mod, sdat$data, re.form = NA)
fofit <- predict(mod, re.form = NA)
fe <- drop(X %*% lme4::fixef(mod))

test_that("nauf_lmer works", {
  expect_equal(preds, fitt)
  expect_equal(pfit, fitt)
  expect_equal(fopreds, fe)
  expect_equal(fofit, fe)
})

rm(list = ls())
