

#' @importClassesFrom lme4 lmerMod
nauf.lmerMod <- setClass("nauf.lmerMod", contains = "lmerMod")

#' @importClassesFrom lme4 glmerMod
nauf.glmerMod <- setClass("nauf.glmerMod", contains = "glmerMod")

# S3 classes
# nauf.formula
# nauf.terms
# nauf.frame
# nauf.glm
# nauf.ref.grid

