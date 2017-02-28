
#' @importClassesFrom lme4 merMod lmerMod
setClass("nauf_lmermod",
  contains = "lmerMod",
  slots = list(terms = "nauf")
)

#' @importClassesFrom lme4 merMod glmerMod
setClass("nauf_glmermod",
  contains = "glmerMod",
  slots = list(terms = "nauf")
)
