# nauf
An R package for fitting regressions with not applicable values in unordered factors

## Under Construction
This package is still under construction.  I will make an initial release
when I have finished the documentation and unit tests for basic (fixed effects
only) regressions using nauf.

## Package Description
The *nauf* package provides tools for fitting regressions where
unordered factors may have NA values in addition to their
regular levels.  In this package, NA is used to encode when an
unordered factor is truly *not applicable*.  This is different
than "not available" or "missing at random".  The concept applies only to
unordered factors, and indicates that the factor is simply not meaningful
for an observation, or that while the observation may technically be
definable by one of the factor levels, the interpretation of its belonging
to that group isn't the same.
 
For example, in a linguistic study examining
three dialects of the same language, where two dialects have both monolingual
and bilingual speakers, but the third has only monolingual speakers, a binary
factor *bilingual* with levels *FALSE* and *TRUE* doesn't
apply to the third dialect.  While you could code all observations from
the third dialect as *FALSE*, the meaning would not be the same,
since in the first two dialects *FALSE* means
*monolingual as opposed to bilingual*, while in the third it does not
have this meaning.  In this case, the factor could be coded as *NA*
for observations from the third dialect, and the tools in the *nauf*
package could be used to include them in a regression.  They are included
by using sum contrasts for unordered factors, allowing *NA* values
to pass to the model matrix, and then setting all *NA* values to zero
(in sociolinguistics this is called factor slashing).  The treatment of
unordered factors interacting with these *NA* values is more complex
and can often require contrast changes.

The *nauf* package works by implementing new methods for the generics
*model.frame* and *model.matrix* from the *stats* package, and
ensuring that these generics are used by already existing regression fitting
functions (i.e. there is as little interference with the actual model fitting
procedures as possible, further allowing already existing generics like
*predict*, *summary*, etc. to be seemlessly applied to the fitted model objects).

The main function, *nauf_reg*, is a wrapper function that assigns
the formula for a regression an extra (first) class attribute of
*nauf*, and then calls an appropriate regression
fitting function from the *stats* or *MASS* package, causing these
functions to use the *nauf* generics for *model.frame* and
*model.matrix* (see nauf_model_frame and nauf_model_matrix in the R folder).
Currently, any regression which can be fit with the functions *lm*, *glm*,
and *glm.nb* can be fit using *nauf_reg* to allow *NA* values in the unordered
factors.  For the treatment of unordered factors, see *named_contr_sum*,
*nauf_interaction*, and *nauf_contrasts* in the *contrast_functions.R* file.

## Future Development
After getting the initial package for fixed-effects-only models up and running,
the next step is to implement *nauf* methods in *lme4*, allowing the fixed
effects to be coded as described above, and also allowing for random effects
grouping factors to be coded as *NA* (i.e. have the random effects structure
apply only to a subset of the data while still including all of the data
in the same regression).  There are many tools in the *nlme* package for fitting
more complex random effects structures, but my goal here is to make these
automated so that all the user has to do is code factors as *NA* when they are
not applicable, and the rest simply follows. Once the *lme4* interfacing is
working, I then plan to implement it into *rstanarm*.
