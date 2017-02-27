# nauf
An R package for fitting regressions with not applicable values in unordered factors

## Package Description
The *nauf* package provides tools for fitting regressions where
unordered factors may have NA values in addition to their
regular levels.  In this package, NA is used to encode when an
unordered factor is truly *not applicable*.  This is different
than "not available" or "missing at random".  The concept applies only to
unordered factors, and indicates that the factor is simply not meaningful
for an observation, or that while the observation may technically be
definable by one of the factor levels, the interpretation of its belonging
to that group isn't the same.  For imbalanced observational data, coding
as NA may also be used to control for a factor which is only
contrastive within a subset of the data due to the sampling scheme.  The
*nauf* package provides functions to implement these values
automatically; the user simply needs to code the relevant observations as
NA.

The *nauf* package works by implementing new methods for the generics
*model.frame* and *model.matrix* from the *stats* package,
ensuring that these generics are used by already existing regression fitting
functions. All unordered factors are coded with sum contrasts using
*named_contr_sum*, which uses the contrasts returned by
*contr.sum*, but names the dummy variables rather than
numbering them so the output is more easily interpreted.  These contrasts
are assigned ignoring NA values, and then in the model matrix
all NA values are set to zero.  Interaction terms
which involve more than one unordered factor, with at least one of these
unordered factors having NA values, may require different contrasts
to be used in the interaction term than in the main effects, as determined
with *nauf_interaction*.  A list of the contrast matrices used
in a regression can be obtained with *nauf_contrasts*.

The function *nauf_reg*, is a wrapper function that assigns
the formula for a regression an extra (first) class attribute of
*nauf*, and then calls an appropriate regression
fitting function from the *stats* or *MASS* packages, causing these
functions to use the *nauf* generics for *model.frame* and
*model.matrix* (see *nauf_model_frame* and
*nauf_model_matrix*).  Currently, any regression which can be
fit with the functions *lm*, *glm*,
and *glm.nb* can be fit using *nauf_reg* to
allow NA values in the unordered factors.  The functions
*nauf_grid* and *nauf_pmmeans* interface with
*ref.grid* and *lsmeans* to
obtain predicted marginal means (also called least-squares means) and
pairwise comparisons for factors, interactions between factors, and
interactions between factors and covariates, allowing for the estimates
to be generated for only the relevant subsets of the data.

## Future Development
In a future release, *nauf* methods will be implemented for *lme4*, allowing the
fixed effects to be coded as described above, and also allowing for random effects
grouping factors to be coded as NA (i.e. have the random effects structure
apply only to a subset of the data while still including all of the data
in the same regression).  There are many tools in the *nlme* package for fitting
more complex random effects structures, but the goal here is to make these
automated so that all the user has to do is code factors as NA when they are
not applicable, and the rest simply follows.  Bayesian regression functionality
will also be implemented for *rstanarm*.
