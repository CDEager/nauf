
<!-- README.md is generated from README.Rmd. Please edit that file -->
nauf 1.1.0
==========

[![Build Status](https://travis-ci.org/CDEager/nauf.svg?branch=master)](https://travis-ci.org/CDEager/nauf)

Installation
------------

To install the *nauf* package, call:

``` r
install.packages("nauf")
```

Package use
-----------

It is often the case that a factor only makes sense in a subset of a dataset (i.e. for some observations a factor may simply not be meaningful), or that with observational datasets there are no observations in some levels of an interaction term. There are also cases where a random effects grouping factor is only applicable in a subset the data, and it is desireable to model the noise introduced by the repeated measures on the group members within the subset of the data where the repeated measures exist. The *nauf* package allows unordered factors and random effects grouping factors to be coded as NA in the subsets of the data where they are not applicable or otherwise not contrastive. Sum contrasts are used for all unordered factors (using **named\_contr\_sum** in the *standardize* package), and then NA values are set to 0. This allows all of the data to be modeled together without creating collinearity or making the output difficult to interpret.

For example, in the **fricatives** dataset, the factor *uvoi* (underlying voicing) is not contrastive in Spanish, and in Catalan it can only be contrastive in for certain word positions, leading to an imbalanced distribution:

``` r
library(nauf)
#> Loading required package: standardize
#> Loading required package: lme4
#> Loading required package: Matrix
#> Loading required package: rstanarm
#> Loading required package: Rcpp
#> rstanarm (Version 2.15.3, packaged: 2017-04-29 06:18:44 UTC)
#> - Do not expect the default priors to remain the same in future rstanarm versions.
#> Thus, R scripts should specify priors explicitly, even if they are just the defaults.
#> - For execution on a local, multicore CPU with excess RAM we recommend calling
#> options(mc.cores = parallel::detectCores())

summary(fricatives)
#>       dur              pvoi             lang        wordpos   
#>  Min.   : 29.71   Min.   :0.0625   Catalan:974   Final  :296  
#>  1st Qu.: 71.14   1st Qu.:0.5263   Spanish:648   Initial:397  
#>  Median : 86.29   Median :0.6471                 Medial :929  
#>  Mean   : 90.35   Mean   :0.6937                              
#>  3rd Qu.:105.53   3rd Qu.:0.9231                              
#>  Max.   :250.78   Max.   :1.0000                              
#>                                                               
#>           uvoi         speaker    
#>  Neutralized: 160   s12    : 142  
#>  Voiced     : 327   s09    : 141  
#>  Voiceless  :1135   s10    : 102  
#>                     s20    :  71  
#>                     s36    :  70  
#>                     s19    :  64  
#>                     (Other):1032

dat <- fricatives

u <- unique(dat[, c("lang", "wordpos", "uvoi")])
u <- u[order(u$lang, u$wordpos, u$uvoi), ]
rownames(u) <- NULL
u
#>      lang wordpos        uvoi
#> 1 Catalan   Final Neutralized
#> 2 Catalan Initial   Voiceless
#> 3 Catalan  Medial      Voiced
#> 4 Catalan  Medial   Voiceless
#> 5 Spanish   Final   Voiceless
#> 6 Spanish Initial   Voiceless
#> 7 Spanish  Medial   Voiceless
```

With *nauf*, we can code *uvoi* as NA when it is not contrastive, and include *uvoi* slopes only for speakers where it is contrastive by creating language-specific speaker columns set to NA for the opposite language:

``` r
u$uvoi[!(u$lang == "Catalan" & u$wordpos == "Medial")] <- NA
u
#>      lang wordpos      uvoi
#> 1 Catalan   Final      <NA>
#> 2 Catalan Initial      <NA>
#> 3 Catalan  Medial    Voiced
#> 4 Catalan  Medial Voiceless
#> 5 Spanish   Final      <NA>
#> 6 Spanish Initial      <NA>
#> 7 Spanish  Medial      <NA>

dat$uvoi[!(dat$lang == "Catalan" & dat$wordpos == "Medial")] <- NA
dat$c_speaker <- dat$s_speaker <- dat$speaker
dat$c_speaker[dat$lang != "Catalan"] <- NA
dat$s_speaker[dat$lang != "Spanish"] <- NA

sdat <- standardize(pvoi ~ lang * wordpos + uvoi +
  (1 + wordpos + uvoi | c_speaker) + (1 + wordpos | s_speaker),
  dat)

mod <- nauf_lmer(sdat$formula, sdat$data)

summary(mod)
#> Linear mixed model fit by REML ['nauf.lmerMod']
#> Formula: 
#> pvoi ~ lang * wordpos + uvoi + (1 + wordpos + uvoi | c_speaker) +  
#>     (1 + wordpos | s_speaker)
#>    Data: sdat$data
#> 
#> REML criterion at convergence: 3584.7
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -4.4721 -0.6148  0.0511  0.5195  3.1344 
#> 
#> Random effects:
#>  Groups    Name           Variance Std.Dev. Corr             
#>  c_speaker (Intercept)    0.154115 0.39258                   
#>            wordposFinal   0.042090 0.20516  -0.67            
#>            wordposInitial 0.034222 0.18499   0.36 -0.91      
#>            uvoiVoiced     0.028117 0.16768  -0.21  0.76 -0.73
#>  s_speaker (Intercept)    0.123763 0.35180                   
#>            wordposFinal   0.011092 0.10532   0.92            
#>            wordposInitial 0.006788 0.08239  -0.38  0.01      
#>  Residual                 0.482529 0.69464                   
#> Number of obs: 1622, groups:  c_speaker, 26; s_speaker, 16
#> 
#> Fixed effects:
#>                            Estimate Std. Error t value
#> (Intercept)                 0.02215    0.06339   0.349
#> langCatalan                 0.28058    0.06339   4.426
#> wordposFinal                0.46883    0.04183  11.208
#> wordposInitial             -0.40194    0.03791 -10.603
#> uvoiVoiced                  0.64890    0.04716  13.761
#> langCatalan:wordposFinal    0.28252    0.04183   6.754
#> langCatalan:wordposInitial -0.39017    0.03791 -10.292
#> 
#> Correlation of Fixed Effects:
#>             (Intr) lngCtl wrdpsF wrdpsI uvoVcd lngC:F
#> langCatalan -0.150                                   
#> wordposFinl  0.105 -0.416                            
#> wordposIntl  0.045  0.155 -0.669                     
#> uvoiVoiced  -0.098 -0.098  0.282 -0.236              
#> lngCtln:wrF -0.416  0.105  0.076 -0.184  0.282       
#> lngCtln:wrI  0.155  0.045 -0.184  0.057 -0.236 -0.669
```

Predicted marginal means can be calculated for specific subsets where a factor is contrastive. For example, *uvoi* is only contrastive for word-medial Catalan fricatives, so we could call:

``` r
rg <- nauf_ref.grid(mod)

nauf_pmmeans(rg, "uvoi", pairwise = TRUE,
  subset = list(lang = "Catalan", wordpos = "Medial")
)
#> 
#> Predicted marginal means for 'uvoi'
#> NA not considered a level for: 'uvoi'
#> 
#> Factors conditioned on: 'lang' 'wordpos' 
#> 
#> See the 'subset' element of the 'specs' attribute for subsetted groups
#> 
#> $pmmeans
#>  uvoi          pmmean        SE    df   lower.CL    upper.CL
#>  Voiced     0.9923974 0.1005598 25.74  0.7855905  1.19920437
#>  Voiceless -0.3054125 0.1163456 23.48 -0.5458214 -0.06500364
#> 
#> Confidence level used: 0.95 
#> 
#> $contrasts
#>  contrast           estimate         SE    df t.ratio p.value
#>  Voiced - Voiceless  1.29781 0.09431177 20.33  13.761  <.0001
```
