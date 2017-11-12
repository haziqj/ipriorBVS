R/ipriorBVS: Bayesian Variable Selection for Linear Models using I-priors
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
Bayesian variable selection for linear models using I-priors in R. This work is part of the PhD project entitled *Regression Modelling with Priors using Fisher Information Covariance Kernels (I-priors)*. Visit <http://phd.haziqj.ml> for details.

Benchmark data (Tibshirani, 1996)
---------------------------------

A toy data set designed by [Tibshirani (1996)](https://statweb.stanford.edu/~tibs/lasso/lasso.pdf), often used to compare variable selection methods. `n = 20` data points are generated from a linear model with parameters `beta = c(3, 1.5, 0, 0, 2, 0, 0, 0)` and `sigma = 3`. The `X` are generated from a normal distribution with mean zero, and the correlation between the `i`th and `j`th variable is `0.5 ^ abs(i - j)`. This is implemented in the `gen_benchmark()` function included in the package.

``` r
(dat <- gen_benchmark(n = 50, sd = 3, seed = 123))
## n    =  50 
## p    =  8 
## SNR2 =  1.787658
```

### Model fit

The model fitted either using formula or non-formula syntax. We are then able to obtain posterior inclusion probabilities (PIPs) for the each variable, and also posterior model probabilities (PMPs). For comparison, Bayes factors and deviances are reported as well.

``` r
runjags::runjags.options(silent.jags = TRUE, silent.runjags = TRUE)
(mod <- ipriorBVS(y ~ ., dat))
##            PIP     1     2     3     4     5
## X.1      1.000     x     x     x     x     x
## X.2      0.840     x     x     x     x     x
## X.3      0.568           x     x            
## X.4      0.524                 x     x     x
## X.5      0.644     x     x                 x
## X.6      0.294                              
## X.7      0.480                 x            
## X.8      0.238                              
## PMP            0.061 0.048 0.041 0.040 0.037
## BF             1.000 0.785 0.662 0.648 0.604
## Deviance       93.76 92.07 91.42 96.29 94.16
```

### Coefficients

The model coefficients are averaged across all probable sub-models, which yields a kind of "model-averaged" coefficients.

``` r
coef(mod)
##               PIP   Mean   S.D.   2.5%  97.5%
## (Intercept) 1.000 -0.128  0.459 -1.069  0.739
## X.1         1.000  2.707  0.636  1.588  4.053
## X.2         0.840  1.547  0.878  0.000  2.787
## X.3         0.568  0.607  0.705  0.000  2.119
## X.4         0.524  0.468  0.585 -0.002  1.727
## X.5         0.644  0.858  0.903 -0.110  2.672
## X.6         0.294 -0.158  0.399 -1.273  0.324
## X.7         0.480  0.373  0.523 -0.054  1.582
## X.8         0.238  0.054  0.246 -0.333  0.815
```

------------------------------------------------------------------------

Copyright (C) 2017 [Haziq Jamil](http://haziqj.ml).
