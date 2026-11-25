
# forestcast

<!-- badges: start -->
<!-- badges: end -->

The 'forestcast' package will let users compare multiple models that predict future
forest metrics at the stand level, with a focus on intensively managed loblolly 
pine in the Southeastern US.

Currently it only has one model (Genlob) published by Gyawali and Burkhart in
2015. This model calculates trees per hectare, basal area, and dominant height
at future ages. It also allows for growth responses to thinning and fertilization
to be specified. To run multiple scenarios in one function call, provide vectors
as inputs. All vectors must be the same length or length 1.

## Installation

You can install 'forestcast' from github using the 'remotes' package:

``` r
remotes::install_github("elijahy357/forestcast")
```

## Example

You can use the genlob() function to predict future forest metrics:

``` r
library(forestcast)

# Basic stand with no management
genlob(age1 = 10, age2 = 30, n = 1000, ht = 10, ba = 10, site_index = 18)

# Specifying multiple scenarios:
# Scenario 1 no thin, 112 kg of nitrogen applied per hectare at age 15.
# Scenario 2 no thin, 224 kg of nitrogen applied per hectare at age 15.
# Scenario 3 thin down to 50% of the basal area at age 15, no fertilization.
# Scenario 4 thin down to 40% of the basal area at age 15.
#     336 kg of kg of nitrogen applied per hectare at age 15
genlob(age1 = 10, age2 = 30, n = 1000, ht = 10, ba = 10, site_index = 18,
thin_age = c(NA, NA, 15, 15), thin_ratio = c(1, 1, 0.5, 0.4),
fert_age = c(15, 15, NA, 15), kgNperha = c(112, 224, 0, 336), P = FALSE)
```

