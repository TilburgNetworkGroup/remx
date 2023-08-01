# remx
A package for mixed-effects relational event models. 

This package fits meta-analytic approximations of mixed-effect relational event models. Full models are also supported by **lme4** (frequentist) and **rstanarm** (Bayesian) packages. 

# Installation 
```r{}
install.packages("devtools")
devtools::install_github("TilburgNetworkGroup/remify") 
devtools::install_github("TilburgNetworkGroup/remstimate")
devtools::install_github("TilburgNetworkGroup/remx")
```

# Usage

## Tie-oriented model

```r{}
#----------------------------#
#     Tie-Oriented model     #
#----------------------------#

#Loading libraries
library(remx)
library(remstats)

#Loading the data
edgelist <- networks
edgelist <- lapply(edgelist, function(x) x[,])
for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}

#Computing statistics
effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
reh_obj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "tie"))
stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, reh_obj[[x]]))

#Running the model

fit <- remx(edgelist = edgelist,
            statistics = stats,
            random = c("baseline"),
            fixed = c("inertia", "reciprocity"),
            method = "meta",
            model = "tie")
print(fit)
```

## Actor-oriented model

```r{}
#----------------------------#
#    Actor-Oriented model    #
#----------------------------#

#Computing statistics
sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
reh_obj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "actor"))
stats <- lapply(1:length(edgelist), function(x) remstats::aomstats(reh_obj[[x]], sender_effects, receiver_effects))

#Running the model
fit <- remx(edgelist = edgelist,
            statistics = stats,
            random_sender = c("baseline"),
            fixed_sender = c("indegreeSender"),
            random_receiver = c("indegreeReceiver", "rrankSend"),
            fixed_receiver = NULL,
            method = "meta",
            model = "actor")

print(fit)
```

# Extracting model parameters

These functions are inspired by the functions contained in the **lme4** package.

## Random effects

```r{}
ranef_metarem(fit)
```

## Fixed effects

```r{}
fixef_metarem(fit)
```

## Random-effect covariance matrix

```r{}
VarCov_metarem(fit)
```
