# remx
A package for mixed-effects relational event models. 

This package fits meta-analytic approximations of mixed-effect relational event models. Full models are also support from **lme4** (frequentist) and **rstanarm** (Bayesian) packages. 

# Installation 

```r
remotes::install_github(repo = "TilburgNetworkGroup/remx")
```

# Usage 

```r
#Loading libraries
library(remx)
library(remstats)

#Loading the data
edgelist <- networks
edgelist <- lapply(edgelist, function(x) x[,])
for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}

#Computing statistics
effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "tie"))
stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, rehObj[[x]]))

#Running the model

fit <- remx(edgelist = edgelist,
            statistics = stats,
            random = c("baseline"),
            fixed = c("inertia", "reciprocity"),
            method = "meta",
            model = "tie")

print(fit)
```
