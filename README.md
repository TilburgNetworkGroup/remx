# remx
A package with meta-analytic approximations for relational event models. 

This package supports mixed-effects models, data streams, and large networks. 

# Installation 
```r{}
install.packages("devtools")
devtools::install_github("TilburgNetworkGroup/remify") 
devtools::install_github("TilburgNetworkGroup/remstimate")
devtools::install_github("TilburgNetworkGroup/remstats")
devtools::install_github("TilburgNetworkGroup/remx")
```

# Mixed-Effect Models

In this case, there are multiple independent networks that need to be analized using a multilevel model.

The **remx** function is the main function of the package. It fits all multilevel models available.

### Tie-oriented Model

```r{}

#----------------------------#
#     Tie-Oriented model     #
#----------------------------#

#Loading libraries
library(remx)
library(remstats)

#Loading the data
edgelist <- networks

#Computing statistics
effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "tie"))
stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, rehObj[[x]]))

#Running the model

fit <- remx(reh = rehObj,
            statistics = stats,
            random = c("baseline", "inertia", "reciprocity"))

print(fit)

```

### Actor-oriented Model

```r{}

#----------------------------#
#    Actor-Oriented model    #
#----------------------------#

#Computing statistics
rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "actor"))
sender_effects <- ~ indegreeSender(scaling = "std")
receiver_effects <- ~ indegreeReceiver(scaling = "std") + outdegreeReceiver(scaling = "std")
statistics <- lapply(1:length(edgelist), function(x) remstats::aomstats(rehObj[[x]], sender_effects, receiver_effects))

fit <- remx(reh = rehObj,
            statistics = statistics,
            random_sender = c("baseline", "indegreeSender"),
            random_receiver = c("indegreeReceiver", "outdegreeReceiver"))

print(fit)

```

## Extracting model parameters

These functions are inspired by the functions contained in the **lme4** package. The can be used to extract the parameters for models fitted using the **remx** function.

### Random effects

```r{}
ranef_metarem(fit)
```

### Fixed effects

```r{}
fixef_metarem(fit)
```

### Random-effect covariance matrix

```r{}
VarCov_metarem(fit)
```

# Relational-event Network Data Streams

In this case, we have a relational-event network that is augmented with additional events over time. Every new batch of events will be used to update the model estimates. 

# Computing network statistics

We need to use the arguments **start** and **stop** from **remstats** in order to compute the statistics for a specific portion of the event sequence. For instance, if we wish to fit the model in batches of 500, we need to declare **start** = 1 and **stop** = 500, then **start** = 501 and **stop** = 1000, and so on.

```r{}

#Now let's simulate a network
edgelist <- stream$edgelist

library(remstats)

#We will devide the network in 10 batches of 500
events <- seq(1, 5001, by = 500)

#Declaring which effects we want remstats to compute
effects <- ~ remstats::inertia(scaling = "std") + 
  remstats::reciprocity(scaling = "std") + 
  remstats::indegreeSender(scaling = "std") +
  remstats::outdegreeSender(scaling = "std")

#Getting the remify object
rehObj <- remify::remify(edgelist, model = "tie", actors = stream$actors)

data <- vector("list")

for(i in 2:length(events)){
  #Computing statistics
  stats <- remstats::tomstats(effects, rehObj, start = events[i-1], stop = events[i]-1)
  edl <- edgelist[events[i-1]:(events[i]-1),]
  #Every piece needs to be stored a in a list with edgelist and statistics
  data[[i-1]] <- list(edgelist = edl,
                      reh = remify::remify(edl, model = "tie", actors = stream$actors),
                      statistics = stats)
}
```

# Fitting the model

```r{}
#Let's compute the effects for the first 7 batches of the networks
fit <- strem(data[1:7])

#printing the parameters
print(fit)
```

The package also contains a generic plot function, that can be used to plot the trends, along with 95\% intervals, of the estimated effects across each batch. This function only works for models fitted with the function **strem**.

```r{}
#Getting a plot of the estimates with confidence intervals
plot(fit)

```

We can update our estimates with new batches by simply passing a model previously fitted with the **strem** function to the argument **model** and use the argument **update** = TRUE.

```r{}
#Now we can update the model with the remaining 3 batches
fit <- strem(data[8:10], update = T, model = fit)

#printing the parameters
print(fit)
```


