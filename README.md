# remx
A package with meta-analytic approximations for relational event models. 

This package supports mixed-effects models, data streams, and large networks. 

# Installation 
```R
install.packages("devtools")
install.packages(c('remstats', 'remstimate', 'remify'))
devtools::install_github("Fabio-Vieira/remx")
```

# Mixed-Effect Models

In this case, there are multiple independent networks that need to be analized using a multilevel model.

The **remx** function is the main function of the package. It fits all multilevel models available.

### Tie-oriented Model

```R

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

```R

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

```R
ranef_metarem(fit)
```

### Fixed effects

```R
fixef_metarem(fit)
```

### Random-effect covariance matrix

```R
VarCov_metarem(fit)
```

# Relational-event Network Data Streams

In this case, we have a relational-event network that is augmented with additional events over time. Every new batch of events will be used to update the model estimates. 

## Computing network statistics

We need to use the arguments **start** and **stop** from **remstats** in order to compute the statistics for a specific portion of the event sequence. For instance, if we wish to fit the model in batches of 500, we need to declare **start** = 1 and **stop** = 500, then **start** = 501 and **stop** = 1000, and so on.

```R

#Now let's simulate a network
edgelist <- stream$edgelist

library(remstats)

#We will devide the network in 10 batches of 500
events <- seq(1, 5001, by = 500)

#Declaring which effects we want remstats to compute
edgelist <- stream$edgelist

# Divide into 10 batches of 500 events
events <- seq(1, 5001, by = 500)

# Effects
effects_stream <- ~ remstats::inertia(scaling = "std") +
  remstats::reciprocity(scaling = "std") +
  remstats::indegreeSender(scaling = "std") +
  remstats::outdegreeSender(scaling = "std")

# Full remify object for the entire sequence
rehObj_stream <- remify::remify(edgelist, model = "tie", actors = stream$actors)

data_stream <- vector("list")
origin_i <- 0  # first batch starts from zero

for(i in 2:length(events)){
  edl <- edgelist[events[i-1]:(events[i]-1), ]
  
  # Build reh with correct origin so intereventTime[1] is small
  reh_i <- remify::remify(edl, model = "tie", 
                          actors = stream$actors,
                          origin = origin_i)
  
  # Compute stats from this reh directly — ensures full compatibility
  stats_i <- remstats::tomstats(effects_stream, reh_i)
  
  data_stream[[i-1]] <- list(edgelist = edl,
                             reh = reh_i,
                             statistics = stats_i)
  
  # Update origin for next batch
  origin_i <- edl[nrow(edl), 1]

```

## Fitting the model

```R
#Let's compute the effects for the first 7 batches of the networks
fit <- strem(data_stream[1:7])

#printing the parameters
print(fit)
```

## Plotting trends of the estimates across batches

The package also contains a generic plot function, that can be used to plot the trends, along with 95\% intervals, of the estimated effects across each batch. This function only works for models fitted with the function **strem**.

```R
#Getting a plot of the estimates with confidence intervals
plot(fit)

```

## Updating the model with new batches of events

We can update our estimates with new batches by simply passing a model previously fitted with the **strem** function to the argument **model** and use the argument **update** = TRUE.

```R
#Now we can update the model with the remaining 3 batches
fit <- strem(data_stream[8:10], update = T, model = fit)

#printing the parameters
print(fit)
```

## Citation

Vieira, F., Leenders, R. & Mulder, J. Fast meta-analytic approximations for relational event models: applications to data streams and multilevel data. J Comput Soc Sc 7, 1823–1859 (2024). https://doi.org/10.1007/s42001-024-00290-7

```R
@article{vieira2024fast,
  title={Fast meta-analytic approximations for relational event models: applications to data streams and multilevel data},
  author={Vieira, Fabio and Leenders, Roger and Mulder, Joris},
  journal={Journal of Computational Social Science},
  pages={1--37},
  year={2024},
  publisher={Springer}
}
```

## Acknowledgement

The development of this package was supported by a Vidi Grant (452-17-006) awarded by the Netherlands Organization for Scientific Research (NWO) Grant and an ERC Starting Grant (758791).
