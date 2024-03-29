library(remx)
library(remstats)
library(testthat)

test_that("Testing the data streaming meta analytic approximation",{

  edgelist <- stream$edgelist

  #We will devide the network in 10 batches of 500
  events <- seq(1, 5001, by = 500)
  #Declaring which effects we want remstats to compute
  effects <- ~ remstats::inertia(scaling = "std") + 
    remstats::reciprocity(scaling = "std") + 
    remstats::indegreeSender(scaling = "std") +
    remstats::outdegreeSender(scaling = "std")
  #Getting the remify object

  rehObj <- remify::remify(edgelist, model = "tie")

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
  #Let's compute the effects for the first 7 batches of the networks
  fit <-  expect_no_error(
    strem(data[1:7])
  )
}
)

