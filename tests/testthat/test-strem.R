library(remx)
library(remstats)
library(testthat)

test_that("Testing the data streaming meta analytic approximation",{
  
  edgelist <- stream$edgelist
  events <- seq(1, 5001, by = 500)
  effects <- ~ remstats::inertia(scaling = "std") +
    remstats::reciprocity(scaling = "std") +
    remstats::indegreeSender(scaling = "std") +
    remstats::outdegreeSender(scaling = "std")
  
  data <- vector("list")
  origin_i <- 0
  for(i in 2:length(events)){
    edl <- edgelist[events[i-1]:(events[i]-1), ]
    reh_i <- remify::remify(edl, model = "tie",
                            actors = stream$actors,
                            origin = origin_i)
    stats_i <- remstats::tomstats(effects, reh_i)
    data[[i-1]] <- list(edgelist = edl,
                        reh = reh_i,
                        statistics = stats_i)
    origin_i <- edl[nrow(edl), 1]
  }
  
  fit <- expect_no_error(strem(data[1:7]))
})
