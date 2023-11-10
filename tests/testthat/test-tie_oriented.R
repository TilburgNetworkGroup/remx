library(remx)
library(remstats)
library(testthat)

test_that("Testing the tie-oriented meta analytic approximation",{
  edgelist <- networks
  edgelist <- lapply(edgelist, function(x) x)
  for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}
  #actors <- lapply(1:length(edgelist), function(x) as.character(1:10))
  rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "tie"))
  effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
  stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, rehObj[[x]]))

  fit <- expect_no_error(
         remx(reh = rehObj,
              statistics = stats,
              random = c("baseline"))
  )
})
