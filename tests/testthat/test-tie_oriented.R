library(remx)
library(remstats)
library(testthat)

test_that("Testing the tie-oriented meta analytic approximation",{
  edgelist <- networks
  edgelist <- lapply(edgelist, function(x) x)
  for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}

  effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
  stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, edgelist[[x]]))

  fit <- expect_no_error(
         remx(edgelist = edgelist,
              statistics = stats,
              random = c("baseline"),
              fixed = c("inertia", "reciprocity"),
              method = "meta",
              model = "tie")
  )
})
