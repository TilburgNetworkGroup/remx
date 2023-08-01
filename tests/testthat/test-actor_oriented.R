library(remx)
library(remstats)
library(testthat)

test_that("Testing the actor-oriented meta analytic approximation",{
  edgelist <- networks
  edgelist <- lapply(edgelist, function(x) x)
  for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}

  sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
  receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
  stats <- lapply(1:length(edgelist), function(x) remstats::aomstats(edgelist[[x]], sender_effects, receiver_effects))

  fit <- expect_no_error(
    remx(edgelist = edgelist,
         statistics = stats,
         random_sender = c("baseline"),
         fixed_sender = c("indegreeSender"),
         random_receiver = c("indegreeReceiver", "rrankSend"),
         fixed_receiver = NULL,
         method = "meta",
         model = "actor")
  )
})
