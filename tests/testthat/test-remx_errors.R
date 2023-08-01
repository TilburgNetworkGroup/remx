library(remx)
library(remstats)
library(testthat)

test_that("Testing whether the remx function throw expected errors.", {

  edgelist <- networks
  edgelist <- lapply(edgelist, function(x) x)
  for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}

  effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
  stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, edgelist[[x]]))

  #Passing whatever object as edgelist
  expect_error(
      remx(edgelist = "a",
           statistics = stats,
           random = c("baseline"),
           fixed = c("inertia", "reciprocity"),
           method = "meta",
           model = "tie")
  )

  #If there is no multiple event sequences
  expect_error(
    remx(edgelist = list(edgelist[[1]]),
         statistics = stats,
         random = c("baseline"),
         fixed = c("inertia", "reciprocity"),
         method = "meta",
         model = "tie")
  )

  #Passing whatever as the statistics list
  expect_error(
    remx(edgelist = edgelist,
         statistics = "a",
         random = c("baseline"),
         fixed = c("inertia", "reciprocity"),
         method = "meta",
         model = "tie")
  )

  #statistics and edgelist objects do not have the same length
  expect_error(
    remx(edgelist = edgelist[[1]],
         statistics = stats,
         random = c("baseline"),
         fixed = c("inertia", "reciprocity"),
         method = "meta",
         model = "tie")
  )

  #Testing whether all statistics objects are from remstats class
  stats2 <- stats
  stats2[[15]] <- "a"
  expect_error(
    remx(edgelist = edgelist,
         statistics = stats2,
         random = c("baseline"),
         fixed = c("inertia", "reciprocity"),
         method = "meta",
         model = "tie")
  )

  #If model is tie and not all objects are tomstats
  sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
  receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
  stats3 <- lapply(1:length(edgelist), function(x) remstats::aomstats(edgelist[[x]], sender_effects, receiver_effects))
  stats2[[15]] <- stats3[[15]]
  expect_error(
    remx(edgelist = edgelist,
         statistics = stats2,
         random = c("baseline"),
         fixed = c("inertia", "reciprocity"),
         method = "meta",
         model = "tie")
  )

  #If model is actor and not all objects are aomstats
  stats2 <- stats3
  stats2[[15]] <- stats[[15]]
  expect_error(
    remx(edgelist = edgelist,
         statistics = stats,
         random_sender = c("baseline"),
         fixed_sender = c("indegreeSender"),
         random_receiver = c("indegreeReceiver", "rrankSend"),
         fixed_receiver = NULL,
         method = "meta",
         model = "actor")
  )

  #If model is actor and methor is bayes or classic
  expect_error(
    remx(edgelist = edgelist,
         statistics = stats3,
         random = c("baseline"),
         fixed = c("inertia", "reciprocity"),
         method = "bayes",
         model = "actor")
  )

  #If model is tie, but no random or fixed effects
  expect_error(
    remx(edgelist = edgelist,
         statistics = stats,
         random = NULL,
         fixed = NULL,
         method = "meta",
         model = "tie")
  )

  #If model is actor, but no random or fixed effects
  expect_error(
    remx(edgelist = edgelist,
         statistics = stats3,
         random_sender = NULL,
         fixed_sender = NULL,
         random_receiver = NULL,
         fixed_receiver = NULL,
         method = "meta",
         model = "actor")
  )

})
