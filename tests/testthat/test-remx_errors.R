library(remx)
library(remify)
library(remstats)
library(testthat)

test_that("Testing whether the remx function throw expected errors.", {

  edgelist <- networks
  edgelist <- lapply(edgelist, function(x) x)
  for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}
  #actors <- lapply(1:length(edgelist), function(x) as.character(1:10))
  rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "tie"))
  effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
  stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, rehObj[[x]]))

  #Passing whatever object as edgelist
  expect_error(
      remx(reh = "a",
           statistics = stats,
           random = c("baseline"))
  )

  #If there is no multiple event sequences
  expect_error(
    remx(reh = list(edgelist[[1]]),
         statistics = stats,
         random = c("baseline", "inertia", "reciprocity"))
  )

  #Passing whatever as the statistics list
  expect_error(
    remx(reh = rehObj,
         statistics = "a",
         random = c("baseline", "inertia", "reciprocity"))
  )

  #statistics and edgelist objects do not have the same length
  expect_error(
    remx(reh = rehObj[[1]],
         statistics = stats,
         random = c("baseline", "inertia", "reciprocity"))
  )

  #Testing whether all statistics objects are from remstats class
  stats2 <- stats
  stats2[[15]] <- "a"
  expect_error(
    remx(reh = rehobj,
         statistics = stats2,
         random = c("baseline", "inertia", "reciprocity"))
  )

  #If model is tie and not all objects are tomstats
  rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "actor"))
  sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
  receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
  stats3 <- lapply(1:length(edgelist), function(x) remstats::aomstats(rehObj[[x]], sender_effects, receiver_effects))
  stats2[[15]] <- stats3[[15]]
  expect_error(
    remx(reh = rehObj,
         statistics = stats2,
         random = c("baseline", "inertia", "reciprocity"))
  )

  #If model is actor and not all objects are aomstats
  stats2 <- stats3
  stats2[[15]] <- stats[[15]]
  expect_error(
    remx(reh = rehObj,
         statistics = stats,
         random_sender = c("baseline"),
         random_receiver = c("indegreeReceiver", "rrankSend"))
  )
})
