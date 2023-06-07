## code to prepare `DATASET` dataset goes here
set.seed(100)
library(remulate)

networks <- vector("list")

#random effects means
mu <- c(-7, 0.01, -0.5)
#random effects variances
sigma <- rep(1, 3)

#Number of networks
K <- 15

for(i in 1:K){
  effects <- ~ baseline(rnorm(1,mu[1],sigma[1])) +
    inertia(rnorm(1,mu[2],sigma[2]), scaling = "std") +
    reciprocity((rnorm(1,mu[3],sigma[3])),scaling = "std")
  dat <- remulateTie(effects, actors = 1:10, time = Inf, events = 100, initial = 200)
  networks[[i]] <- dat$edgelist  
}

usethis::use_data(networks, compress = "xz")
