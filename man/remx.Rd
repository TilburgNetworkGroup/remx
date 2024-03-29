% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remx.R
\name{remx}
\alias{remx}
\title{remx}
\usage{
remx(
  reh,
  statistics,
  random = NULL,
  random_sender = NULL,
  random_receiver = NULL,
  Niter = 10000,
  Nchain = 2,
  NburnIn = 5000,
  Nthin = 10
)
}
\arguments{
\item{reh}{a list containing multiple relational event sequences processed by the remify function}

\item{statistics}{a list containing multiple \code{remstats} objects}

\item{random}{a vector containing the names of the covariates that will be treated as random effects, use it for \code{model = "tie"}. If \code{random = NULL} all covariates will be treated as fixed-effects.}

\item{random_sender}{a vector containing the names of the covariates that will be treated as random effects in the sender model, use it for \code{model = "actor"}. If \code{random_sender = NULL} all covariates in the sender model will be treated as fixed-effects.}

\item{random_receiver}{a vector containing the names of the covariates that will be treated as random effects in the receiver model, use it for \code{model = "actor"}. If \code{random_receiver = NULL} all covariates in the receiver model will be treated as fixed-effects.}

\item{Niter}{number of iterations of the Gibbs sampler. The default is \code{Niter = 100000}.}

\item{Nchain}{number of chains of the Gibbs sampler.The default is \code{Nchain = 2}.}

\item{NburnIn}{number of samples to be discarded.The default is \code{NburnIn = 5000}.}

\item{Nthin}{number of samples to discarded in between every sample kept.The default is \code{NburnIn = 10}.}
}
\value{
Returns a \code{metarem} S3 object

For tie-oriented model:

\itemize{
\item \code{delta} a list containing arrays of random effects. Each element is a array with dimensions corresponding to samples, effects, and networks, respectively
\item \code{mu} fixed- and random-effect means, an array with dimensions corresponding to samples, effects and chains
\item \code{sigma} random-effect covariance matrix
\item \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
\item \code{MLE} A list containing the MLE estimates and their respective covariance matrices for each network
\item \code{random} If \code{TRUE}, that means the model contain random effects
\item \code{Niter} Number of iterations ran in the MCMC
\item \code{NburnIn} Number of samples discarded at the beginning of each chain
\item \code{Nchain} Number of chains
\item \code{Nthin} number of samples discarded in between each remaining sample
\item \code{run_time} model run time
}

For actor-oriented model:

\itemize{
\item \itemize{
Inside the sender_model list:
\item \code{gamma} a list containing random effects in the sender model. Each array has dimensions samples, effects and networks
\item \code{mu} fixed- and random-effect means, the array has dimensions samples, effects and chains
\item \code{sigma} random-effect covariance matrix
\item \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
}
\item \itemize{
In the receiver_model list:
\item \code{beta} a list containing arrays of random effects. Each element is a array with dimensions corresponding to samples, effects, and networks, respectively
\item \code{mu} fixed- and random-effect means, an array with dimensions corresponding to samples, effects and chains
\item \code{sigma} random-effect covariance matrix
\item \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
}
\item \code{MLE} A list containing the MLE estimates and their respective covariance matrices for each network
\item \code{random_sender} If \code{TRUE}, that means the sender model contain random effects
\item \code{random_receiver} If \code{TRUE}, that means the receiver model contain random effects
\item \code{Niter} Number of iterations ran in the MCMC
\item \code{NburnIn} Number of samples discarded at the beginning of each chain
\item \code{Nchain} Number of chains
\item \code{Nthin} number of samples discarded in between each remaining sample
\item \code{run_time} model run time
}
}
\description{
A function to fit mixed-effect relational event models. It runs meta-analytic approximations for actor- and tie-oriented relational event models. This function also supports full multilevel models from \code{lme4} and \code{rstanarm}.
}
\examples{
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
           random = c("baseline", "inertia"))

print(fit)

#----------------------------#
#    Actor-Oriented model    #
#----------------------------#

#Computing statistics
rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "actor"))
sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
stats <- lapply(1:length(edgelist), function(x) remstats::aomstats(rehObj[[x]],
sender_effects, receiver_effects))

#Running the model
fit <- remx(reh = rehObj,
            statistics = stats,
            random_sender = c("baseline", "outdegreeSender", "indegreeSender"),
            random_receiver = c("indegreeReceiver", "rrankSend"))

print(fit)

}
