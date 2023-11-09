#' remx
#'
#' A function to fit mixed-effect relational event models. It runs meta-analytic approximations for actor- and tie-oriented relational event models. This function also supports full multilevel models from \code{lme4} and \code{rstanarm}.
#'
#' @param edgelist a list containing multiple relational event sequences
#' @param statistics a list containing multiple \code{remstats} objects
#' @param random a vector containing the names of the covariates that will be treated as random effects, use it for \code{model = "tie"}. If \code{random = NULL} all covariates will be treated as fixed-effects.
#' @param random_sender a vector containing the names of the covariates that will be treated as random effects in the sender model, use it for \code{model = "actor"}. If \code{random_sender = NULL} all covariates in the sender model will be treated as fixed-effects.
#' @param random_receiver a vector containing the names of the covariates that will be treated as random effects in the receiver model, use it for \code{model = "actor"}. If \code{random_receiver = NULL} all covariates in the receiver model will be treated as fixed-effects.
#' @param fixed a character vector containing the names of the covariables that will be treated as fixed. Only use this when \code{method = "classic"} or \code{method = "bayes"}.
#' @param intercept logical value indicating whether an intercept will be included in the model. The default is \code{intercept = TRUE}.
#' @param directed logical value indicating whether the model is for directed network. The default is \code{directed = TRUE}.
#' @param ordinal logical value indicating whether the timestamps of the events are know. The default is \code{ordinal = FALSE}. If \code{ordinal = TRUE}, that means that only the order of the events is known.
#' @param method the estimation method to be used. Choose \code{method = "meta"} for meta-analytic approximations, \code{method - "classic"} for frequentist and \code{method = "bayes"} for Bayesian.
#' @param model the type of model. Choose \code{model = "tie"} for tie oriented model and \code{model = "actor"} for the DyNaM.
#' @param Niter number of iterations of the Gibbs sampler. The default is \code{Niter = 100000}. This parameter is only used for \code{method = "meta"}.
#' @param Nchain number of chains of the Gibbs sampler.The default is \code{Nchain = 2}. This parameter is only used for \code{method = "meta"}.
#' @param NburnIn number of samples to be discarded.The default is \code{NburnIn = 5000}. This parameter is only used for \code{method = "meta"}.
#' @param Nthin number of samples to discarded in between every sample kept.The default is \code{NburnIn = 10}. This parameter is only used for \code{method = "meta"}.
#' @param actors a list containing the names of the actors in every network, should be the same lenght as \code{edgelist} and \code{statistics}
#' @param ... additional parameters. See the help of \code{lme4} for \code{method = "classic"} and \code{rstanarm} for \code{method = "bayes"}.
#'
#' @return  for \code{method = "meta"} returns a metarem S3 object
#' @return For tie-oriented model:
#' @return \code{delta} a list containing arrays of random effects. Each element is a array with dimensions corresponding to samples, effects, and networks, respectively
#' @return \code{mu} fixed- and random-effect means, an array with dimensions corresponding to samples, effects and chains
#' @return \code{sigma} random-effect covariance matrix
#' @return \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
#' @return \code{omega} scale-mixture parameter that allows beta to have a student T likelihood
#' @return \code{MLE} A list containing the MLE estimates and their respective covariance matrices for each network
#' @return \code{random} If \code{TRUE}, that means the model contain random effects
#' @return \code{Niter} Number of iterations ran in the MCMC
#' @return \code{NburnIn} Number of samples discarded at the beginning of each chain
#' @return \code{Nchain} Number of chains
#' @return \code{Nthin} number of samples discarded in between each remaining sample
#' @return \code{run_time} model run time
#'
#' @return For actor-oriented model:
#' @return Inside the sender_model list:
#' @return \code{gamma} a list containing random effects in the sender model. Each array has dimensions samples, effects and networks
#' @return \code{mu} fixed- and random-effect means, the array has dimensions samples, effects and chains
#' @return \code{sigma} random-effect covariance matrix
#' @return \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
#' @return \code{omega} scale-mixture parameter that allows beta to have a student T likelihood
#' @return In the receiver_model list:
#' @return \code{beta} a list containing arrays of random effects. Each element is a array with dimensions corresponding to samples, effects, and networks, respectively
#' @return \code{mu} fixed- and random-effect means, an array with dimensions corresponding to samples, effects and chains
#' @return \code{sigma} random-effect covariance matrix
#' @return \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
#' @return \code{omega} scale-mixture parameter that allows beta to have a student T likelihood
#' @return \code{MLE} A list containing the MLE estimates and their respective covariance matrices for each network
#' @return \code{random_sender} If \code{TRUE}, that means the sender model contain random effects
#' @return \code{random_receiver} If \code{TRUE}, that means the receiver model contain random effects
#' @return \code{Niter} Number of iterations ran in the MCMC
#' @return \code{NburnIn} Number of samples discarded at the beginning of each chain
#' @return \code{Nchain} Number of chains
#' @return \code{Nthin} number of samples discarded in between each remaining sample
#' @return \code{run_time} model run time
#'
#' @examples
#' #----------------------------#
#' #     Tie-Oriented model     #
#' #----------------------------#
#'
#' #Loading libraries
#' library(remx)
#' library(remstats)
#'
#' #Loading the data
#' edgelist <- networks
#' for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}
#'
#' #Computing statistics
#' effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "tie"))
#' stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, rehObj[[x]]))
#'
#' #Running the model
#'
#'fit <- remx(edgelist = edgelist,
#'            statistics = stats,
#'            random = c("baseline", "inertia"),
#'            method = "meta",
#'            model = "tie")
#'
#' print(fit)
#'
#' #----------------------------#
#' #    Actor-Oriented model    #
#' #----------------------------#
#'
#' #Computing statistics
#' rehObj <- lapply(1:length(edgelist), function(x) remify::remify(edgelist[[x]], model = "actor"))
#' sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
#' receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
#' stats <- lapply(1:length(edgelist), function(x) remstats::aomstats(rehObj[[x]],
#' sender_effects, receiver_effects))
#'
#' #Running the model
#' fit <- remx(edgelist = edgelist,
#'             statistics = stats,
#'             random_sender = c("baseline", "outdegreeSender", "indegreeSender"),
#'             random_receiver = c("indegreeReceiver", "rrankSend"),
#'             method = "meta",
#'             model = "actor")
#'
#' print(fit)
#'
#' @export
remx <- function(edgelist,
                 statistics,
                 random = NULL,
                 fixed = NULL,
                 random_sender = NULL,
                 random_receiver = NULL,
                 intercept = TRUE,
                 directed = TRUE,
                 ordinal = FALSE,
                 method = c("meta", "classic", "bayes"),
                 model = c("tie", "actor"),
                 Niter = 10000,
                 Nchain = 2,
                 NburnIn = 5000,
                 Nthin = 10,
                 actors = NULL, ...){
  #Running a few test to check some important arguments
  if(!is.list(edgelist)){
    stop("The edgelist argument should be a list containing multiple relational event sequences.")
  } else {
    if(length(edgelist) == 1){
      stop("The argument edgelist should be a list containing multiple relational event sequences. Try remstimate::remstimate() instead.")
    }
  }
  if(!is.list(statistics)){
    stop("The statistics argument should be a list containing multiple remstats objects.")
  } else{
    if(length(edgelist) != length(statistics)){
     stop("edgelist and statistics must have the same length!")
    }
    if(sum(sapply(statistics, function(x) class(x)[2]) == "remstats") != length(statistics)){
      stop("statistics must be a list containing remstats objects.")
    }
    if((model == "actor") & (sum(sapply(statistics, function(x) class(x)[1]) == "aomstats") != length(statistics))){
      stop("For model = 'actor' all arrays should be of class aomstats")
    }
    if((model == "tie") & (sum(sapply(statistics, function(x) class(x)[1]) == "tomstats") != length(statistics))){
      stop("For model = 'tie' all arrays should be of class tomstats")
    }
  }
  if(model == "actor" & (method == "classic" | method == "bayes")){
    stop("For model = 'actor', only method = 'meta' is supported.")
  }

  if(!is.null(fixed) & (method == "meta")){
    stop("The parameter fixed can only be declared for methods 'classic' and 'bayes'.")
  }

  #If the user enters either classic or bayes, they'll be prompted to answer
  #a yes or no question
  if(method == "classic" | method == "bayes"){
    if(ordinal){
      return(warning("For ordinal = TRUE, only method = 'meta' is supported."))
    }
    proceed <- utils::askYesNo("We do not advise method = 'classic' or method = 'bayes', it is going to take an enormous amount of time to fit the model. Do you wish to proceed?")
    if(is.na(proceed)){proceed <- FALSE}

    if(proceed){
      cat("As you wish! We are proceeding. You will pay with your time!\nRemember: Time is the only thing you won't get back.")
      #First steps is creating the formulas for the lme4 or rstanarm
      #Creating the formula to start the function
      form <- createFormula(random, fixed, intercept)
      #Reshaping the data (we need the cpp function rehPoisson)
      df <- lapply(1:length(statistics), function(x) cbind(rehPoisson(statistics[[x]]), x))
      df <- do.call("rbind", df)
      names(df) <- c("dyad", "logtimediff", dimnames(statistics[[1]]$statistics)[[3]], "cluster")

      if(method == "classic"){
        #Fitting the full frequentist model
        t1 <- Sys.time()
        fit <- lme4::glmer(formula = form,
                           data = df,
                           family = stats::poisson(link = "log"), ...)
        t2 <- Sys.time() - t1
        cat(paste("\nThe model took", format(round(t2,2), units = "mins"), "of your life."))
        #return(fit)
      }

      if(method == "bayes"){
        #Fitting the full Bayesian model
        t1 <- Sys.time()
        fit <- rstanarm::stan_glmer(formula = form,
                                    data = df,
                                    family = stats::poisson(link = "log"), ...)
        (t2 <- Sys.time() - t1)
        cat(paste("\nThe model took", format(round(t2,2), units = "mins"), "of your life."))
        #return(fit)
      }
    } else {
      return(warning("The model didn't run."))
    }
  } else {
    #This runs the meta-analytic model
    if(model == "tie"){
      #statistics <- lapply(statistics, function(x) x$statistics)
      #This runs the dyadic relational event model
      fit <- hrem(edgelist = edgelist,
                  statistics = statistics,
                  random = random,
                  Niter = Niter,
                  Nchain = Nchain,
                  NburnIn = NburnIn,
                  Nthin = Nthin,
                  actors = actors,
                  directed = directed,
                  ordinal = ordinal)
      fit$model <- "tie"
      #cat(paste("The model took", format(round(fit$run_time,2), units = "mins"), "out of your life."))
      #return(fit)
    } else if(model == "actor"){
      #This runs the actor-oriented model
      #statistics <- lapply(statistics, function(x) x$statistics)
      #for(i in 1:length(statistics)){names(statistics[[i]]) <- c("rate", "choice")}
      fit <- hremActor(edgelist = edgelist,
                       statistics = statistics,
                       random_sender = random_sender,
                       random_receiver = random_receiver,
                       actors = actors,
                       directed = directed,
                       ordinal = ordinal,
                       Niter = Niter,
                       Nchain = Nchain,
                       NburnIn = NburnIn,
                       Nthin = Nthin)
      fit$model <- "actor"
      #cat(paste("The model took", format(round(fit$run_time,2), units = "mins"), "out of your life."))
      #return(fit)
    }
  }

  #So in this part we will define the "metarem" class, this code has the same structure as remstimate
  output <- NULL

  if(method %in% c("classic", "bayes")){
    #In this case we want to preserve the original class, from lme4 or rstanarm,
    #so the user can get all the functionalities from those packages
    return(fit)
  } else { #Now we define the "metarem" class
    output <- structure(fit, class = "metarem")
    return(output)
  }
}

#########################################################################################
#########################################################################################
#########################################################################################

#' hrem
#'
#' A function to fit tie-oriented mixed-effect relational event models using a meta-analytix approximation. This function is called inside the \code{remx} function, so the user should avoid calling this function and using \code{remx} instead.
#'
#' @param edgelist a list containing multiple relational event sequences
#' @param statistics a list containing multiple \code{remstats} objects
#' @param random a vector containing the names of the covariates that will be treated as random effects.
#' @param Niter number of iterations of the Gibbs sampler. The default is \code{Niter = 100000}.
#' @param Nchain number of chains of the Gibbs sampler.The default is \code{Nchain = 2}.
#' @param NburnIn number of samples to be discarded.The default is \code{NburnIn = 5000}.
#' @param Nthin number of samples to discarded in between every sample kept.The default is \code{NburnIn = 10}.
#' @param actors a list containing the names of the actors in every network, should be the same lenght as \code{edgelist} and \code{statistics}.
#' @param ordinal logical value indicating whether the timestamps of the events are know. The default is \code{ordinal = FALSE}. If \code{ordinal = TRUE}, that means that only the order of the events is known.
#' @param directed logical value indicating whether the model is for directed network. The default is \code{directed = TRUE}.
#'
#' @return  A list containing the output of the Gibbs sampler for this model.
hrem <- function(edgelist,
                 statistics,
                 random = NULL,
                 Niter = 10000,
                 Nchain = 2,
                 NburnIn = 5000,
                 Nthin = 10,
                 actors = NULL,
                 ordinal = FALSE,
                 directed = TRUE){

  #edgelist: this will be a list with multiple relational event sequences
  #statistics: This is list containing multiple arrays with statistics
  #reh and statistics must be of the same length

  #This look gets the estimates from the remstimate function and store them in a
  #list with the MLE's and their standard errors
  
  #Checking whether all covariables declared in the random parameter are in the data matrix
  covs <- dimnames(statistics[[1]])[[3]]
  cnms <- sum(!(random %in% covs))
  if(cnms > 0){
    not_in_covs <- which(!(random %in% covs))
    stop(paste("Covariables", random[not_in_covs], "might not be in the data matrix.", "Check parameter random for typos."))
  }
  
  estimates <- vector("list")
  t1 <- Sys.time()
  for(i in 1:length(edgelist)){

    edgelistREH <- remify::remify(edgelist[[i]], model = "tie", actors = actors[[i]],
                               directed = directed, ordinal = ordinal)

    rem <- remstimate::remstimate(reh = edgelistREH,
                                  stats = statistics[[i]],
                                  method = "MLE",
                                  #model = "tie",
    )

    estimates[[i]] <- list(coef = rem$coefficients, var = rem$vcov)

    print(paste("MLE for cluster", i, "has been obtained."))

  }
  #getting the names for the fixed-effects
  fixed <- dimnames(statistics[[1]])[3][[1]]

  names(estimates) <- paste0("cluster_", 1:length(edgelist))


  #checking whether all statistics array contain the same number of covariates
  K <- length(edgelist) #number of groups

  if(!is.null(random)){
    p <- length(random) #number random effects
    #We assume that all statistics array will contain the same statistics
    beta_hat <- as.matrix(sapply(estimates, function(x) (x$coef[random])))
    omega_hat <- lapply(estimates, function(x) as.matrix(x$var[random,random]))
    omega_hat <- simplify2array(omega_hat, except = NULL)

    #Decomposing the covariance matrix to compute the conditional distribution of the random effects
    omega_11 <- lapply(estimates, function(x) matrix(x$var[!(fixed %in% random),!(fixed %in% random)], nrow = sum(!(fixed %in% random)), ncol = sum(!(fixed %in% random))))
    omega_22 <- lapply(estimates, function(x) matrix(x$var[random,random], nrow = p, ncol = p))
    omega_21 <- lapply(estimates, function(x) matrix(x$var[!(fixed %in% random),(fixed %in% random)], ncol = sum(fixed %in% random), nrow = sum(!(fixed %in% random))))
    omega_11 <- simplify2array(omega_11, except = NULL)
    omega_22 <- simplify2array(omega_22, except = NULL)
    omega_21 <- simplify2array(omega_21, except = NULL)
    if(p == 1){
      beta_hat <- t(beta_hat)
    }
    randomModel <- TRUE
  } else {
    p <- 0
    beta_hat <- matrix(0, nrow = 1, ncol = 1)
    omega_hat <- array(0, dim = c(1,1,K))
    omega_11 <- array(0, dim = c(1,1,K))
    omega_22 <- array(0, dim = c(1,1,K))
    omega_21 <- array(0, dim = c(1,1,K))
    delta <- array(0, dim = c(1, 1, K))
    randomModel <- FALSE
  }

  q <- length(fixed) #number of fixed effects
  psi <- as.matrix(sapply(estimates, function(x) x$coef))
  SigmaPsi <- lapply(estimates, function(x) x$var)
  if(q == 1){
    SigmaPsi <- array(simplify2array(SigmaPsi), dim = c(q,q,K))
    psi <- t(psi)
  } else {
    SigmaPsi <- simplify2array(SigmaPsi)
  }

  if(p == 0){
    mu <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,1))) #random-effect means
    #I'm making a confusing, using lambda here, but on the paper, these are called alphas
    Lambda <- matrix(stats::rgamma(1*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    #initial values for MCMC
    delta <- array(NA, dim = c(1, K, Nchain)) #group-specific effects
    sigma <- array(NA, dim = c(1, 1, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      delta[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,1)))
      sigma[,,i] <- diag(stats::rpois(1,100), nrow = 1)
    }
  } else {
    mu <- t(mvtnorm::rmvnorm(Nchain, mean = rowMeans(beta_hat))) #random-effect means
    #I'm making a confusing, using lambda here, but on the paper, these are called alphas
    Lambda <- matrix(stats::rgamma(p*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    #initial values for MCMC
    delta <- array(NA, dim = c(p, K, Nchain)) #group-specific effects
    sigma <- array(NA, dim = c(p, p, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      delta[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,p)))
      sigma[,,i] <- diag(stats::rpois(p,100), nrow = p)
    }
  }

  #Prior for Mu and Sigma
  eta <- 2
  xi <- 10 #it can be an arbitrarily large number (Gelman, 2006)
  TauMean <- diag(10, q) #scale matrix of mu
  if(!is.null(random)){
    #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
    muMat <- matrix(0, nrow = p, ncol = q)
    ind <- which(fixed %in% random)
    for(i in 1:length(ind)){muMat[i,ind[i]] <- 1}

    if(p != q){
      #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
      psiMat <- matrix(0, nrow = q-p, ncol = q)
      ind <- which(!(fixed %in% random))
      for(i in 1:length(ind)){psiMat[i,ind[i]] <- 1}
      #initial values
    } else {
      #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
      psiMat <- muMat
    }
    #initial values
    psi_in <- t(mvtnorm::rmvnorm(Nchain, mean =  psiMat %*% rowMeans(psi))) #random-effect means

  } else {
    psi_in <- t(mvtnorm::rmvnorm(Nchain, mean = rowMeans(psi))) #random-effect means
  }

  print(paste("Running the Gibbs sampler."))

  MCMC <- sampler(Niter = Niter,
                  Nchain = Nchain,
                  NburnIn = NburnIn,
                  Nthin = Nthin,
                  p = p,
                  q = q,
                  K = K,
                  delta = delta,
                  beta_hat = beta_hat,
                  omega_hat = omega_hat,
                  omega_hat_11 = omega_11,
                  omega_hat_21 = omega_21,
                  omega_hat_22 = omega_22,
                  Mu = mu,
                  Sigma = sigma,
                  eta = eta,
                  Lambda = Lambda,
                  xi = xi,
                  SigmaMean = SigmaPsi,
                  psi = psi_in,
                  psi_hat = psi,
                  random_effect = muMat,
                  fixed_effect = psiMat,
                  randomModel = randomModel)

  if(!is.null(fixed)){
    colnames(MCMC$mu) <- fixed
  }
  if(!is.null(random)){
    #colnames(MCMC$mu) <- random
    colnames(MCMC$alpha) <- random
    MCMC$MLE <- estimates
    for(i in 1:Nchain){
      dimnames(MCMC$delta[[i]])[[2]] <- random
      dimnames(MCMC$sigma[[i]])[[1]] <- dimnames(MCMC$sigma[[i]])[[2]] <- random
      dimnames(MCMC$delta[[i]])[[3]] <- paste0("cluster_", 1:K)
    }
    names(MCMC$delta) <- paste0("chain_", 1:Nchain)
    names(MCMC$sigma) <- paste0("chain_", 1:Nchain)
  }
  if(is.null(random)){
    MCMC$random <- F
  } else {
    MCMC$random <- T
  }
  MCMC$Niter <- Niter
  MCMC$NburnIn <- NburnIn
  MCMC$Nchain <- Nchain
  MCMC$Nthin <- Nthin
  print(paste("The Gibbs sampler has finished."))
  t <- Sys.time() - t1

  MCMC$run_time <- t

  return(MCMC) #returns the list with posterior samples from sampler function

}

#########################################################################################
#########################################################################################
#########################################################################################
#' hremActor
#'
#' A function to fit DyNaM mixed-effect relational event models using a meta-analytic approximation. This function is called inside the \code{remx} function, so the user should avoid calling this function and using \code{remx} instead.
#'
#' @param edgelist a list containing multiple relational event sequences
#' @param statistics a list containing multiple \code{remstats} objects
#' @param random_sender a vector containing the names of the covariates that will be treated as random effects in the sender model.
#' @param random_receiver a vector containing the names of the covariates that will be treated as random effects in the receiver model.
#' @param Niter number of iterations of the Gibbs sampler. The default is \code{Niter = 100000}.
#' @param Nchain number of chains of the Gibbs sampler.The default is \code{Nchain = 2}.
#' @param NburnIn number of samples to be discarded.The default is \code{NburnIn = 5000}.
#' @param Nthin number of samples to discarded in between every sample kept.The default is \code{NburnIn = 10}.
#' @param actors a list containing the names of the actors in every network, should be the same lenght as \code{edgelist} and \code{statistics}.
#' @param ordinal logical value indicating whether the timestamps of the events are know. The default is \code{ordinal = FALSE}. If \code{ordinal = TRUE}, that means that only the order of the events is known.
#' @param directed logical value indicating whether the model is for directed network. The default is \code{directed = TRUE}.
#'
#' @return  A list containing the output of the Gibbs sampler for this model.
hremActor <- function(edgelist,
                      statistics,
                      Niter = 10000,
                      Nchain = 2,
                      NburnIn = 5000,
                      Nthin = 10,
                      random_sender = NULL,
                      random_receiver = NULL,
                      actors = NULL,
                      ordinal = FALSE,
                      directed = TRUE){
  #HERE THE FUNCTION STARTS
  
  #Checking whether all covariables declared in the random parameter are in the data matrix
  covs_sender <- dimnames(statistics[[1]]$sender_stats)[[3]]
  cnms_sender <- sum(!(random_sender %in% covs_sender))
  if(cnms_sender > 0){
    not_in_covs_sender <- which(!(random_sender %in% covs_sender))
    stop(paste("Covariables", random_sender[not_in_covs_sender], "might not be in the sender data matrix.", "Check parameter random for typos."))
  }
  
  covs_receiver <- dimnames(statistics[[1]]$receiver_stats)[[3]]
  cnms_receiver <- sum(!(random_receiver %in% covs_receiver))
  if(cnms_receiver > 0){
    not_in_covs_receiver <- which(!(random_receiver %in% covs_receiver))
    stop(paste("Covariables", random_receiver[not_in_covs_receiver], "might not be in the receiver data matrix.", "Check parameter random for typos."))
  }
  
  estimates <- vector("list")
  t1 <- Sys.time()
  for(i in 1:length(edgelist)){
    edgelistREH <- remify::remify(edgelist[[i]], model = "actor", actors = actors[[i]],
                                  directed = directed, ordinal = ordinal)
    
    rem <- remstimate::remstimate(reh = edgelistREH,
                                  stats = statistics[[i]],
                                  method = "MLE",
                                  #model = "actor",
    )
    
    estimates[[i]] <- list(sender = list(coef = rem$sender_model$coefficients, var = rem$sender_model$vcov),
                           receiver = list(coef = rem$receiver_model$coefficients, var = rem$receiver_model$vcov))
    
    print(paste("MLE for cluster", i, "has been obtained."))
    
  }
  
  fixed_sender <- dimnames(statistics[[1]]$sender_stats)[3][[1]]
  fixed_receiver <- dimnames(statistics[[1]]$receiver_stats)[3][[1]]
  
  #checking whether all statistics array contain the same number of covariates
  K <- length(edgelist) #number of groups
  
  #Checking for random effects
  #SENDER MODEL
  if(!is.null(random_sender)){
    p <- length(random_sender) #number random effects sender model
    #Sender model
    gamma_hat <- as.matrix(sapply(estimates, function(x) (x$sender$coef[random_sender])))
    zeta_hat <- lapply(estimates, function(x) as.matrix(x$sender$var[random_sender,random_sender]))
    zeta_hat <- simplify2array(zeta_hat, except = NULL)
    
    #Decomposing the covariance matrix to compute the conditional distribution of the random effects
    zeta_11 <- lapply(estimates, function(x) matrix(x$sender$var[!(fixed_sender %in% random_sender),!(fixed_sender %in% random_sender)], ncol = sum(!(fixed_sender %in% random_sender)), nrow = sum(!(fixed_sender %in% random_sender))))
    zeta_22 <- lapply(estimates, function(x) matrix(x$sender$var[random_sender,random_sender], nrow = p, ncol = p))
    zeta_21 <- lapply(estimates, function(x) matrix(x$sender$var[!(fixed_sender %in% random_sender),(fixed_sender %in% random_sender)], ncol = sum(fixed_sender %in% random_sender), nrow = sum(!(fixed_sender %in% random_sender))))
    zeta_11 <- simplify2array(zeta_11, except = NULL)
    zeta_22 <- simplify2array(zeta_22, except = NULL)
    zeta_21 <- simplify2array(zeta_21, except = NULL)
    #Fixing the structure of the random effects in the sender model
    if(p == 1){
      gamma_hat <- t(gamma_hat)
    }
    randomModelSnd <- TRUE
  } else {
    p <- 0
    gamma_hat <- matrix(0, nrow = 1, ncol = 1)
    zeta_hat <- array(0, dim = c(1,1,K))
    zeta_11 <- array(0, dim = c(1,1,K))
    zeta_22 <- array(0, dim = c(1,1,K))
    zeta_21 <- array(0, dim = c(1,1,K))
    gamma <- array(0, dim = c(1, 1, K))
    randomModelSnd <- FALSE
  }
  
  #RECEIVER MODEL
  if(!is.null(random_receiver)){
    v <- length(random_receiver) #random effects receiver model
    #Receiver model
    beta_hat <- as.matrix(sapply(estimates, function(x) (x$receiver$coef[random_receiver])))
    omega_hat <- lapply(estimates, function(x) as.matrix(x$receiver$var[random_receiver,random_receiver]))
    omega_hat <- simplify2array(omega_hat, except = NULL)
    
    #Decomposing the covariance matrix to compute the conditional distribution of the random effects
    omega_11 <- lapply(estimates, function(x) matrix(x$receiver$var[!(fixed_receiver %in% random_receiver),!(fixed_receiver %in% random_receiver)], ncol = sum(!(fixed_receiver %in% random_receiver))))
    omega_22 <- lapply(estimates, function(x) matrix(x$receiver$var[random_receiver,random_receiver], nrow = v, ncol = v))
    omega_21 <- lapply(estimates, function(x) matrix(x$receiver$var[!(fixed_receiver %in% random_receiver),(fixed_receiver %in% random_receiver)], ncol = sum(fixed_receiver %in% random_receiver), nrow = sum(!(fixed_receiver %in% random_receiver))))
    omega_11 <- simplify2array(omega_11, except = NULL)
    omega_22 <- simplify2array(omega_22, except = NULL)
    omega_21 <- simplify2array(omega_21, except = NULL)
    #Fixing the structure of the random effects in the receiver model
    if(v == 1){
      beta_hat <- t(beta_hat)
    }
    randomModelRec <- TRUE
  } else {
    v <- 0
    beta_hat <- matrix(0, nrow = 1, ncol = 1)
    omega_hat <- array(0, dim = c(1,1,K))
    omega_11 <- array(0, dim = c(1,1,K))
    omega_22 <- array(0, dim = c(1,1,K))
    omega_21 <- array(0, dim = c(1,1,K))
    beta <- array(0, dim = c(1, 1, K))
    randomModelRec <- FALSE
  }
  
  
  #Checking if the receiver model contains fixed effects
  if(!is.null(fixed_receiver)){
    u <- length(fixed_receiver) #number of fixed effects
    psi <- as.matrix(sapply(estimates, function(x) x$receiver$coef))
    SigmaPsi <- lapply(estimates, function(x) as.matrix(x$receiver$var))
    SigmaPsi <-  simplify2array(SigmaPsi, except = NULL)
    if(u == 1){
      psi <- t(psi)
    }
  } else {
    psi <- matrix(0, nrow = 1, ncol = 1)
    SigmaPsi <- matrix(0, nrow = 1, ncol = 1)
    u <- 0
  }
  
  #Checking if the sender model contains fixed effects
  if(!is.null(fixed_sender)){
    q <- length(fixed_sender) #number of fixed effects
    phi <- as.matrix(sapply(estimates, function(x) x$sender$coef))
    SigmaPhi <- lapply(estimates, function(x) as.matrix(x$sender$var))
    SigmaPhi <- simplify2array(SigmaPhi, except = NULL)
    if(q == 1){
      phi <- t(phi)
    }
  } else {
    phi <- matrix(0, nrow = 1, ncol = 1)
    SigmaPhi <- matrix(0, nrow = 1, ncol = 1)
    q <- 0
  }
  
  ###########################################################################################
  ###########################################################################################
  ###########################################################################################
  
  #Setting up the additional quantities needed
  
  if(v == 0){
    #initial values for MCMC
    Lambda_beta <- matrix(stats::rgamma(1*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    mu_beta <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,1))) #random-effect means
    beta <- array(NA, dim = c(1, K, Nchain)) #group-specific effects
    sigma <- array(NA, dim = c(1, 1, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      beta[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,1)))
      sigma[,,i] <- diag(stats::rpois(1,100), nrow = 1)
    }
  } else {
    #initial values for MCMC
    Lambda_beta <- matrix(stats::rgamma(v*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    mu_beta <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,v))) #random-effect means
    beta <- array(NA, dim = c(v, K, Nchain)) #group-specific effects
    sigma <- array(NA, dim = c(v, v, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      beta[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,v)))
      sigma[,,i] <- diag(stats::rpois(v,100), nrow = v)
    }
  }
  if(p == 0){
    #I'm making a confusing, using lambda here, but on the paper, these are called alphas
    Lambda_gamma <- matrix(stats::rgamma(1*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    mu_gamma <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,1))) #random-effect means
    gamma <- array(NA, dim = c(1, K, Nchain)) #group-specific effects
    zeta <- array(NA, dim = c(1, 1, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      gamma[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,1)))
      zeta[,,i] <- diag(stats::rpois(1,100), nrow = 1)
    }
  } else {
    Lambda_gamma <- matrix(stats::rgamma(p*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    mu_gamma <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,p))) #random-effect means
    gamma <- array(NA, dim = c(p, K, Nchain)) #group-specific effects
    zeta <- array(NA, dim = c(p, p, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      gamma[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,p)))
      zeta[,,i] <- diag(stats::rpois(p,100), nrow = p)
    }
  }
  
  #Prior for Mu and Sigma
  eta <- 2
  xi <- 10 #it can be an arbitrarily large number (Gelman, 2006)
  
  #Building the matrices that will separate random and fixed effects for receiver model
  if(!is.null(random_sender)){
    #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
    muGMat <- matrix(0, nrow = p, ncol = q)
    ind <- which(fixed_sender %in% random_sender)
    for(i in 1:length(ind)){muGMat[i,ind[i]] <- 1}
    
    if(p != q){
      #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
      phiMat <- matrix(0, nrow = q-p, ncol = q)
      ind <- which(!(fixed_sender %in% random_sender))
      for(i in 1:length(ind)){phiMat[i,ind[i]] <- 1}
      #initial values
    } else {
      #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
      phiMat <- muGMat
    }
    
    #initial values
    phi_in <- t(mvtnorm::rmvnorm(Nchain, mean =  phiMat %*% as.matrix(rowMeans(phi)))) #random-effect means
    
  } else {
    phi_in <- t(mvtnorm::rmvnorm(Nchain, mean = rowMeans(phi))) #random-effect means
  }
  
  #Building the matrices that will separate random and fixed effects for receiver model
  if(!is.null(random_receiver)){
    #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
    muBMat <- matrix(0, nrow = v, ncol = u)
    ind <- which(fixed_receiver %in% random_receiver)
    for(i in 1:length(ind)){muBMat[i,ind[i]] <- 1}
    
    if(u != v){
      #This is a matrix used to separate the random and the fixed effects, in the conditional distribution of delta
      psiMat <- matrix(0, nrow = u-v, ncol = u)
      ind <- which(!(fixed_receiver %in% random_receiver))
      for(i in 1:length(ind)){psiMat[i,ind[i]] <- 1}
    } else {
      psiMat <- muBMat
    }
    #initial values
    psi_in <- t(mvtnorm::rmvnorm(Nchain, mean =  psiMat %*% as.matrix(rowMeans(psi)))) #random-effect means
    
  } else {
    psi_in <- t(mvtnorm::rmvnorm(Nchain, mean = rowMeans(psi))) #random-effect means
  }
  
  print(paste("Running the Gibbs sampler."))
  
  MCMC <- samplerActor(Niter = Niter,
                       Nchain = Nchain,
                       NburnIn = NburnIn,
                       Nthin = Nthin,
                       v = v,
                       u = u,
                       K = K,
                       beta = beta,
                       beta_hat = beta_hat,
                       omega_hat = omega_hat,
                       omega_hat_11 = omega_11,
                       omega_hat_21 = omega_21,
                       omega_hat_22 = omega_22,
                       MuB = mu_beta,
                       Sigma = sigma,
                       eta = eta,
                       LambdaB =  Lambda_beta,
                       xi = xi,
                       SigmaPsi = SigmaPsi,
                       psi = psi_in,
                       psi_hat = psi,
                       p = p,
                       q = q,
                       gamma = gamma,
                       gamma_hat = gamma_hat,
                       zeta_hat = zeta_hat,
                       zeta_hat_11 = zeta_11,
                       zeta_hat_21 = zeta_21,
                       zeta_hat_22 = zeta_22,
                       MuG = mu_gamma,
                       Zeta = zeta,
                       LambdaG = Lambda_gamma,
                       SigmaPhi = SigmaPhi,
                       phi = phi_in,
                       phi_hat = phi,
                       fixed_effect_rec = psiMat,
                       random_effect_rec = muBMat,
                       fixed_effect_snd = phiMat,
                       random_effect_snd = muGMat,
                       randomModelSnd = randomModelSnd,
                       randomModelRec = randomModelRec)
  
  print(paste("The Gibbs sampler has finished."))
  
  #Receiver model
  colnames(MCMC$receiver_model$mu) <- fixed_receiver
  if(!is.null(random_receiver)){
    colnames(MCMC$receiver_model$alpha) <- random_receiver
    MCMC$MLE <- estimates
    for(i in 1:Nchain){
      colnames(MCMC$receiver_model$beta[[i]]) <- random_receiver
      colnames(MCMC$receiver_model$sigma[[i]]) <- rownames(MCMC$receiver_model$sigma[[i]]) <- random_receiver
      dimnames(MCMC$receiver_model$beta[[i]])[[3]] <- paste0("cluster_", 1:K)
    }
    names(MCMC$receiver_model$beta) <- paste0("chain_", 1:Nchain)
    names(MCMC$receiver_model$sigma) <- paste0("chain_", 1:Nchain)
  }
  
  #Sender model
  colnames(MCMC$sender_model$mu) <- fixed_sender
  if(!is.null(random_sender)){
    colnames(MCMC$sender_model$alpha) <- random_sender
    MCMC$MLE <- estimates
    for(i in 1:Nchain){
      dimnames(MCMC$sender_model$gamma[[i]])[[2]] <- random_sender
      colnames(MCMC$sender_model$sigma[[i]]) <- rownames(MCMC$sender_model$sigma[[i]]) <- random_sender
      dimnames(MCMC$sender_model$gamma[[i]])[[3]] <- paste0("cluster_", 1:K)
    }
    names(MCMC$sender_model$gamma) <- paste0("chain_", 1:Nchain)
    names(MCMC$sender_model$sigma) <- paste0("chain_", 1:Nchain)
  }
  if(is.null(random_sender)){
    MCMC$random_sender <- F
  } else {
    MCMC$random_sender <- T
  }
  if(is.null(random_receiver)){
    MCMC$random_receiver <- F
  } else {
    MCMC$random_receiver <- T
  }
  MCMC$Niter <- Niter
  MCMC$NburnIn <- NburnIn
  MCMC$Nchain <- Nchain
  MCMC$Nthin <- Nthin
  MCMC$run_time <- Sys.time() - t1
  return(MCMC)
}


#########################################################################################
#########################################################################################
#########################################################################################

#This function will create a formula from the names of the variables that are included in the
#fixed and random parameters of the remx function

#' createFormula
#'
#' A function only used inside the \code{remx} function to create a formula that will be passed to \code{glmer} or \code{stan_glmer} functions. It should never be called outside \code{remx}
#'
#' @param random a vector containing the names of variables to be treated as random effects
#' @param fixed a vector containing the names of variables to be treated as fixed effects
#' @param intercept a logical value stating whether the formula should contain an intercept or not, the default is \code{intercept = TRUE}
#'
#' @return a formula object
createFormula <- function(random = NULL, fixed = NULL, intercept = T){

  if(!is.null(random)){
    rnd <- paste0(random, collapse = "+")
  }
  if(!is.null(fixed)){
    #We always want to have the random as fixed effect too!
    #Since the distribution of the random effects is N(0, Sigma)
    fix <- paste0(c(fixed, random), collapse = "+")
  }

  #Creating the formula
  if(intercept){
    if(is.null(random)){
      if("baseline" %in% fixed){
        form <- stats::as.formula(paste0("dyad", "~", fix, "-1+offset(logtimediff)"))
      } else {
        form <- stats::as.formula(paste0("dyad", "~", fix, "+offset(logtimediff)"))
      }
    } else if(is.null(fixed)){
      if("baseline" %in% random){
        form <- stats::as.formula(paste0("dyad", "~", "offset(logtimediff)+", "(", rnd, "-1|cluster)"))
      } else {
        form <- stats::as.formula(paste0("dyad", "~", "offset(logtimediff)+", "(", rnd, "|cluster)"))
      }
    } else {
      if("baseline" %in% random){
        form <- stats::as.formula(paste0("dyad", "~", fix, "-1+offset(logtimediff)+", "(", rnd, "-1|cluster)"))
      } else {
        form <- stats::as.formula(paste0("dyad", "~", fix, "+offset(logtimediff)+", "(", rnd, "|cluster)"))
      }
    }
  } else {
    if(is.null(random)){
      form <- stats::as.formula(paste0("dyad", "~", fix, "-1+offset(logtimediff)"))
    } else if(is.null(fixed)){
      form <- stats::as.formula(paste0("dyad", "~ -1+offset(logtimediff)+", "(", rnd, "-1|cluster)"))
    } else {
      form <- stats::as.formula(paste0("dyad", "~", fix, "-1+offset(logtimediff)+", "(", rnd, "-1|cluster)"))
    }
  }

  return(form)

}

#########################################################################################
#########################################################################################
#########################################################################################

#From here on, we have the auxiliary functions that extract some useful quantities from
#the models, there are functions similar to this ones in rstanarm and lme4

#####################################################################################
#####################################################################################
#####################################################################################

#' fixef_metarem
#'
#' A function to extract fixed effects of a \code{metarem} object. Similar to the function available in \code{lme4} and \code{rstanarm}.
#'
#' @param object a metarem object
#'
#' @return a list containing the posterior means of the fixed effects and their standard deviations
#' @export
#This function extracts the fixed effects from the metarem object
fixef_metarem <- function(object) {
  #x must be a metarem object
  if(!inherits(object, "metarem")){
    stop("object must be an object of metarem class.")
  } else {
    #Here I have a colMeans inside a rowMeans, just in case we have
    #multiple chains in the model

    #The reason for the following lines is that the mean of the random effects
    #should also be seen as a fixed effect
    if(object$model == "tie"){
      a <- object$mu
      fe <- rowMeans(colMeans(a))
      std <- rowMeans(apply(a, c(2,3), stats::sd))
    } else { #model == "actor
      a <- abind::abind(object$sender_model$mu, object$sender_model$phi, along = 2)
      snd <- rowMeans(colMeans(a))
      stdA <- rowMeans(apply(a, c(2,3), stats::sd))
      b <- object$receiver_model$mu
      rec <- rowMeans(colMeans(b))
      stdB <- rowMeans(apply(b, c(2,3), stats::sd))
      fe <- list(sender_model = snd,
                 receiver_model = rec)
      std <- list(sender_model = stdA,
                  receiver_model = stdB)
    }
    return(list(fixed_effects = fe,
                std_error = std))
  }
}

#####################################################################################
#####################################################################################
#####################################################################################

#' ranef_metarem
#'
#' A function to extract random effects of a \code{metarem} object. Similar to the function available in \code{lme4} and \code{rstanarm}.
#'
#' @param object a metarem object
#'
#' @return a list containing the posterior means of the random effects and their standard deviations
#' @export
#This function extracts the random effects from the metarem object
ranef_metarem <- function(object) {
  if(!inherits(object, "metarem")){
    stop("object must be an object of metarem class.")
  } else {
    #Here I have a colMeans inside a rowMeans, just in case we have
    #multiple chains in the model
    if(object$model == "tie"){
      if(is.na(object$random)){
        stop("Model does not contain random effects.")
      } else {
        variable <- dimnames(object$delta[[1]])[[2]]
        number_rnd_eff <- dim(object$delta[[1]])[2]
        if(number_rnd_eff == 1){
          rnd <- matrix(rowMeans(sapply(object$delta, function(y) apply(y, 3, colMeans))), ncol = 1)
          std <- matrix(rowMeans(sapply(object$delta, function(y) apply(y, 3, stats::sd))), ncol = 1)
          colnames(rnd) <- colnames(std) <- variable
        } else {
          if(is.na(object$random)){
            stop("Model does not contain random effects.")
          } else{
            a <- lapply(object$delta, function(y) t(colMeans(y)))
            a <- simplify2array(a)
            rnd <- t(apply(a, 1, rowMeans))
            std <- lapply(object$delta, function(y) t(apply(y, c(2,3), stats::sd)))
            std <- simplify2array(std)
            std <- t(apply(std, 1, rowMeans))
          }
        }
      }
    } else { #model == "actor"
      if(is.na(object$random_sender) & is.na(object$random_receiver)){
        stop("Model does not contain random effects.")
      } else {
        if(!object$random_sender){
          rnd_snd <- NULL
          std_snd <- NULL
        } else {
          variable <- dimnames(object$sender_model$gamma[[1]])[[2]]
          number_rnd_snd <- dim(object$sender_model$gamma[[1]])[2]
          if(number_rnd_snd == 1){
            rnd_snd <- matrix(rowMeans(sapply(object$sender_model$gamma, function(y) apply(y, 3, colMeans))), ncol = 1)
            std_snd <- matrix(rowMeans(sapply(object$sender_model$gamma, function(y) apply(y, 3, stats::sd))), ncol = 1)
            colnames(rnd_snd) <- colnames(std_snd) <- variable
          } else {
            #here comes the case of multiple random effects
            #here comes the case of multiple random effects
            a <- lapply(object$sender_model$gamma, function(y) t(colMeans(y)))
            a <- simplify2array(a)
            rnd_snd <- t(apply(a, 1, rowMeans))
            std <- lapply(object$sender_model$gamma, function(y) t(apply(y, c(2,3), stats::sd)))
            std <- simplify2array(std)
            std_snd <- t(apply(std, 1, rowMeans))
          }
        }
        if(!object$random_receiver){
          rnd_rec <- NULL
          std_rec <- NULL
        } else {
          variable <- dimnames(object$receiver_model$beta[[1]])[[2]]
          number_rnd_rec <- dim(object$receiver_model$beta[[1]])[2]
          if(number_rnd_rec == 1){
            rnd_rec <- matrix(rowMeans(sapply(object$receiver_model$beta, function(y) apply(y, 3, colMeans))), ncol = 1)
            std_rec <- matrix(rowMeans(sapply(object$receiver_model$beta, function(y) apply(y, 3, stats::sd))), ncol = 1)
            colnames(rnd_rec) <- colnames(std_rec) <- variable
          } else {
            #here comes the case of multiple random effects
            a <- lapply(object$receiver_model$beta, function(y) t(colMeans(y)))
            a <- simplify2array(a)
            rnd_rec <- t(apply(a, 1, rowMeans))
            std <- lapply(object$receiver_model$beta, function(y) t(apply(y, c(2,3), stats::sd)))
            std <- simplify2array(std)
            std_rec<- t(apply(std, 1, rowMeans))
          }
        }
      }
      rnd <- list(sender_model = rnd_snd,
                  receiver_model = rnd_rec)
      std <- list(sender_model = std_snd,
                  receiver_model = std_rec)
    }

    return(list(random_effects = rnd,
                std_error = std))
  }
}

#####################################################################################
#####################################################################################
#####################################################################################

#' VarCov_metarem
#'
#' A function to extract the covariance matrix of the random effects of a \code{metarem} object.
#'
#' @param object a metarem object
#'
#' @return a matrix of the posterior means of the covariance matrix of the random effects.
#' @export
VarCov_metarem <- function(object){
  if(!inherits(object, "metarem")){
    stop("object must be an object of metarem class.")
  } else {
    if(object$model == "tie"){
      if(object$random){
        sigma <- lapply(object$sigma, function(x) rowMeans(x, dims = 2))
        sigma <- simplify2array(sigma)
        if(is.null(dim(sigma))){
          VarCov <- as.matrix(mean(sigma))
          nome <- dimnames(object$sigma[[1]])[[1]]
          rownames(VarCov) <- colnames(VarCov) <- nome
        } else {
          VarCov <- rowMeans(sigma, dims = 2)
          VarCov <- (VarCov + t(VarCov))/2 #forcing it to be symmetric
        }
      } else {
        warning("Model does not contain random effects.")
      }
    } else { #model = "actor"
      if(object$random_sender){
        sigma <- lapply(object$sender_model$sigma, function(x) rowMeans(x, dims = 2))
        sigma <- simplify2array(sigma)
        if(dim(object$sender_model$sigma[[1]])[1] == 1){
          VarCovSnd <- as.matrix(mean(sigma))
          nome <- dimnames(object$sender_model$sigma[[1]])[[1]]
          rownames(VarCovSnd) <- colnames(VarCovSnd) <- nome
        } else {
          VarCovSnd <- rowMeans(sigma, dims = 2)
          VarCovSnd <- (VarCovSnd + t(VarCovSnd))/2
        }
      } else {
        VarCovSnd <- NULL
      }
      if(object$random_receiver){
        sigma <- lapply(object$receiver_model$sigma, function(x) rowMeans(x, dims = 2))
        sigma <- simplify2array(sigma)
        if(dim(object$receiver_model$sigma[[1]])[1] == 1){
          VarCovRec <- as.matrix(mean(sigma))
          nome <- dimnames(object$receiver_model$sigma[[1]])[[1]]
          rownames(VarCovRec) <- colnames(VarCovRec) <- nome
        } else {
          VarCovRec <- rowMeans(sigma, dims = 2)
          VarCovRec <- (VarCovRec + t(VarCovRec))/2
        }
      } else {
        VarCovRec <- NULL
      }
    }
  }
  if(object$model == "tie"){
    return(VarCov)
  } else {
    return(list(sender_model = as.matrix(VarCovSnd),
                receiver_model = as.matrix(VarCovRec)))
  }

}

#####################################################################################
#####################################################################################
#####################################################################################

#In this part, we have the summary functions

#####################################################################################
#####################################################################################
#####################################################################################


#metarem <- function(remx, typeName = NULL){
#  UseMethod("metarem")
#}

#' summary.metarem
#' @title summary.metarem
#' @rdname summary.metarem
#' @description A function that arranges a summary of a 'metarem' object
#' @param object is a \code{metarem} object
#' @param ... further arguments to be passed to the 'summary' method
#' @method summary metarem
#' @export
summary.metarem <- function(object, ...){
  if (!inherits(object, "metarem"))
    warning("calling summary.metarem(<fake-metarem-object>) ...")
  summary_out <- list()
  if(object$model == "tie"){
    #Preparing a table with the fixed effects
    fe <- fixef_metarem(object)
    fixed_effects <- fe$fixed_effects
    std_fe <- fe$std_error
    samples <- abind::abind(object$mu, object$psi, along = 2)
    lower <- rowMeans(apply(samples, c(2, 3), stats::quantile, probs = c(0.025)))
    upper <- rowMeans(apply(samples, c(2, 3), stats::quantile, probs = c(0.975)))
    coefsTab <- cbind("post. mean" = fixed_effects,
                      "post. sd" = std_fe,
                      "2.5% CI" = lower,
                      "97.5% CI" = upper)
    summary_out$coefsTab <- coefsTab
    if(object$random){
      summary_out$VarCov <- VarCov_metarem(object)
    }
    attr(summary_out, "model") <- "tie"
    attr(summary_out, "Niter") <- object$Niter
    attr(summary_out, "NburnIn") <- object$NburnIn
    attr(summary_out, "Nchain") <- object$Nchain
    attr(summary_out, "Nthin") <- object$Nthin
    attr(summary_out, "run_time") <- object$run_time
  } else { #model = "actor"
    #Fixed-effects in the sender model
    fe <- fixef_metarem(object)
    fixed_effects <- fe$fixed_effects
    std_fe <- fe$std_error
    samplesSnd <- object$sender_model$mu
    lowerSnd <- rowMeans(apply(samplesSnd, c(2, 3), stats::quantile, probs = c(0.025)))
    upperSnd <- rowMeans(apply(samplesSnd, c(2, 3), stats::quantile, probs = c(0.975)))
    coefsTabSnd <- cbind("post. mean" = fixed_effects$sender_model,
                         "post. sd" = std_fe$sender_model,
                         "2.5% CI" = lowerSnd,
                         "97.5% CI" = upperSnd)
    #Fixed-effects in the receiver model
    samplesRec <- object$receiver_model$mu
    lowerRec <- rowMeans(apply(samplesRec, c(2, 3), stats::quantile, probs = c(0.025)))
    upperRec <- rowMeans(apply(samplesRec, c(2, 3), stats::quantile, probs = c(0.975)))
    coefsTabRec <- cbind("post. mean" = fixed_effects$receiver_model,
                         "post. sd" = std_fe$receiver_model,
                         "2.5% CI" = lowerRec,
                         "97.5% CI" = upperRec)
    summary_out$coefsTabSnd <- coefsTabSnd
    summary_out$coefsTabRec <- coefsTabRec
    if(object$random_sender | object$random_receiver){
      VarCov <- VarCov_metarem(object)
      summary_out$VarCovSnd <- VarCov$sender_model
      summary_out$VarCovRec <- VarCov$receiver_model
    }
    attr(summary_out, "model") <- "actor"
    attr(summary_out, "Niter") <- object$Niter
    attr(summary_out, "NburnIn") <- object$NburnIn
    attr(summary_out, "Nchain") <- object$Nchain
    attr(summary_out, "Nthin") <- object$Nthin
    attr(summary_out, "run_time") <- object$run_time
    attr(summary_out, "random_sender") <- object$random_sender
    attr(summary_out, "random_receiver") <- object$random_receiver
  }
  class(summary_out) <- "summary.metarem"
  return(summary_out)
}


#' print.summary.metarem
#' @title print.summary.metarem
#' @rdname print.summary.metarem
#' @description A function that prints out a summary of a 'metarem' object.
#' @param x is a \code{summary.metarem} object
#' @param ... further arguments to be passed to the 'print.summary' method
#' @method print summary.metarem
#' @export
print.summary.metarem <- function(x, ...){

  if(attr(x, "model") == "tie"){
    cat("Meta-analytic relation event model", paste0("(", attr(x, "model"), " oriented)"), "\n")
    cat("\n")
    cat("\n")
    cat("Fixed effects:\n\n")
    stats::printCoefmat(round(x$coefsTab,3))
    cat("\n")
    cat("Covariance matrix of random effects (posterior mean):\n\n")
    stats::printCoefmat(round(x$VarCov,3))
    cat("\n")
    cat("Statistics generated based on", attr(x, "Niter"), "MCMC samples.\n")
    cat("The first", attr(x, "NburnIn"), "were discared and one in every", attr(x, "Nthin"), "were kept.\n")
    cat(paste("The model took", format(round(attr(x, "run_time"),2), units = "mins"), "to run."))
  }
  if(attr(x, "model") == "actor"){ #x$model = "actor"
    cat("Meta-analytic relation event model", paste0("(", attr(x, "model"), " oriented)"), "\n")
    cat("\n")
    cat("\n")
    cat("Fixed effects (Sender Rate):\n\n")
    stats::printCoefmat(round(x$coefsTabSnd,3))
    cat("\n")
    cat("Fixed effects (Receiver Choice):\n\n")
    stats::printCoefmat(round(x$coefsTabRec,3))
    if(attr(x, "random_sender")){
      cat("\n")
      cat("Covariance matrix of random effects (Sender Rate) (posterior mean):\n\n")
      stats::printCoefmat(round(x$VarCovSnd,3))
    }
    if(attr(x, "random_receiver")){
      cat("\n")
      cat("Covariance matrix of random effects (Receiver Choice) (posterior mean):\n\n")
      stats::printCoefmat(round(x$VarCovRec,3))
    }
    cat("\n")
    cat("Statistics generated based on", attr(x, "Niter"), "MCMC samples.\n")
    cat("The first", attr(x, "NburnIn"), "were discared and one in every", attr(x, "Nthin"), "were kept.\n")
    cat(paste("The model took", format(round(attr(x, "run_time"),2), units = "mins"), "to run."))
  }

}

#' print.metarem
#' @title print.metarem
#' @rdname print.metarem
#' @description A function that prints out a summary of a 'metarem' object.
#' @param x is a \code{summary.metarem} object
#' @param ... further arguments to be passed to the 'print.summary' method
#' @method print metarem
#' @export
print.metarem <- function(x, ...){

  #If only print is called, then we call summary
  return(summary.metarem(x))

}

#####################################################################################
#####################################################################################
#####################################################################################

#In this part, we have the streamrem code

#####################################################################################
#####################################################################################
#####################################################################################

#'streamrem
#'
#'A function to estimate relational event models on data streams.
#'
#'@param data This is either a list containing edgelist and statistics for the batches of the networks (see example below), or a character string containing the paths where the data is saved in the hard drive. In case the data is in HD the batches should be saved separately (in format .rds) and each should be a list with edgelist and statistics.
#'@param actors A character vector containing the actors in the network, this should be particularly used when not every actor has been observed. See \code{?remify::remify}.
#'@param directed If \code{TRUE} then the network is treated as directed.
#'@param ordinal Set to \code{TRUE} in case only the order of the events is known.
#'@param update Set to \code{TRUE} to update the estimates with new batches. You should provide a model previously fitted with the \code{streamrem} function.
#'@param model A model previously fitted with the \code{streamrem} function.
#'
#'@examples
#'
#'#First install remulate, remstats and remify
#'#library(devtools)
#'#install_github("TilburgNetworkGroup/remulate")
#'#devtools::install_github("TilburgNetworkGroup/remify")
#'#devtools::install_github("TilburgNetworkGroup/remstats")
#'
#' edgelist <- stream
#'
#' library(remstats)
#'
#' names(edgelist) <- c("time", "actor1", "actor2")
#'
#' actors <- as.character(1:10)
#'
#' #We will devide the network in 10 batches of 500
#' events <- seq(1, 5001, by = 500)
#'
#' #Declaring which effects we want remstats to compute
#' effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#'
#' #Getting the remify object
#' rehObj <- remify::remify(edgelist, model = "tie", actors = actors)
#'
#' data <- vector("list")
#'
#' for(i in 2:length(events)){
#'   #Computing statistics
#'   stats <- remstats::tomstats(effects, rehObj, start = events[i-1], stop = events[i]-1)
#'   edl <- edgelist[events[i-1]:(events[i]-1),]
#'   #Every piece needs to be stored a in a list with edgelist and statistics
#'   data[[i-1]] <- list(edgelist = edl,
#'                       statistics = stats)
#' }
#
#' #Let's compute the effects for the first 7 batches of the networks
#' fit <- streamrem(data[1:7], actors)
#'
#' #printing the parameters
#' print(fit)
#' #Plotting a trend line of the estimates
#' plot(fit)
#'
#' #Now we can update the model with the remaining 3 batches
#' fit2 <- streamrem(data[8:10], actors, update = TRUE, model = fit)
#'
#' #printing the parameters
#' print(fit2)
#' #Plot a trend line with the estimates
#' plot(fit2)
#'
#'@export
streamrem <- function(data,
                      actors = NULL,
                      directed = TRUE,
                      ordinal = FALSE,
                      update = F,
                      model = NULL){
  #Storing the results
  estimates <- vector("list")
  #Checking whether we need to update the model
  if(update & is.null(model)){
    stop("update = TRUE requires a model object previously fitted with streamrem.")
  }
  #Number of batches
  L <- length(data)

  #Looping through the partitions
  tt <- Sys.time()
  for(i in 1:L){
    #data has to be either a character vector containing the path where the data is saved
    #or a list containing the batches, every batch has to be a list with edgelist and statistics
    if(is.character(data)){
      dat <- readRDS(data[i])
    } else {
      dat <- data[[i]]
    }
    #Checking if the list is names
    if(is.null(names(dat))){
      stop("The list must be named, containing edgelist and statistics.")
    }
    #Extracting the edgelist, the list has to be named
    edgelist <- dat$edgelist
    if(i == 1){
      #store the time of the last event to make the subtraction
      tempo <- edgelist[nrow(edgelist),1]
      #Creating the matrices to store the estimates
      p <- dim(dat$statistics)[3]
      #Storing streaming estimates
      beta <- matrix(NA, nrow = p, ncol = L)
      omega <- array(NA, dim = c(p, p, L))
      if(update){
        edgelist[,1] <- edgelist[,1] - model$timeL
        edgelist[,1] <- as.numeric(edgelist[,1])
      }
    } else {
      #store the time of the last event to make the subtraction
      tempo <- edgelist[nrow(edgelist),1]
      edgelist[,1] <- edgelist[,1] - timeF
      edgelist[,1] <- as.numeric(edgelist[,1])
    }
    #Extracting statistics, the list have to be named
    stats <- dat$statistics
    class(stats) <- c("tomstats", "remstats")
    edgelistREH <- suppressWarnings(
      remify::remify(edgelist, ordinal = F, actors = actors, model = "tie", directed = directed)
    )

    #class(dataList[[2]]) <- c("tomstats", "remstats")
    print(paste("Obtaining the MLE estimates."))
    rem <- suppressWarnings(
      remstimate::remstimate(reh = edgelistREH,
                             stats = stats,
                             method = "MLE")
    )
    print(paste("The MLE estimates were obtained."))
    #Changing the name of tempo, so we don't f*ck it up in the next iteration
    timeF <- tempo

    estimates[[i]] <- list(coef = rem$coefficients,
                           var = (rem$vcov + t(rem$vcov)/2))

    #Getting the regression coefficients and covariance matrices
    if(i == 1){
      if(update){
        l <- dim(model$omega)[3] #last batch of the previous fit
        omega[,,i] <- solve(solve(estimates[[i]]$var) + solve(model$omega[,,l]))
        beta[,i] <- omega[,,i] %*% (solve(estimates[[i]]$var) %*% estimates[[i]]$coef + solve(model$omega[,,l]) %*% model$beta[,l])
      } else {
        omega[,,i] <- estimates[[i]]$var
        beta[,i] <- estimates[[i]]$coef
      }
    } else{
      omega[,,i] <- solve(solve(estimates[[i]]$var) + solve(omega[,,i-1]))
      beta[,i] <- omega[,,i] %*% (solve(estimates[[i]]$var) %*% estimates[[i]]$coef + solve(omega[,,i-1]) %*% beta[,i-1])
    }
    print(paste("This is iteration", i, "out of", length(data)))
  }
  (tt2 <- Sys.time() - tt)
  #Giving names to the estimates
  statsnames <- dimnames(stats)[[3]]
  if(update){
    beta <- cbind(model$beta, beta)
    omega <- abind::abind(model$omega, omega, along = 3)
    colnames(beta) <- dimnames(omega)[[3]] <- paste0("batch_", 1:ncol(beta))
  } else {
    colnames(beta) <- dimnames(omega)[[3]] <- paste0("batch_", 1:L)
  }
  rownames(beta) <- statsnames
  dimnames(omega)[[1]] <- dimnames(omega)[[2]] <- statsnames
  #Output
  final_list <- list(beta = beta,
                     omega = omega,
                     run_time = tt2,
                     timeL = timeF)
  #Computing intervals
  conf_intervals <- intervals(final_list)
  final_list$confint <- conf_intervals
  class(final_list) <- "streamrem"
  return(final_list)
}


#################################################################################
#################################################################################
#################################################################################

#'interval
#'
#'@param fit A streamrem object.
#'@param prob A probability to compute the confidence interval, default is 0.95.
#'@export
intervals <- function(fit,
                      prob = 0.95){
  #Extracting the standard errors of from the covariance matrices
  stddev <- sapply(1:dim(fit$omega)[3], function(x) diag((fit$omega[,,x])^.5))
  #Computing upper bounds
  up <- fit$beta - stats::qnorm((1-prob)/2) * stddev
  #Computing lower bounds
  low <- fit$beta + stats::qnorm((1-prob)/2) * stddev
  #number of covariates in the model
  n_covs <- nrow(low)
  covs <- rownames(low)
  #creating the output list
  output <- lapply(1:n_covs, function(x) data.frame(lower = low[x,], upper = up[x,]))
  names(output) <- covs
  return(output)
}

#################################################################################
#################################################################################
#################################################################################

#'plot.streamrem
#'
#'Plots lines with the trends of each parameter estimates and 95% confidence intervals.
#'@param x A streamrem object.
#'@param same_page The default is FALSE. If TRUE, the function will create all graphs in one plot.
#'@param ... further arguments to be passed to the 'summary' method
#'@export
plot.streamrem <- function(x, same_page = FALSE, ...){
  if (!inherits(x, "streamrem"))
    warning("calling streamrem(<fake-streamrem-object>) ...")
  beta <- x$beta
  intervals <- x$confint
  covs <- rownames(beta)
  if(same_page){
    graphics::par(mfrow = c(nrow(beta)/2, 2))
  }
  for(i in 1:nrow(beta)){
    xlim <- c(1, ncol(beta))
    scl <- c(x$confint[[i]][,1],x$confint[[i]][,2])
    ylim <- range(scl)
    graphics::plot.new()
    graphics::plot.window(xlim, ylim, main = covs[i])
    graphics::lines(beta[i,])
    graphics::lines(x$confint[[i]][,1], lty = 2)
    graphics::lines(x$confint[[i]][,2], lty = 2)
    graphics::axis(side = 1, at = seq(xlim[1], xlim[2]))
    graphics::axis(side = 2, at = round(c(min(ylim), beta[i,1], max(ylim)),4))
    graphics::box(which = "plot")
    graphics::title(main = covs[i], ylab = "values", xlab = "batches")
  }
}

#################################################################################
#################################################################################
#################################################################################

#' summary.streamrem
#' @title summary.streamrem
#' @rdname summary.streamrem
#' @description A function that arranges a summary of a 'streamrem' object
#' @param object is a \code{streamrem} object
#' @param ... further arguments to be passed to the 'summary' method
#' @method summary streamrem
#' @export
summary.streamrem <- function(object, ...){
  if (!inherits(object, "streamrem"))
    warning("calling streamrem(<fake-streamrem-object>) ...")
  summary_out <- vector("list")
  summary_out$beta <- object$beta
  summary_out$run_time <- object$run_time
  class(summary_out) <- "summary.streamrem"
  return(summary_out)
}

#################################################################################
#################################################################################
#################################################################################

#' print.summary.streamrem
#' @title print.summary.streamrem
#' @rdname print.summary.streamrem
#' @description A function that prints out a summary of a 'metarem' object.
#' @param x is a \code{summary.streamrem} object
#' @param ... further arguments to be passed to the 'print.summary' method
#' @method print summary.streamrem
#' @export
print.summary.streamrem <- function(x, ...){
  cat("Meta-analytic relational event model. \n")
  cat("\n")
  cat(paste("The model ran", nrow(x$beta), "covariates in total."))
  cat("\n")
  stats::printCoefmat(t(x$beta),3)
  cat("\n")
  cat(paste("The model took",  format(round(x$run_time,2), units = "secs"), "to run."))
}

#################################################################################
#################################################################################
#################################################################################

#' print.streamrem
#' @title print.streamrem
#' @rdname print.streamrem
#' @description A function that prints out a summary of a 'streamrem' object.
#' @param x is a \code{summary.streamrem} object
#' @param ... further arguments to be passed to the 'print.summary' method
#' @method print streamrem
#' @export
print.streamrem <- function(x, ...){

  #If only print is called, then we call summary
  return(summary.streamrem(x))

}
