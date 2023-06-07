#' remx
#'
#' A function to fit mixed-effect relational event models. It runs meta-analytic approximations for actor- and tie-oriented relational event models. This function also supports full multilevel models from \code{lme4} and \code{rstanarm}.
#'
#' @param edgelist a list containing multiple relational event sequences
#' @param statistics a list containing multiple \code{remstats} objects
#' @param random a vector containing the names of the covariates that will be treated as random effects, use it for \code{model = "tie"}.
#' @param fixed a vector containing the names of the covariates that will be treated as fixed effects, use it for \code{model = "tie"}.
#' @param random_sender a vector containing the names of the covariates that will be treated as random effects in the sender model, use it for \code{model = "actor"}.
#' @param fixed_sender a vector containing the names of the covariates that will be treated as fixed effects in sender model, use it for \code{model = "actor"}.
#' @param random_receiver a vector containing the names of the covariates that will be treated as random effects in the receiver model, use it for \code{model = "actor"}.
#' @param fixed_receiver a vector containing the names of the covariates that will be treated as fixed effects in the receiver model, use it for \code{model = "actor"}.
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
#' @return \code{beta} a list containing arrays of random effects. Each element is a array with dimensions corresponding to samples, effects, and networks, respectively
#' @return \code{mu} random-effect means, an array with dimensions corresponding to samples, effects and chains
#' @return \code{sigma} random-effect covariance matrix
#' @return \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
#' @return \code{psi} fixed effects, an array with dimensions samples, effects and chains
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
#' @return Inside the sender_rate list:
#' @return \code{gamma} a list containing random effects in the sender model. Each array has dimensions samples, effects and networks
#' @return \code{mu} random-effect means, the array has dimensions samples, effects and chains
#' @return \code{sigma} random-effect covariance matrix
#' @return \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
#' @return \code{phi} fixed effects, array dimensions are samples, effects, and chains
#' @return \code{omega} scale-mixture parameter that allows beta to have a student T likelihood
#' @return In the receiver_choice list:
#' @return \code{beta} a list containing arrays of random effects. Each element is a array with dimensions corresponding to samples, effects, and networks, respectively
#' @return \code{mu} random-effect means, an array with dimensions corresponding to samples, effects and chains
#' @return \code{sigma} random-effect covariance matrix
#' @return \code{alpha} mixture parameter that allows \code{sigma} to have a half-T prior
#' @return \code{psi} fixed effects, an array with dimensions samples, effects and chains
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
#' edgelist <- lapply(edgelist, function(x) x[,])
#' for(i in 1:length(edgelist)){names(edgelist[[i]]) <- c("time", "actor1", "actor2")}
#'
#' #Computing statistics
#' effects <- ~ remstats::inertia(scaling = "std") + remstats::reciprocity(scaling = "std")
#' stats <- lapply(1:length(edgelist), function(x) remstats::tomstats(effects, edgelist[[x]]))
#'
#' #Running the model
#'
#' fit <- remx(edgelist = edgelist,
#'             statistics = stats,
#'             random = c("baseline"),
#'             fixed = c("inertia", "reciprocity"),
#'             method = "meta",
#'             model = "tie")
#'
#' print(fit)
#'
#' #----------------------------#
#' #    Actor-Oriented model    #
#' #----------------------------#
#'
#' #Computing statistics
#' sender_effects <- ~ indegreeSender(scaling = "std") + outdegreeSender(scaling = "std")
#' receiver_effects <- ~ indegreeReceiver(scaling = "std") + rrankSend()
#' stats <- lapply(1:length(edgelist), function(x) remstats::aomstats(edgelist[[x]], sender_effects, receiver_effects))
#'
#' #Running the model
#' fit <- remx(edgelist = edgelist,
#'             statistics = stats,
#'             random_sender = c("baseline"),
#'             fixed_sender = c("indegreeSender"),
#'             random_receiver = c("indegreeReceiver", "rrankSend"),
#'             fixed_receiver = NULL,
#'             method = "meta",
#'             model = "actor")
#'
#' print(fit)
#' @export
remx <- function(edgelist,
                 statistics,
                 random = NULL,
                 fixed = NULL,
                 random_sender = NULL,
                 fixed_sender = NULL,
                 random_receiver = NULL,
                 fixed_receiver = NULL,
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
  if(model == "tie"){
    if(is.null(random) & is.null(fixed)){
      stop("Both random and fixed parameters are NULL. At least one of them must contain the name of one or more covariables.")
    }
  }
  if(model == "actor"){
    if((is.null(random_sender) & is.null(fixed_sender)) | is.null(random_receiver) & is.null(fixed_receiver)){
      stop("Either random_sender and fixed_sender or random_receiver and fixed_receiver are both NULL. At least one of each must contain the name of one or more covariables.")
    }
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
                  fixed = fixed,
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
                       fixed_sender = fixed_sender,
                       fixed_receiver = fixed_receiver,
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
#' @param fixed a vector containing the names of the covariates that will be treated as fixed effects.
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
                 fixed = NULL,
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

  if(is.null(random) & is.null(fixed)){
    return(
      warning("Both random and fixed parameters are NULL. At least one of them must contain the name of one or more covariables.")
    )
  }

  estimates <- vector("list")
  t1 <- Sys.time()
  for(i in 1:length(edgelist)){

    edgelistREH <- remify::reh(edgelist[[i]], model = "tie", actors = actors[[i]],
                               directed = directed, ordinal = ordinal)

    rem <- remstimate::remstimate(reh = edgelistREH,
                                  stats = statistics[[i]],
                                  method = "MLE",
                                  model = "tie",
    )

    estimates[[i]] <- list(coef = rem$coefficients, var = rem$vcov)

    print(paste("MLE for cluster", i, "has been obtained."))

  }

  names(estimates) <- paste0("cluster_", 1:length(edgelist))


  #checking whether all statistics array contain the same number of covariates
  K <- length(edgelist) #number of groups

  if(!is.null(random)){
    p <- length(random) #number random effects
    #We assume that all statistics array will contain the same statistics
    beta_hat <- as.matrix(sapply(estimates, function(x) (x$coef[random])))
    omega_hat <- lapply(estimates, function(x) as.matrix(x$var[random,random]))
    if(p == 1){
      omega_hat <- array(simplify2array(omega_hat), dim = c(p,p,K))
      beta_hat <- t(beta_hat)
    } else {
      omega_hat <- simplify2array(omega_hat)
    }
    fixedModel <- FALSE
  } else {
    p <- 0
    beta_hat <- matrix(0, nrow = 1, ncol = 1)
    omega_hat <- array(0, dim = c(1,1,K))
    fixedModel <- TRUE
  }


  if(!is.null(fixed)){
    q <- length(fixed) #number of fixed effects
    psi <- as.matrix(sapply(estimates, function(x) x$coef[fixed]))
    SigmaPsi <- lapply(estimates, function(x) x$var[fixed,fixed])
    if(q > 1){
      SigmaPsi <- sapply(SigmaPsi, diag)
    } else if(q == 1){
      SigmaPsi <- matrix(simplify2array(SigmaPsi), ncol = K)
      psi <- t(psi)
    }
    randomModel <- FALSE
  } else {
    psi <- matrix(0, nrow = 1, ncol = 1)
    SigmaPsi <- matrix(0, nrow = 1, ncol = 1)
    q <- 0
    randomModel <- TRUE
  }

  nu <- sapply(edgelist, nrow) - p #vectors of df's
  if(sum(nu < 0) != 0){
    return(paste("Nu parameter must be greater than zero for all groups! You either need more events or less covariates."))
  }

  if(p == 0){
    mu <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,1))) #random-effect means
    #I'm making a confusing, using lambda here, but on the paper, these are called alphas
    Lambda <- matrix(stats::rgamma(1*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    omega <- matrix(stats::rgamma(K*Nchain,1,1),ncol = Nchain) #scale mixture for likelihood
    #initial values for MCMC
    beta <- array(NA, dim = c(1, K, Nchain)) #group-specific effects
    sigma <- array(NA, dim = c(1, 1, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      beta[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,1)))
      sigma[,,i] <- diag(stats::rpois(1,100), nrow = 1)
    }
  } else {
    mu <- t(mvtnorm::rmvnorm(Nchain, mean = rep(0,p))) #random-effect means
    #I'm making a confusing, using lambda here, but on the paper, these are called alphas
    Lambda <- matrix(stats::rgamma(p*Nchain, 1, 1), ncol = Nchain) #scale mixture for cov mat
    omega <- matrix(stats::rgamma(K*Nchain,1,1),ncol = Nchain) #scale mixture for likelihood
    #initial values for MCMC
    beta <- array(NA, dim = c(p, K, Nchain)) #group-specific effects
    sigma <- array(NA, dim = c(p, p, Nchain)) #random-effect covariance matrix
    for(i in 1:Nchain){
      beta[,,i] <- t(mvtnorm::rmvnorm(K, mean = rep(0,p)))
      sigma[,,i] <- diag(stats::rpois(p,100), nrow = p)
    }
  }

  #Prior for Mu and Sigma
  eta <- 2
  xi <- 10 #it can be an arbitrarily large number (Gelman, 2006)
  TauMu <- diag(10, p) #scale matrix of mu
  TauPsi <- rep(10, q) #scale matrix of psi

  print(paste("Running the Gibbs sampler."))

  MCMC <- sampler(Niter = Niter,
                  Nchain = Nchain,
                  NburnIn = NburnIn,
                  Nthin = Nthin,
                  p = p,
                  q = q,
                  K = K,
                  nu = nu,
                  beta = beta,
                  beta_hat = beta_hat,
                  omega_hat = omega_hat,
                  Mu = mu,
                  Sigma = sigma,
                  eta = eta,
                  Lambda = Lambda,
                  xi = xi,
                  TauMu = TauMu,
                  TauPsi = TauPsi,
                  SigmaPsi = SigmaPsi,
                  omega = omega,
                  psi = psi,
                  randomModel = randomModel,
                  fixedModel =  fixedModel)

  if(!is.null(fixed)){
    colnames(MCMC$psi) <- fixed
  }
  if(!is.null(random)){
    colnames(MCMC$mu) <- random
    colnames(MCMC$alpha) <- random
    MCMC$MLE <- estimates
    for(i in 1:Nchain){
      colnames(MCMC$beta[[i]]) <- random
      colnames(MCMC$sigma[[i]]) <- rownames(MCMC$sigma[[i]]) <- random
      dimnames(MCMC$beta[[i]])[[3]] <- paste0("cluster_", 1:K)
    }
    names(MCMC$beta) <- paste0("chain_", 1:Nchain)
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
#' @param fixed_sender a vector containing the names of the covariates that will be treated as fixed effects in sender model.
#' @param random_receiver a vector containing the names of the covariates that will be treated as random effects in the receiver model.
#' @param fixed_receiver a vector containing the names of the covariates that will be treated as fixed effects in the receiver model.
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
                      random_sender = NULL,
                      fixed_sender = NULL,
                      random_receiver = NULL,
                      fixed_receiver = NULL,
                      actors = NULL,
                      ordinal = FALSE,
                      directed = TRUE,
                      Niter = 10000,
                      Nchain = 2,
                      NburnIn = 5000,
                      Nthin = 10){
  #HERE THE FUNCTION STARTS

  if(is.null(random_sender) & is.null(fixed_sender)){
    return(
      warning("Both random_sender and fixed_sender parameters are NULL. At least one of them must contain the name of one or more covariables.")
    )
  }
  if(is.null(random_receiver) & is.null(fixed_receiver)){
    return(
      warning("Both random_receiver and fixed_receiver parameters are NULL. At least one of them must contain the name of one or more covariables.")
    )
  }

  estimates <- vector("list")
  t1 <- Sys.time()
  for(i in 1:length(edgelist)){
    edgelistREH <- remify::reh(edgelist[[i]], model = "actor", actors = actors[[i]],
                               directed = directed, ordinal = ordinal)

    rem <- remstimate::remstimate(reh = edgelistREH,
                                  stats = statistics[[i]],
                                  method = "MLE",
                                  model = "actor",
    )

    estimates[[i]] <- list(sender = list(coef = rem$sender_rate$coefficients, var = rem$sender_rate$vcov),
                           receiver = list(coef = rem$receiver_choice$coefficients, var = rem$receiver_choice$vcov))

    print(paste("MLE for cluster", i, "has been obtained."))

  }

  names(estimates) <- paste0("cluster_", 1:length(edgelist))

  #checking whether all statistics array contain the same number of covariates
  K <- length(edgelist) #number of groups

  #Checking for random effects
  #SENDER MODEL
  if(!is.null(random_sender)){
    p <- length(random_sender) #number random effects sender model
    #Sender model
    gamma_hat <- as.matrix(sapply(estimates, function(x) (x$sender$coef[random_sender])))
    zeta_hat <- lapply(estimates, function(x) as.matrix(x$sender$var[random_sender,random_sender]))
    #Fixing the structure of the random effects in the sender model
    if(p == 1){
      zeta_hat <- array(simplify2array(zeta_hat), dim = c(p,p,K))
      gamma_hat <- t(gamma_hat)
    } else {
      zeta_hat <- simplify2array(zeta_hat)
    }
    fixModelSnd <- FALSE
  } else {
    p <- 0
    gamma_hat <- matrix(0, nrow = 1, ncol = 1)
    zeta_hat <- array(0, dim = c(1,1,K))
    fixModelSnd <- TRUE
  }

  #RECEIVER MODEL
  if(!is.null(random_receiver)){
    v <- length(random_receiver) #random effects receiver model
    #Receiver model
    beta_hat <- as.matrix(sapply(estimates, function(x) (x$receiver$coef[random_receiver])))
    omega_hat <- lapply(estimates, function(x) as.matrix(x$receiver$var[random_receiver,random_receiver]))
    #Fixing the structure of the random effects in the receiver model
    if(v == 1){
      omega_hat <- array(simplify2array(omega_hat), dim = c(v,v,K))
      beta_hat <- t(beta_hat)
    } else {
      omega_hat <- simplify2array(omega_hat)
    }
    fixModelRec <- FALSE
  } else {
    v <- 0
    beta_hat <- matrix(0, nrow = 1, ncol = 1)
    omega_hat <- array(0, dim = c(1,1,K))
    fixModelRec <- TRUE
  }


  #Checking if the receiver model contains fixed effects
  if(!is.null(fixed_receiver)){
    u <- length(fixed_receiver) #number of fixed effects
    psi <- as.matrix(sapply(estimates, function(x) x$receiver$coef[fixed_receiver]))
    SigmaPsi <- lapply(estimates, function(x) x$receiver$var[fixed_receiver,fixed_receiver])
    if(u > 1){
      SigmaPsi <- sapply(SigmaPsi, diag)
    } else if(u == 1){
      SigmaPsi <- matrix(simplify2array(SigmaPsi), ncol = K)
      psi <- t(psi)
    }
    randomModelRec <- FALSE
  } else {
    psi <- matrix(0, nrow = 1, ncol = 1)
    SigmaPsi <- matrix(0, nrow = 1, ncol = 1)
    u <- 0
    randomModelRec <- TRUE
  }

  #Checking if the sender model contains fixed effects
  if(!is.null(fixed_sender)){
    q <- length(fixed_sender) #number of fixed effects
    phi <- as.matrix(sapply(estimates, function(x) x$sender$coef[fixed_sender]))
    SigmaPhi <- lapply(estimates, function(x) x$sender$var[fixed_sender,fixed_sender])
    if(q > 1){
      SigmaPhi <- sapply(SigmaPhi, diag)
    } else if(q == 1){
      SigmaPhi <- matrix(simplify2array(SigmaPhi), ncol = K)
      phi <- t(phi)
    }
    randomModelSnd <- FALSE
  } else {
    phi <- matrix(0, nrow = 1, ncol = 1)
    SigmaPhi <- matrix(0, nrow = 1, ncol = 1)
    q <- 0
    randomModelSnd <- TRUE
  }

  ###########################################################################################
  ###########################################################################################
  ###########################################################################################

  #Setting up the additional quantities needed

  nu <- sapply(edgelist, nrow) - p #vectors of df's
  if(sum(nu < 0) != 0){
    return(paste("Nu parameter must be greater than zero for all groups! You either need more events or less covariates."))
  }

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

  #Initial values for the scale mixture parameter
  omega_gamma <- matrix(stats::rgamma(K*Nchain,1,1),ncol = Nchain) #scale mixture for likelihood
  omega_beta <- matrix(stats::rgamma(K*Nchain,1,1),ncol = Nchain) #scale mixture for likelihood

  #Prior for Mu and Sigma
  eta <- 2
  xi <- 10 #it can be an arbitrarily large number (Gelman, 2006)
  TauMu <- diag(10, v) #scale matrix of mu
  TauPsi <- rep(10, u) #scale matrix of psi
  TauMuG <- diag(10, p) #scale matrix of mu
  TauPhi <- rep(10, q) #scale matrix of psi

  print(paste("Running the Gibbs sampler."))

  MCMC <- samplerActor(Niter = Niter,
                       Nchain = Nchain,
                       NburnIn = NburnIn,
                       Nthin = Nthin,
                       v = v,
                       u = u,
                       K = K,
                       nu = nu,
                       beta = beta,
                       beta_hat = beta_hat,
                       omega_hat = omega_hat,
                       MuB = mu_beta,
                       Sigma = sigma,
                       eta = eta,
                       LambdaB =  Lambda_beta,
                       xi = xi,
                       TauMuB = TauMu,
                       TauPsi = TauPsi,
                       SigmaPsi = SigmaPsi,
                       omegaB = omega_beta,
                       psi = psi,
                       p = p,
                       q = q,
                       gamma = gamma,
                       gamma_hat = gamma_hat,
                       zeta_hat = zeta_hat,
                       MuG = mu_gamma,
                       Zeta = zeta,
                       LambdaG = Lambda_gamma,
                       TauMuG = TauMuG,
                       TauPhi = TauPhi,
                       SigmaPhi = SigmaPhi,
                       omegaG = omega_gamma,
                       phi = phi,
                       randomModelSnd = randomModelSnd,
                       randomModelRec = randomModelRec,
                       fixModelSnd = fixModelSnd,
                       fixModelRec = fixModelRec)

  print(paste("The Gibbs sampler has finished."))

  #Receiver model
  if(!is.null(fixed_receiver)){colnames(MCMC$receiver_choice$psi) <- fixed_receiver}
  if(!is.null(random_receiver)){
    colnames(MCMC$receiver_choice$mu) <- random_receiver
    colnames(MCMC$receiver_choice$alpha) <- random_receiver
    MCMC$MLE <- estimates
    for(i in 1:Nchain){
      colnames(MCMC$receiver_choice$beta[[i]]) <- random_receiver
      colnames(MCMC$receiver_choice$sigma[[i]]) <- rownames(MCMC$receiver_choice$sigma[[i]]) <- random_receiver
      dimnames(MCMC$receiver_choice$beta[[i]])[[3]] <- paste0("cluster_", 1:K)
    }
    names(MCMC$receiver_choice$beta) <- paste0("chain_", 1:Nchain)
    names(MCMC$receiver_choice$sigma) <- paste0("chain_", 1:Nchain)
  }

  #Sender model
  if(!is.null(fixed_sender)){colnames(MCMC$sender_rate$phi) <- fixed_sender}
  if(!is.null(random_sender)){
    colnames(MCMC$sender_rate$mu) <- random_sender
    colnames(MCMC$sender_rate$alpha) <- random_sender
    MCMC$MLE <- estimates
    for(i in 1:Nchain){
      colnames(MCMC$sender_rate$gamma[[i]]) <- random_sender
      colnames(MCMC$sender_rate$sigma[[i]]) <- rownames(MCMC$sender_rate$sigma[[i]]) <- random_sender
      dimnames(MCMC$sender_rate$gamma[[i]])[[3]] <- paste0("cluster_", 1:K)
    }
    names(MCMC$sender_rate$gamma) <- paste0("chain_", 1:Nchain)
    names(MCMC$sender_rate$sigma) <- paste0("chain_", 1:Nchain)
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
      a <- abind::abind(object$mu, object$psi, along = 2)
      fe <- rowMeans(colMeans(a))
      std <- rowMeans(apply(a, c(2,3), stats::sd))
    } else { #model == "actor
      a <- abind::abind(object$sender_rate$mu, object$sender_rate$phi, along = 2)
      snd <- rowMeans(colMeans(a))
      stdA <- rowMeans(apply(a, c(2,3), stats::sd))
      b <- abind::abind(object$receiver_choice$mu, object$receiver_choice$psi, along = 2)
      rec <- rowMeans(colMeans(b))
      stdB <- rowMeans(apply(b, c(2,3), stats::sd))
      fe <- list(sender_rate = snd,
                 receiver_choice = rec)
      std <- list(sender_rate = stdA,
                  receiver_choice = stdB)
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
        variable <- dimnames(object$beta[[1]])[[2]]
        number_rnd_eff <- dim(object$beta[[1]])[2]
        if(number_rnd_eff == 1){
          rnd <- matrix(rowMeans(sapply(object$beta, function(y) apply(y, 3, colMeans))), ncol = 1)
          std <- matrix(rowMeans(sapply(object$beta, function(y) apply(y, 3, stats::sd))), ncol = 1)
          colnames(rnd) <- colnames(std) <- variable
        } else {
          if(is.na(object$random)){
            stop("Model does not contain random effects.")
          } else{
            a <- lapply(object$beta, function(y) t(colMeans(y)))
            a <- simplify2array(a)
            rnd <- t(apply(a, 1, rowMeans))
            std <- lapply(object$beta, function(y) t(apply(y, c(2,3), stats::sd)))
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
          variable <- dimnames(object$sender_rate$gamma[[1]])[[2]]
          number_rnd_snd <- dim(object$sender_rate$gamma[[1]])[2]
          if(number_rnd_snd == 1){
            rnd_snd <- matrix(rowMeans(sapply(object$sender_rate$gamma, function(y) apply(y, 3, colMeans))), ncol = 1)
            std_snd <- matrix(rowMeans(sapply(object$sender_rate$gamma, function(y) apply(y, 3, stats::sd))), ncol = 1)
            colnames(rnd_snd) <- colnames(std_snd) <- variable
          } else {
            #here comes the case of multiple random effects
            #here comes the case of multiple random effects
            a <- lapply(object$sender_rate$gamma, function(y) t(colMeans(y)))
            a <- simplify2array(a)
            rnd_snd <- t(apply(a, 1, rowMeans))
            std <- lapply(object$sender_rate$gamma, function(y) t(apply(y, c(2,3), stats::sd)))
            std <- simplify2array(std)
            std_snd <- t(apply(std, 1, rowMeans))
          }
        }
        if(!object$random_receiver){
          rnd_rec <- NULL
          std_rec <- NULL
        } else {
          variable <- dimnames(object$receiver_choice$beta[[1]])[[2]]
          number_rnd_rec <- dim(object$receiver_choice$beta[[1]])[2]
          if(number_rnd_rec == 1){
            rnd_rec <- matrix(rowMeans(sapply(object$receiver_choice$beta, function(y) apply(y, 3, colMeans))), ncol = 1)
            std_rec <- matrix(rowMeans(sapply(object$receiver_choice$beta, function(y) apply(y, 3, stats::sd))), ncol = 1)
            colnames(rnd_rec) <- colnames(std_rec) <- variable
          } else {
            #here comes the case of multiple random effects
            a <- lapply(object$receiver_choice$beta, function(y) t(colMeans(y)))
            a <- simplify2array(a)
            rnd_rec <- t(apply(a, 1, rowMeans))
            std <- lapply(object$receiver_choice$beta, function(y) t(apply(y, c(2,3), stats::sd)))
            std <- simplify2array(std)
            std_rec<- t(apply(std, 1, rowMeans))
          }
        }
      }
      rnd <- list(sender_rate = rnd_snd,
                  receiver_choice = rnd_rec)
      std <- list(sender_rate = std_snd,
                  receiver_choice = std_rec)
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
        sigma <- lapply(object$sender_rate$sigma, function(x) rowMeans(x, dims = 2))
        sigma <- simplify2array(sigma)
        if(dim(object$sender_rate$sigma[[1]])[1] == 1){
          VarCovSnd <- as.matrix(mean(sigma))
          nome <- dimnames(object$sender_rate$sigma[[1]])[[1]]
          rownames(VarCovSnd) <- colnames(VarCovSnd) <- nome
        } else {
          VarCovSnd <- rowMeans(sigma, dims = 2)
          VarCovSnd <- (VarCovSnd + t(VarCovSnd))/2
        }
      } else {
        VarCovSnd <- NULL
      }
      if(object$random_receiver){
        sigma <- lapply(object$receiver_choice$sigma, function(x) rowMeans(x, dims = 2))
        sigma <- simplify2array(sigma)
        if(dim(object$receiver_choice$sigma[[1]])[1] == 1){
          VarCovRec <- as.matrix(mean(sigma))
          nome <- dimnames(object$receiver_choice$sigma[[1]])[[1]]
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
    return(list(sender_rate = as.matrix(VarCovSnd),
                receiver_choice = as.matrix(VarCovRec)))
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
    samplesSnd <- abind::abind(object$sender_rate$mu, object$sender_rate$phi, along = 2)
    lowerSnd <- rowMeans(apply(samplesSnd, c(2, 3), stats::quantile, probs = c(0.025)))
    upperSnd <- rowMeans(apply(samplesSnd, c(2, 3), stats::quantile, probs = c(0.975)))
    coefsTabSnd <- cbind("post. mean" = fixed_effects$sender_rate,
                         "post. sd" = std_fe$sender_rate,
                         "2.5% CI" = lowerSnd,
                         "97.5% CI" = upperSnd)
    #Fixed-effects in the receiver model
    samplesRec <- abind::abind(object$receiver_choice$mu, object$receiver_choice$psi, along = 2)
    lowerRec <- rowMeans(apply(samplesRec, c(2, 3), stats::quantile, probs = c(0.025)))
    upperRec <- rowMeans(apply(samplesRec, c(2, 3), stats::quantile, probs = c(0.975)))
    coefsTabRec <- cbind("post. mean" = fixed_effects$receiver_choice,
                         "post. sd" = std_fe$receiver_choice,
                         "2.5% CI" = lowerRec,
                         "97.5% CI" = upperRec)
    summary_out$coefsTabSnd <- coefsTabSnd
    summary_out$coefsTabRec <- coefsTabRec
    if(object$random_sender | object$random_receiver){
      VarCov <- VarCov_metarem(object)
      summary_out$VarCovSnd <- VarCov$sender_rate
      summary_out$VarCovRec <- VarCov$receiver_choice
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
