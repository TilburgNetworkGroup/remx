//#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <RcppDist.h>
// [[Rcpp::depends(RcppArmadillo, RcppDist)]]
using namespace Rcpp;
using namespace arma;

//' mat2DF
//' 
//' A function that converts a matrix to data frame
//' 
//' @param x a matrix to be converted into a data frame
//' 
//' @return x converted to a data frame
//'
// [[Rcpp::export]]
Rcpp::DataFrame mat2DF(const arma::mat& x) {
  Function asDF("as.data.frame");
  return asDF(x);
}

//' rehPoisson
//' 
//' A function that reshapes the data set to allow implementation of the REM writen as Poisson regression
//' 
//' @param stats a remstats object, the objects inside the list must be in the same order as they would be in remstats 
//' 
//' @return a data frame with events stacked 
//'
// [[Rcpp::export]]
Rcpp::DataFrame rehPoisson(Rcpp::List stats){
  
  //stats is a list from remstats
  
  //extracting statistics array
  arma::cube statistics = stats[0];
  //extracting edgelist
  Rcpp::DataFrame edgelist = Rcpp::as<Rcpp::DataFrame>(stats[1]);
  //extracting riskset
  Rcpp::DataFrame riskset = Rcpp::as<Rcpp::DataFrame>(stats[2]);
  //extracting evls object
  arma::mat evls = stats[5];
  //dyad indicators
  arma::vec ind = evls.col(0);
  //vector of event times
  arma::vec time = evls.col(1);
  //number of events
  int M = edgelist.nrows();
  //number of dyads
  int D = riskset.nrows();
  //including a zero at the beginning of event time vector
  arma::vec time_zero(M + 1);
  for(int i = 0; i < M; ++i){
    time_zero(i+1) = time(i);
  }
  //computing intervent time
  for(int i = 0; i < M; ++i){
    time(i) = log(time_zero(i+1) - time_zero(i));
  }
  
  //creating vector to store repeated inter-event times
  arma::vec times;
  //creating vector to store observed dyads indicator
  arma::vec indicator;
  //creating matrix to reshape the statistics array
  arma::mat st;
  //looping through events
  for(int i = 0; i < M; ++i){
    //repeating inter-event times
    arma::vec t = rep(time(i), D);
    times = join_cols(times, t);
    //creating dyads indicator
    arma::vec ind2(D, fill::zeros);
    ind2(ind(i)-1) = 1;
    indicator = join_cols(indicator, ind2);
    //creating an auxiliary matrix to store the slice belonging to each event
    arma::mat slc = statistics.row(i);
    st = join_cols(st, slc);
  }
  
  //creating final output
  arma::mat final = join_rows(indicator, times, st);
  
  return(mat2DF(final));
  
}

//' burnInRows
//'
//' A function that applies row-wise burn-in to the MCMC samples
//' 
//' @param Nthin number of samples to skip before keeping one MCMC sample
//' @param Nburn number of samples to discard at the beginning of the chain
//' @param samples a 3-d array containing MCMC samples
//' 
//' @return 3-d array with altered dimension, leaving out Nburn samples and keeping one out of every Nthin
//'
// [[Rcpp::export]]
arma::cube burnInRows(double Nthin, double Nburn, const arma::cube& samples){
  
  double rows = (samples.n_rows - Nburn)/Nthin;
  
  arma::cube final(rows, samples.n_cols, samples.n_slices);
  
  for(arma::uword i = 0; i < final.n_slices; ++i){
    
    int num = Nburn;
    
    for(arma::uword j = 0; j < final.n_rows; ++j){
      
      final.slice(i).row(j) = samples.slice(i).row(num);
      
      num += Nthin;
      
    }
    
  }
  
  return(final);
  
}

//' burnInSlices
//' 
//' A function that applies slice-wise burn-in to the MCMC samples
//' 
//' @param Nthin number of samples to skip before keeping one MCMC sample
//' @param Nburn number of samples to discard at the beginning of the chain
//' @param samples a 3-d array containing MCMC samples
//' 
//' @return 3-d array with altered dimension, leaving out Nburn samples and keeping one out of every Nthin
//'
// [[Rcpp::export]]
arma::cube burnInSlices(double Nthin, double Nburn, const arma::cube& samples){
  
  double slices = (samples.n_slices - Nburn)/Nthin;
  
  arma::cube final(samples.n_rows, samples.n_cols, slices);
  
  int num = Nburn;
  
  for(arma::uword i = 0; i < slices; ++i){
      
    final.slice(i) = samples.slice(num);
      
    num += Nthin;
    
  }
  
  return(final);
  
}


////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

//' UpdateBeta
//' 
//' A function to update the group specific parameters of the multilevel model
//' 
//' @param p an integer that determines the dimension of the vector of parameters in each group
//' @param K an integer representing the number of groups in the analysis
//' @param beta_hat a matrix with dimensions p x K representing the MLE estimates for each independent groups
//' @param omega_hat a 3-d array with dimensions p x p x K representing the covariance matrix of the MLEs
//' @param Mu a vector with dimension p representing the random-effect means
//' @param Sigma a matrix p x p representing the covariance matrix of the random effects
//' @param omega a vector with dimension K representing the scale mixture variable that allows beta to have a student-t likelihood
//' 
//' @return a matrix with dimensions p x K representing one sample of the posterior distribution of the beta
//'
// [[Rcpp::export]]
arma::mat UpdateBeta(int p, int K, const arma::mat& beta_hat, const arma::cube& omega_hat,
                     const arma::vec& Mu, const arma::mat& Sigma, const arma::vec& omega){
  
  //output
  arma::mat output(p, K);
  
  //posterior mean
  arma::vec post_mean(p);
  
  //posterior covariance matrix
  arma::mat post_var(p, p);
  
  for(int i = 0; i < K; ++i){
      
      post_var = inv(inv(Sigma) + (1/omega(i) * inv(omega_hat.slice(i))));
      
      post_mean = post_var * ((1/omega(i) * inv(omega_hat.slice(i))) * beta_hat.col(i) + inv(Sigma) * Mu); 
      
      output.col(i) = rmvnorm(1, post_mean, post_var).t();
    
  }
  
  return(output);
  
}



//' UpdateOmega
//'
//' A function to update the scale mixture variable that allows beta to have a student-t likelihood
//'
//' @param p an integer that determines the dimension of the vector of parameters in each group
//' @param K an integer representing the number of groups in the analysis
//' @param beta_hat a matrix with dimensions p x K representing the MLE estimates for each independent groups
//' @param omega_hat a 3-d array with dimensions p x p x K representing the covariance matrix of the MLEs
//' @param Mu a vector with dimension p representing the random-effect means
//' @param Sigma a matrix p x p representing the covariance matrix of the random effects
//' @param nu a vector of dimension K representing the degress of freedom of the student-t distribution for each group
//' @param beta a matrix with dimensions p x K representing one MCMC sample of the group-specific parameters, aka random effects
//' 
//' @return a vector with dimension K representing one MCMC sample of the scale mixture variable omega
//'
// [[Rcpp::export]]
arma::vec UpdateOmega(int p, int K, const arma::mat& beta_hat, const arma::cube& omega_hat,
                      const arma::vec& Mu, const arma::mat& Sigma, const arma::vec& nu, 
                      const arma::mat& beta){
  
  //output
  arma::vec output(K);
  
  for(int i = 0; i < K; ++i){
    
    output(i) = 1/R::rgamma((nu(i) + p)/2, 0.5 * sum(((beta_hat.col(i) - beta.col(i)).t() * inv(omega_hat.slice(i)) * (beta_hat.col(i) - beta.col(i)))) + nu(i)/2);
    
  }
  
  return(output);
  
}


//' UpdateMu
//' 
//' A function to update the vector of random-effect means
//' @param p an integer representing the number of random effects in the model
//' @param K is the number of groups in the multilevel analysis
//' @param TauMu is a matrix with dimensions p x p representing the covariance matrix of the prior of Mu
//' @param Sigma is a p x p matrix representing the covariance matrix of the random effects
//' @param beta is a p x K matrix representing one MCMC sample of the p random effects for each of the K groups
//' 
//' @return a matrix with dimensions 1 x p representing one MCMC sample of the random-effect mean vector
//'
// [[Rcpp::export]]
arma::mat UpdateMu(int p, int K, const arma::mat& TauMu, const arma::mat& Sigma,
                   const arma::mat& beta){
  
  arma::rowvec beta_sum(p);
  
  for(int i = 0; i < p; ++i){
    
    beta_sum(i) = sum(beta.row(i)); //suming the rows of beta matrix
    
  }
  
  arma::mat variance = inv(inv(TauMu) + K * inv(Sigma)); //posterior cov mat
  
  arma::vec mean = variance * (beta_sum * inv(Sigma)).t(); //posterior mean vec
  
  arma::mat mu = rmvnorm(1, mean, variance); //generating posterior samples
  
  return(mu);
  
}

//' UpdateSigma
//' 
//' A function to update the covariance matrix of the random effects, for details see Huang & Wand (2013)
//' 
//' @param eta an integer representing a prior hyperparameter, this will allow Sigma to have a Half-t prior
//' @param K an integer representing the number of groups
//' @param p an integer representing the number of the covariates in the models, Sigma in the end will have dimensions p x p
//' @param Lambda a vector with dimension p, this is an auxiliary variable that allows Sigma to have a Half-t prior
//' @param beta is a matrix with dimensions p x K representing one MCMC sample of the random effects
//' @param Mu is a vector with dimension p representing the mean of the random effects
//' 
//' @return a p x p matrix representing one MCMC sample of the covariance matrix of the random effects
//'
// [[Rcpp::export]]
arma::mat UpdateSigma(int eta, int K, int p, const arma::vec& Lambda,
                      const arma::mat& beta, const arma::vec & Mu){
  
  arma::mat sqsum(p,p, fill::zeros); //storing the sum of squares
  
  for(int i = 0; i < K; ++i){
    
    sqsum += (beta.col(i) - Mu) * (beta.col(i) - Mu).t();
    
  }
  
  arma::mat sig = riwish(eta + K + p - 1, (2 * eta * diagmat(1/Lambda)) + sqsum);
  
  return(sig);
  
}

//' UpdateLambda
//' 
//' A function to update the mixture parameters that allow Sigma to have a Halt-t prior, see Huang & Wand (2013)
//' 
//' @param p an integer denoting the dimension of the covariance matrix Sigma
//' @param eta an integer that we set equal to 2, so the covariances have an uniform prior, see Huang & Wand (2013)
//' @param Sigma a p x p matrix representing one MCMC sample of the random-effect covariance matrix
//' @param xi a scalar that can be an arbitrarily large number, see Gelman (2006)
//' 
//' @return a row vector representing one MCMC sample of the mixture parameter that allows the covariance matrix of the random effects to have a Half-t prior
//'
// [[Rcpp::export]]
arma::rowvec UpdateLambda(int p, int eta, const arma::mat& Sigma, double xi){
  
  arma::mat invSigma = inv(Sigma);
  
  int shape = (eta + p)/2;
  
  arma::vec diagonal = invSigma.diag();
  
  arma::vec rate(p);
  
  for(int j = 0; j < p; ++j){
    
    rate(j) = eta * diagonal(j) + 1/pow(xi, 2);
    
  }
  
  arma::rowvec lamb(p);
  
  for(int i = 0; i < p; ++i){
    
    lamb(i) = 1/R::rgamma(shape, rate(i));
    
  }
  
  return(lamb);
  
}



//' UpdatePsi
//' 
//' A function that updates the fixed-effect parameters in the mixed-effect relational event model 
//'
//' @param q an integer representing the number of fixed effects in the model
//' @param K an integer representing the number of groups (or networks) in the multilevel analysis
//' @param TauPsi a vector of dimension q, which represents the prior variance of the fixed effects
//' @param SigmaPsi a matrix with dimensions q x K representing the standard errors of the MLEs for the fixed effects
//' @param psi a matrix with dimensions q x K representing the MLEs for the fixed effects independently estimated from every group
//' 
//' @return a row vector with dimension q representing one MCMC sample of the fixed effects
//'
// [[Rcpp::export]]
arma::rowvec UpdatePsi(int q, int K, const arma::vec& TauPsi, const arma::mat& SigmaPsi, 
                       const arma::mat& psi){
  
  arma::rowvec fixed(q); //store samples in the end
  
  for(int i = 0; i < q; ++i){
    
    double variance = 1/(sum(1/SigmaPsi.row(i)) + 1/TauPsi(i));
    
    double mean = sum(psi.row(i)/SigmaPsi.row(i)) * variance;
    
    fixed(i) = R::rnorm(mean, sqrt(variance)); //sampling
    
  }
  
  return(fixed);
  
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////

//' sampler
//' 
//' A function that runs the Gibbs sampler for a tie-oriented mixed-effect relational event model
//' 
//' @param Niter number of iterations to run the Gibbs sampler
//' @param Nchain number of MCMC chains
//' @param NburnIn number of samples to discard at the beginning of every chain
//' @param Nthin number of samples to skip for every MCMC sample kept 
//' @param p number of random effects in the model
//' @param q number of fixed effects in the model
//' @param K number of groups (or networks) in the multilevel analysis
//' @param nu a vector with dimension K, representing the degrees of freedom of the student-t distribution, it is set to (number of events in every group) - p
//' @param beta a 3-d array with dimension p x K x Nchain containing the starting values for the random effects parameters
//' @param beta_hat a matrix with dimensions p X K containing the MLEs of the random effects estimated from the independent networks
//' @param omega_hat a 3-d array with dimensions p x p x K containing the covariance matrices of the MLEs estimated from the independent groups
//' @param Mu a matrix with dimensions Nchain x p containing the initial values for the random-effect means
//' @param Sigma a 3-d array with dimensions p x p x Nchain containing the initial values for the random-effect covariance matrix
//' @param eta an integer that we set equal to 2, so the covariances have an uniform prior, see Huang & Wand (2013)
//' @param Lambda a matrix with dimensions Nchain x p containing initial values for the Lambda parameter
//' @param xi a scalar that can be an arbitrarily large number, see Gelman (2006)
//' @param TauMu is a matrix with dimensions p x p representing the covariance matrix of the prior of Mu
//' @param TauPsi a vector of dimension q, which represents the prior variance of the fixed effects 
//' @param SigmaPsi a matrix with dimensions q x K representing the standard errors of the MLEs for the fixed effects
//' @param omega a matrix with dimensions p x Nchain containing the initial values for the scale mixture parameters that allow the random-effects to have a student-t likelihood
//' @param psi a matrix with dimensions q x K representing the MLEs for the fixed effects independently estimated from every group
//' @param randomModel a boolean, if equals to TRUE it states that the model contain only random effects
//' @param fixedModel a boolean, if equals TRUE it states that the model contain only fixed effects
//' 
//' @return a list containing the MCMC samples for all parameters in the model
//'
// [[Rcpp::export]]
Rcpp::List sampler(const int Niter, const int Nchain, const int NburnIn, const int Nthin, const int p, const int q, const int K,
                   const arma::vec& nu, const arma::cube beta, const arma::mat& beta_hat, const arma::cube& omega_hat,
                   const arma::mat& Mu, const arma::cube Sigma, double eta, const arma::mat& Lambda, const double& xi,
                   const arma::mat& TauMu, const arma::vec& TauPsi, const arma::mat& SigmaPsi, const arma::mat& omega,
                   const arma::mat& psi, bool randomModel, bool fixedModel){
  
  //Niter: is the number of MCMC iterations
  //Nchain: is the number of Markov chains you want to run
  //NburnIn: is the number of samples you want to discard as burn-in phase
  //p: is the number of covariates in the model
  //K: is the number of clusters or groups in the multilevel model
  //nu: is a vector of degrees of freedom for the multivariate t distribution (use number of events - p)
  //beta: is a p x K matrix with starting values for the regression coefficients of each cluster
  //beta_hat: also a p x K matrix with MLE's obtained from remstimate
  //omega_hat: is a cube of p x p x K cov matrices for the MLE's obtained from remstimate
  //Mu: initial values for the grand mean vector in the hierarchical prior
  //Sigma: initial value for covariance matrix in the hierarchial prior
  //eta: is the prior hyperparameter for the covariance matrix, it is the degrees of freedom of the Inverse Wishart distribution
  //Lambda: prior hyperparameter for the covariance matrix, scale matrix of the Inverse Wishart distribution
  //TauMu: it is the prior covariance matrix for the grand mean vector Mu
  //TauPsi: is the standard error of the covariates that will be treated as fixed-effects
  //SigmaPsi: is the prior covariance matrix of the fixed effects
  //omega: is the scale-mixture variable
  //psi: is the point estimates of the fixed-effects
  //randomModel: is a boolean, if TRUE the model is a random-effects model
  
  Rcpp::List betaOut(Nchain); //store beta samples for every cluster
  
  arma::cube muOut(Niter, p, Nchain); //store grand mean samples
  
  arma::cube psiOut(Niter, q, Nchain); //store fixed effect samples 

  Rcpp::List sigmaOut(Nchain); //store cov matrix samples
  
  arma::cube lambdaOut(Niter, p, Nchain); //store scale mixture for the cov matrix
  
  arma::cube omegaOut(Niter, K, Nchain); //store scale mixture for the likelihood
  
  /////////////////////////////////////////////////////////////////////////////
  
  //Now we get into the sampler algorithm
  
  for(int i = 0; i < Nchain; ++i){
    
    arma::mat betaSample(p, K); //this an intermediate list to store the output of the posterior for the group-specific effects
    arma::cube betaCube(Niter, p, K); //This is a cube that stores beta samples for each cluster in its slices
    arma::cube sigmaCube(p, p, Niter); //Stores covariance matrix in each slice
    
    for(int j = 0; j < Niter; ++j){
      if(!fixedModel){
        if(j == 0){ //using initial values
          
          betaSample = UpdateBeta(p, K, beta_hat, omega_hat, Mu.col(i), Sigma.slice(i), omega.col(i));
          muOut.slice(i).row(j) = UpdateMu(p, K, TauMu, Sigma.slice(i), beta.slice(i));
          sigmaCube.slice(j) = UpdateSigma(eta, K, p, Lambda.col(i), beta.slice(i), Mu.col(i));
          lambdaOut.slice(i).row(j) = UpdateLambda(p, eta, Sigma.slice(i), xi);
          omegaOut.slice(i).row(j) = UpdateOmega(p, K, beta_hat, omega_hat, Mu.col(i), Sigma.slice(i), nu, beta.slice(i)).t();
          
        } else {
          
          betaSample = UpdateBeta(p, K, beta_hat, omega_hat, muOut.slice(i).row(j-1).t(), sigmaCube.slice(j-1), omegaOut.slice(i).row(j-1).t());
          muOut.slice(i).row(j) = UpdateMu(p, K, TauMu, sigmaCube.slice(j-1), betaSample);
          sigmaCube.slice(j) = UpdateSigma(eta, K, p, lambdaOut.slice(i).row(j-1).t(), betaSample, muOut.slice(i).row(j-1).t());
          lambdaOut.slice(i).row(j) = UpdateLambda(p, eta, sigmaCube.slice(j-1), xi);
          omegaOut.slice(i).row(j) = UpdateOmega(p, K, beta_hat, omega_hat, muOut.slice(i).row(j-1).t(), sigmaCube.slice(j-1), nu, betaSample).t();
          
        }
      }
      
      if(!randomModel){
        psiOut.slice(i).row(j) = UpdatePsi(q, K, TauPsi, SigmaPsi, psi); //sampling fixed effects
      }
      
      if(!fixedModel){
        betaSample = betaSample.t();
        for(int l = 0; l < K; ++l){
          betaCube.slice(l).row(j) = betaSample.row(l);
        }
      } 
      
    }
    if(!fixedModel){
      betaOut[i] = burnInRows(Nthin, NburnIn, betaCube); //already burning in
      sigmaOut[i] = burnInSlices(Nthin, NburnIn, sigmaCube); //already burning in 
    }
    
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  
  if(!fixedModel){
    muOut = burnInRows(Nthin, NburnIn, muOut); //burning in the grand mean
    lambdaOut = burnInRows(Nthin, NburnIn, lambdaOut); //burning in the scale mixture of the cov matrix
    omegaOut = burnInRows(Nthin, NburnIn, omegaOut); //burning in the scale mixture of the likelihood
  }

  if(!randomModel){  
  psiOut = burnInRows(Nthin, NburnIn, psiOut); //burning in the fixed effects 
  }
  
    
  if(!randomModel){//Mixed-effect model
    if(!fixedModel){//Fixed-effect model
      return(Rcpp::List::create(Rcpp::Named("beta") = betaOut, //random effects
                                Rcpp::Named("mu") = muOut, //random-effect means
                                Rcpp::Named("sigma") = sigmaOut, //random-effect covariance matrix
                                Rcpp::Named("alpha") = lambdaOut, //mixture variable of cov-matrix prior
                                Rcpp::Named("psi") = psiOut, //fixed effects
                                Rcpp::Named("omega") = omegaOut)); //scale mixture of random-effects likelihood
    } else {
      return(Rcpp::List::create(Rcpp::Named("psi") = psiOut)); 
    }
  } else {//Random-effect model
    return(Rcpp::List::create(Rcpp::Named("beta") = betaOut,
                              Rcpp::Named("mu") = muOut,
                              Rcpp::Named("sigma") = sigmaOut,
                              Rcpp::Named("alpha") = lambdaOut,
                              Rcpp::Named("omega") = omegaOut));
  }

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//' samplerActor
//' 
//' A function that runs the Gibbs sampler for the DyNaM (actor-oriented) mixed-effect relational event model
//' 
//' @param Niter number of iterations to run the Gibbs sampler
//' @param Nchain number of MCMC chains
//' @param NburnIn number of samples to discard at the beginning of every chain
//' @param Nthin number of samples to skip for every MCMC sample kept 
//' @param v an integer representing the number of random effects in the receiver model
//' @param u an integer representing the number of fixed effects in the receiver model 
//' @param K number of groups (or networks) in the multilevel analysis
//' @param nu a vector with dimension K, representing the degrees of freedom of the student-t distribution, it is set to (number of events in every group) - p
//' @param beta a 3-d array with dimension v x K x Nchain containing the starting values for the random effects parameters in the receiver model
//' @param beta_hat a matrix with dimensions v X K containing the MLEs of the random effects in the receiver model estimated from the independent networks
//' @param omega_hat a 3-d array with dimensions v x v x K containing the covariance matrices of the MLEs for the random effects in the receiver model estimated from the independent groups
//' @param MuB a matrix with dimensions v x Nchain containing the initial values of the random-effect means in the receiver model
//' @param Sigma a 3-d array with dimensions v x v x Nchain containing the initial values for the covariance matrix of the random-effects in the receiver model
//' @param eta an integer that we set equal to 2, so the covariances have an uniform prior, see Huang & Wand (2013)
//' @param LambdaB a matrix with dimensions Nchain x v containing initial values for the Lambda parameter in the receiver model
//' @param xi a scalar that can be an arbitrarily large number, see Gelman (2006)
//' @param TauMuB is a matrix with dimensions v x v representing the covariance matrix of the prior of Mu in the receiver model
//' @param TauPsi a vector of dimension u, which represents the prior variance of the fixed effects in the receiver model
//' @param SigmaPsi a matrix with dimensions u x K representing the standard errors of the MLEs for the fixed effects in the receiver model
//' @param omegaB a matrix with dimensions v x Nchain containing the initial values for the scale mixture parameters that allow the random-effects in the receiver model to have a student-t likelihood
//' @param psi a matrix with dimensions u x K representing the MLEs for the fixed effects in the receiver model independently estimated from every group
//' @param p an integer representing the number of random effects in the sender model
//' @param q an integer representing the number of fixed effects in the sender model
//' @param gamma a 3-d array with dimension p x K x Nchain containing the starting values for the random effects parameters in the sender
//' @param gamma_hat a matrix with dimensions p X K containing the MLEs of the random effects in the sender model estimated from the independent networks
//' @param zeta_hat a 3-d array with dimensions p x p x K containing the covariance matrices of the MLEs for the random effects in the sender model estimated from the independent groups
//' @param MuG a matrix with dimensions p x Nchain containing the initial values of the random-effect means in the sender model
//' @param Zeta a 3-d array with dimensions p x p x Nchain containing the initial values for the covariance matrix of the random-effects in the sender model
//' @param LambdaG a matrix with dimensions Nchain x p containing initial values for the Lambda parameter in the sender model
//' @param TauMuG is a matrix with dimensions p x p representing the covariance matrix of the prior of Mu in the sender model
//' @param TauPhi a vector of dimension q, which represents the prior variance of the fixed effects in the sender model
//' @param SigmaPhi a matrix with dimensions q x K representing the standard errors of the MLEs for hte fixed effects in the sender model
//' @param omegaG matrix with dimensions p x Nchain containing the initial values for the scale mixture parameters that allow the random-effects in the receiver model to have a student-t likelihood
//' @param phi a matrix with dimensions q x K representing the MLEs for the fixed effects in the sender model independently estimated from every group
//' @param randomModelSnd a boolean, if equals to TRUE it states that the sender model contain only random effects
//' @param fixModelSnd a boolean, if equals TRUE it states that the sender model contain only fixed effects
//' @param randomModelRec a boolean, if equals to TRUE it states that the receiver model contain only random effects
//' @param fixModelRec a boolean, if equals TRUE it states that the receiver model contain only fixed effects
//' 
//' @return a list containing two other lists with the MCMC sample for all parameters in the sender and receiver models
//'
// [[Rcpp::export]]
Rcpp::List samplerActor(const int Niter, const int Nchain, const int NburnIn, const int Nthin, const int v, const int u, const int K,
                        const arma::vec& nu, const arma::cube beta, const arma::mat& beta_hat, const arma::cube& omega_hat,
                        const arma::mat& MuB, const arma::cube Sigma, double eta, const arma::mat& LambdaB, const double& xi,
                        const arma::mat& TauMuB, const arma::vec& TauPsi, const arma::mat& SigmaPsi, const arma::mat& omegaB,
                        const arma::mat& psi, const int p, const int q, const arma::cube gamma, 
                        const arma::mat& gamma_hat, const arma::cube& zeta_hat, const arma::mat& MuG, const arma::cube Zeta,
                        const arma::mat& LambdaG, const arma::mat& TauMuG, const arma::vec& TauPhi,
                        const arma::mat& SigmaPhi, const arma::mat& omegaG, const arma::mat& phi,
                        bool randomModelSnd, bool randomModelRec, bool fixModelSnd, bool fixModelRec){
  
  //Niter: is the number of MCMC iterations
  //Nchain: is the number of Markov chains you want to run
  //NburnIn: is the number of samples you want to discard as burn-in phase
  //p: is the number of covariates in the model
  //K: is the number of clusters or groups in the multilevel model
  //nu: is a vector of degrees of freedom for the multivariate t distribution (use number of events - p)
  //beta: is a p x K matrix with starting values for the regression coefficients of each cluster
  //beta_hat: also a p x K matrix with MLE's obtained from remstimate
  //omega_hat: is a cube of p x p x K cov matrices for the MLE's obtained from remstimate
  //Mu: initial values for the grand mean vector in the hierarchical prior
  //Sigma: initial value for covariance matrix in the hierarchial prior
  //eta: is the prior hyperparameter for the covariance matrix, it is the degrees of freedom of the Inverse Wishart distribution
  //Lambda: prior hyperparameter for the covariance matrix, scale matrix of the Inverse Wishart distribution
  //TauMu: it is the prior covariance matrix for the grand mean vector Mu
  //TauPsi: is the standard error of the covariates that will be treated as fixed-effects
  //SigmaPsi: is the prior covariance matrix of the fixed effects
  //omega: is the scale-mixture variable
  //psi: is the point estimates of the fixed-effects
  //randomModel: is a boolean, if TRUE the model is a random-effects model
  
  //RECEIVER MODEL
  Rcpp::List betaOut(Nchain); //store beta samples for every cluster
  
  arma::cube muBOut(Niter, v, Nchain); //store grand mean samples
  
  arma::cube psiOut(Niter, u, Nchain); //store fixed effect samples 
  
  Rcpp::List sigmaOut(Nchain); //store cov matrix samples
  
  arma::cube lambdaBOut(Niter, v, Nchain); //store scale mixture for the cov matrix
  
  arma::cube omegaBOut(Niter, K, Nchain); //store scale mixture for the likelihood
  
  //SENDER MODEL
  Rcpp::List gammaOut(Nchain); //store beta samples for every cluster
  
  arma::cube muGOut(Niter, p, Nchain); //store grand mean samples
  
  arma::cube phiOut(Niter, q, Nchain); //store fixed effect samples 
  
  Rcpp::List zetaOut(Nchain); //store cov matrix samples
  
  arma::cube lambdaGOut(Niter, p, Nchain); //store scale mixture for the cov matrix
  
  arma::cube omegaGOut(Niter, K, Nchain); //store scale mixture for the likelihood
  
  
  /////////////////////////////////////////////////////////////////////////////
  
  //Now we get into the sampler algorithm
  
  for(int i = 0; i < Nchain; ++i){
    
    //RECEIVER MODEL
    arma::mat betaSample(v, K); //this an intermediate list to store the output of the posterior for the group-specific effects
    
    arma::cube betaCube(Niter, v, K); //This is a cube that stores beta samples for each cluster in its slices
    
    arma::cube sigmaCube(v, v, Niter); //Stores covariance matrix in each slice
    
    //SENDER MODEL
    arma::mat gammaSample(p, K); //this an intermediate list to store the output of the posterior for the group-specific effects
    
    arma::cube gammaCube(Niter, p, K); //This is a cube that stores beta samples for each cluster in its slices
    
    arma::cube zetaCube(p, p, Niter); //Stores covariance matrix in each slice
    
    for(int j = 0; j < Niter; ++j){
      
      if(j == 0){ //using initial values
        
        if(!fixModelRec){
          //RECEIVER MODEL
          betaSample = UpdateBeta(v, K, beta_hat, omega_hat, MuB.col(i), Sigma.slice(i), omegaB.col(i));
          muBOut.slice(i).row(j) = UpdateMu(v, K, TauMuB, Sigma.slice(i), beta.slice(i));
          sigmaCube.slice(j) = UpdateSigma(eta, K, v, LambdaB.col(i), beta.slice(i), MuB.col(i));
          lambdaBOut.slice(i).row(j) = UpdateLambda(v, eta, Sigma.slice(i), xi);
          omegaBOut.slice(i).row(j) = UpdateOmega(v, K, beta_hat, omega_hat, MuB.col(i), Sigma.slice(i), nu, beta.slice(i)).t();
        }
        if(!fixModelSnd){
          //SENDER MODEL
          gammaSample = UpdateBeta(p, K, gamma_hat, zeta_hat, MuG.col(i), Zeta.slice(i), omegaG.col(i));
          muGOut.slice(i).row(j) = UpdateMu(p, K, TauMuG, Zeta.slice(i), gamma.slice(i));
          zetaCube.slice(j) = UpdateSigma(eta, K, p, LambdaG.col(i), gamma.slice(i), MuG.col(i));
          lambdaGOut.slice(i).row(j) = UpdateLambda(p, eta, Zeta.slice(i), xi);
          omegaGOut.slice(i).row(j) = UpdateOmega(p, K, gamma_hat, zeta_hat, MuG.col(i), Zeta.slice(i), nu, gamma.slice(i)).t();
        }
      } else {
        if(!fixModelRec){
          //RECEIVER MODEL
          betaSample = UpdateBeta(v, K, beta_hat, omega_hat, muBOut.slice(i).row(j-1).t(), sigmaCube.slice(j-1), omegaBOut.slice(i).row(j-1).t());
          muBOut.slice(i).row(j) = UpdateMu(v, K, TauMuB, sigmaCube.slice(j-1), betaSample);
          sigmaCube.slice(j) = UpdateSigma(eta, K, v, lambdaBOut.slice(i).row(j-1).t(), betaSample, muBOut.slice(i).row(j-1).t());
          lambdaBOut.slice(i).row(j) = UpdateLambda(v, eta, sigmaCube.slice(j-1), xi);
          omegaBOut.slice(i).row(j) = UpdateOmega(v, K, beta_hat, omega_hat, muBOut.slice(i).row(j-1).t(), sigmaCube.slice(j-1), nu, betaSample).t();
        }
        if(!fixModelSnd){
          //SENDER MODEL
          gammaSample = UpdateBeta(p, K, gamma_hat, zeta_hat, muGOut.slice(i).row(j-1).t(), zetaCube.slice(j-1), omegaGOut.slice(i).row(j-1).t());
          muGOut.slice(i).row(j) = UpdateMu(p, K, TauMuG, zetaCube.slice(j-1), gammaSample);
          zetaCube.slice(j) = UpdateSigma(eta, K, p, lambdaGOut.slice(i).row(j-1).t(), gammaSample, muGOut.slice(i).row(j-1).t());
          lambdaGOut.slice(i).row(j) = UpdateLambda(p, eta, zetaCube.slice(j-1), xi);
          omegaGOut.slice(i).row(j) = UpdateOmega(p, K, gamma_hat, zeta_hat, muGOut.slice(i).row(j-1).t(), zetaCube.slice(j-1), nu, gammaSample).t();
        }
      }
      
      if(!randomModelRec){//RECEIVER MODEL
        psiOut.slice(i).row(j) = UpdatePsi(u, K, TauPsi, SigmaPsi, psi); //sampling fixed effects
      }
      
      if(!randomModelSnd){//SENDER MODEL
        phiOut.slice(i).row(j) = UpdatePsi(q, K, TauPhi, SigmaPhi, phi); //sampling fixed effects
      }
      
      if(!fixModelRec){
        betaSample = betaSample.t(); //RECEIVER
        for(int l = 0; l < K; ++l){
          betaCube.slice(l).row(j) = betaSample.row(l); //RECEIVER
        }
      }
      
      if(!fixModelSnd){
        gammaSample = gammaSample.t(); //SENDER 
        for(int l = 0; l < K; ++l){
          gammaCube.slice(l).row(j) = gammaSample.row(l); //SENDER
        }
      }
      
    }
    
    if(!fixModelRec){
      //RECEIVER MODEL
      betaOut[i] = burnInRows(Nthin, NburnIn, betaCube); //already burning in
      sigmaOut[i] = burnInSlices(Nthin, NburnIn, sigmaCube); //already burning in 
    }
    
    if(!fixModelSnd){
      //SENDER MODEL
      gammaOut[i] = burnInRows(Nthin, NburnIn, gammaCube); //already burning in
      zetaOut[i] = burnInSlices(Nthin, NburnIn, zetaCube); //already burning in 
    }
    
  }
  
  //////////////////////////////////////////////////////////////////////////////////
  
  if(!fixModelRec){//RECEIVER MODEL
    muBOut = burnInRows(Nthin, NburnIn, muBOut); //burning in the grand mean
    lambdaBOut = burnInRows(Nthin, NburnIn, lambdaBOut); //burning in the scale mixture of the cov matrix
    omegaBOut = burnInRows(Nthin, NburnIn, omegaBOut); //burning in the scale mixture of the likelihood
  }
  
  if(!fixModelSnd){//SENDER MODEL
    muGOut = burnInRows(Nthin, NburnIn, muGOut); //burning in the grand mean
    lambdaGOut = burnInRows(Nthin, NburnIn, lambdaGOut); //burning in the scale mixture of the cov matrix
    omegaGOut = burnInRows(Nthin, NburnIn, omegaGOut); //burning in the scale mixture of the likelihood
  }
  
  if(!randomModelRec){//RECEIVER  
    psiOut = burnInRows(Nthin, NburnIn, psiOut); //burning in the fixed effects 
  }
  
  if(!randomModelSnd){//SENDER  
    phiOut = burnInRows(Nthin, NburnIn, phiOut); //burning in the fixed effects 
  }
  
  ///OUTPUT OF THE FUNCTIONS
  
  Rcpp::List S1; //sender output
  Rcpp::List S2; //receiver output
  
  //Sender model contains fixed effects
  if(!randomModelSnd){ //model contains 
    if(!fixModelSnd){
      S1 = Rcpp::List::create(Rcpp::Named("gamma") = gammaOut,
                              Rcpp::Named("mu") = muGOut,
                              Rcpp::Named("sigma") = zetaOut,
                              Rcpp::Named("alpha") = lambdaGOut,
                              Rcpp::Named("phi") = phiOut,
                              Rcpp::Named("omega") = omegaGOut);
    } else {
      S1 = Rcpp::List::create(Rcpp::Named("phi") = phiOut);
    }
  } else {
    S1 = Rcpp::List::create(Rcpp::Named("gamma") = gammaOut,
                            Rcpp::Named("mu") = muGOut,
                            Rcpp::Named("sigma") = zetaOut,
                            Rcpp::Named("alpha") = lambdaGOut,
                            Rcpp::Named("omega") = omegaGOut);
  }
  //Receiver model contains fixed effects
  if(!randomModelRec){ //model contains 
    if(!fixModelRec){
      S2 = Rcpp::List::create(Rcpp::Named("beta") = betaOut, //random effects
                              Rcpp::Named("mu") = muBOut, //random-effect means
                              Rcpp::Named("sigma") = sigmaOut, //random-effect covariance matrix
                              Rcpp::Named("alpha") = lambdaBOut, //mixture variable of cov-matrix prior
                              Rcpp::Named("psi") = psiOut, //fixed effects
                              Rcpp::Named("omega") = omegaBOut); //scale mixture of random-effects likelihood
    } else {
      S2 = Rcpp::List::create(Rcpp::Named("psi") = psiOut);
    }
  } else {
    S2 = Rcpp::List::create(Rcpp::Named("beta") = betaOut, //random effects
                            Rcpp::Named("mu") = muBOut, //random-effect means
                            Rcpp::Named("sigma") = sigmaOut, //random-effect covariance matrix
                            Rcpp::Named("alpha") = lambdaBOut, //mixture variable of cov-matrix prior
                            Rcpp::Named("omega") = omegaBOut); //scale mixture of random-effects likelihood
  }
  
  return(Rcpp::List::create(Rcpp::Named("sender_rate") = S1,
                           Rcpp::Named("receiver_choice") = S2));
  
}