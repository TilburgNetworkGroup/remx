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
arma::cube burnInRows(double Nthin,
                      double Nburn,
                      const arma::cube& samples){

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
arma::cube burnInSlices(double Nthin,
                        double Nburn,
                        const arma::cube& samples){

  double slices = (samples.n_slices - Nburn)/Nthin;

  arma::cube final(samples.n_rows, samples.n_cols, slices);

  int num = Nburn;

  for(arma::uword i = 0; i < slices; ++i){

    final.slice(i) = samples.slice(num);

    num += Nthin;

  }

  return(final);

}

//' rndMeans
//' This function only returns the components of the fixed effects that are random effects means
//'@param x matrix specifying which parameters are to be extracted from the vector
//'@param psi a vector containing fixed- and random-effects estimates
//[[Rcpp::export]]
arma::vec rndMeans(const arma::mat& x,
                   const arma::vec& psi){

 arma::vec mu = x * psi;

 arma::uvec ind = find(mu != 0);

 arma::vec mean(ind.size());

 for(arma::uword i = 0; i < ind.size(); ++i){
   mean(i) = mu(ind(i));
 }

 return(mean);
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
//' @param Mu a vector p x 1 representing the mean of the random effects
//' @param Sigma a matrix p x p representing the covariance matrix of the random effects
//'
//' @return a matrix with dimensions p x K representing one sample of the posterior distribution of the beta
//'
// [[Rcpp::export]]
arma::mat UpdateBeta(int p,
                     int K,
                     const arma::mat& beta_hat,
                     const arma::cube& omega_hat,
                     const arma::vec& Mu,
                     const arma::mat& Sigma){

  //output
  arma::mat output(p, K);

  //posterior mean
  arma::vec post_mean(p);

  //posterior covariance matrix
  arma::mat post_var(p, p);

  for(int i = 0; i < K; ++i){

      post_var = inv(inv(Sigma) + inv(omega_hat.slice(i)));

      post_mean = post_var * (inv(omega_hat.slice(i)) * beta_hat.col(i) + inv(Sigma) * Mu);

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
arma::vec UpdateOmega(int p,
                      int K,
                      const arma::mat& beta_hat,
                      const arma::cube& omega_hat,
                      const arma::vec& Mu,
                      const arma::mat& Sigma,
                      const arma::vec& nu,
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
arma::mat UpdateMu(int p,
                   int K,
                   const arma::mat& TauMu,
                   const arma::mat& Sigma,
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
//'
//' @return a p x p matrix representing one MCMC sample of the covariance matrix of the random effects
//'
// [[Rcpp::export]]
arma::mat UpdateSigma(int eta,
                      int K,
                      int p,
                      const arma::vec& Lambda,
                      const arma::mat& beta){

  arma::mat sqsum(p,p, fill::zeros); //storing the sum of squares

  for(int i = 0; i < K; ++i){

    sqsum += (beta.col(i)) * (beta.col(i)).t();

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
arma::rowvec UpdateLambda(int p,
                          int eta,
                          const arma::mat& Sigma,
                          double xi){

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
arma::rowvec UpdatePsi(int q,
                       int K,
                       const arma::vec& TauPsi,
                       const arma::mat& SigmaPsi,
                       const arma::mat& psi){

  arma::rowvec fixed(q); //store samples in the end

  for(int i = 0; i < q; ++i){

    double variance = 1/(sum(1/SigmaPsi.row(i)) + 1/TauPsi(i));

    double mean = sum(psi.row(i)/SigmaPsi.row(i)) * variance;

    fixed(i) = R::rnorm(mean, sqrt(variance)); //sampling

  }

  return(fixed);

}

//'UpdateMean
//'This is a new function that models fixed-effects and random-effects means together so we do not ignore their dependency
//'This function already does what we want
//'
//' @param q number of fixed effects
//' @param p number of random effects
//' @param K number of networks
//' @param SigmaMean covariance matrices of the MLE estimates
//' @param psi MLE estimates
//' @param delta random-effects estimates
//' @param randomModel boolean, if TRUE model contains random effects
//[[Rcpp::export]]
arma::rowvec UpdateMean(int q,
                        int p,
                        int K,
                        const arma::cube& SigmaMean,
                        const arma::mat& psi,
                        const arma::mat& delta,
                        bool randomModel){

  arma::mat variance(q, q, fill::zeros);

  arma::vec mean(q, fill::zeros);

  //creating a vector to compensate the fact that delta needs to be binded with a vector of 0's
  arma::vec zero(q-p, fill::zeros);

  //creating another vector to join them
  arma::vec joint_vec;

  //in this loop we just add the pieces that need to be summed over all networks
  for(int i = 0; i < q; ++i){
    if(randomModel){
      joint_vec = join_vert(zero, delta.col(i));
    } else {
      joint_vec = zero;
    }
    //Making the summation for the posterior variance
    variance += inv(SigmaMean.slice(i));
    //Making the summation for the posterior mean
    mean += inv(SigmaMean.slice(i)) * (psi.col(i));
  }
  //computing posterior variance
  variance = inv(variance);
  //computing posterior mean
  mean = variance * (mean);

  return(rmvnorm(1, mean, variance));
}

//'UpdateDelta
//'This is the function that samples the random effects
//'
//'@param q number of fixed effects
//'@param p number of random effects
//'@param K number of networks
//'@param beta_hat MLE estimates
//'@param Mu mean of the random effects
//'@param omega_hat_11 covariance matrix decomposition of the fixed effects
//'@param omega_hat_21 covariance matrix decomposition of the covariance between fixed and random effects
//'@param omega_hat_22 covariance matrix decomposition of the random effects
//'@param psi_hat fixed effects
//'@param psi initial values or previous MCMC sample of fixed effects
//'@param Sigma covariance matrix of the random effects
//[[Rcpp::export]]
arma::mat UpdateDelta(int q,
                      int p,
                      int K,
                      const arma::mat& beta_hat,
                      const arma::vec& Mu,
                      const arma::cube& omega_hat_11,
                      const arma::cube& omega_hat_21,
                      const arma::cube& omega_hat_22,
                      const arma::mat& psi_hat,
                      const arma::vec& psi,
                      const arma::mat& Sigma){
  //output
  arma::mat delta(p, K);

  //Computing B_hat
  arma::mat B_k;

  //Computing S_hat
  arma::mat S_k;

  //posterior mean
  arma::vec mean;

  //posterior variance
  arma::mat variance;

  for(int i = 0; i < K; ++i){
    //Computing B_k, every network has a different value for this one
    B_k = (omega_hat_21.slice(i).t() * inv(omega_hat_11.slice(i))) * (psi - psi_hat.col(i));
    //Computing S_k
    S_k = inv(omega_hat_22.slice(i) - (omega_hat_21.slice(i).t() * inv(omega_hat_11.slice(i)) * omega_hat_21.slice(i)));
    //Posterior variance
    variance = inv(S_k + inv(Sigma));
    //Posterior mean
    mean = variance * (S_k * (beta_hat.col(i) - Mu + B_k));
    //sampling delta_k
    delta.col(i) = rmvnorm(1, mean, variance).t();
  }

  return(delta);
}

//'UpdateDeltaInd
//'
//'@param p number of random effects
//'@param K number of networks
//'@param beta_hat MLE estimates
//'@param omega_hat_22 covariance matrix of the random effects
//'@param Mu random effects mean
//'@param Sigma covariance matrix of the random effects
//[[Rcpp::export]]
arma::mat UpdateDeltaInd(int p,
                         int K,
                         const arma::mat& beta_hat,
                         const arma::cube& omega_hat_22,
                         const arma::vec& Mu,
                         const arma::mat Sigma){
  arma::mat delta(p,K);
  //posterior mean
  arma::vec mean;
  //posterior variance
  arma::mat variance;

  for(int i = 0; i < K; ++i){
    variance = inv(inv(omega_hat_22.slice(i)) + inv(Sigma));
    mean = variance * (inv(omega_hat_22.slice(i)) * (beta_hat.col(i) - Mu));
    delta.col(i) = rmvnorm(1, mean, variance).t();
  }
  return(delta);
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
//' @param delta a 3-d array with dimension p x K x Nchain containing the starting values for the random effects parameters
//' @param beta_hat a matrix with dimensions p X K containing the MLEs of the random effects estimated from the independent networks
//' @param omega_hat a 3-d array with dimensions p x p x K containing the covariance matrices of the MLEs estimated from the independent groups
//' @param omega_hat_11 covariance matrix decomposition of the fixed effects
//' @param omega_hat_21 covariance matrix decomposition of the covariance between fixed and random effects
//' @param omega_hat_22 covariance matrix decomposition of the random effects
//' @param Mu a matrix with dimensions Nchain x p containing the initial values for the random-effect means
//' @param Sigma a 3-d array with dimensions p x p x Nchain containing the initial values for the random-effect covariance matrix
//' @param eta an integer that we set equal to 2, so the covariances have an uniform prior, see Huang & Wand (2013)
//' @param Lambda a matrix with dimensions Nchain x p containing initial values for the Lambda parameter
//' @param xi a scalar that can be an arbitrarily large number, see Gelman (2006)
//' @param SigmaMean a matrix with dimensions q x K representing the standard errors of the MLEs for the fixed effects
//' @param psi a matrix with dimensions q x K representing the MLEs for the fixed effects independently estimated from every group
//' @param psi_hat a matrix with dimensions q x K representing the MLEs for the fixed effects independently estimated from every group
//' @param random_effect a matrix with 1's and 0's
//' @param fixed_effect a matrix with 1's and 0's
//' @param randomModel a boolean, if equals to TRUE it states that the model contain only random effects
//'
//' @return a list containing the MCMC samples for all parameters in the model
//'
// [[Rcpp::export]]
Rcpp::List sampler(const int Niter,
                  const int Nchain,
                  const int NburnIn,
                  const int Nthin,
                  const int p,
                  const int q,
                  const int K,
                  const arma::cube delta, //if model contain only fixed effects, delta should be a matrix with 0's
                  const arma::mat& beta_hat,
                  const arma::cube& omega_hat,
                  const arma::cube& omega_hat_11,
                  const arma::cube& omega_hat_21,
                  const arma::cube& omega_hat_22,
                  const arma::mat& Mu,
                  const arma::cube Sigma,
                  double eta,
                  const arma::mat& Lambda,
                  const double& xi,
                  const arma::cube& SigmaMean,
                  const arma::mat& psi,
                  const arma::mat& psi_hat,
                  const arma::mat& random_effect, //this will be 1's for random-effects means and 0 otherwise
                  const arma::mat& fixed_effect, //this will be 1's for random-effects means and 0 otherwise
                  bool randomModel){

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

  Rcpp::List deltaOut(Nchain); //store beta samples for every cluster

  arma::cube muOut(Niter, q, Nchain); //store fixed effect samples

  Rcpp::List sigmaOut(Nchain); //store cov matrix samples

  arma::cube lambdaOut(Niter, p, Nchain); //store scale mixture for the cov matrix


  /////////////////////////////////////////////////////////////////////////////

  //Now we get into the sampler algorithm

  for(int i = 0; i < Nchain; ++i){

   arma::mat deltaSample(p, K); //this an intermediate list to store the output of the posterior for the group-specific effects
   arma::cube deltaCube(Niter, p, K); //This is a cube that stores beta samples for each cluster in its slices
   arma::cube sigmaCube(p, p, Niter); //Stores covariance matrix in each slice
   arma::vec rndMean(p);
   arma::vec fixEff(q);
   arma::mat psiAdj;
   if(randomModel){
     psiAdj = fixed_effect * psi_hat;
   } else {
     psiAdj = psi_hat;
   }
   
   for(int j = 0; j < Niter; ++j){
     if(randomModel){
       if(j == 0){ //using initial values
         if(p == q){
           deltaSample = UpdateDeltaInd(p, K, beta_hat, omega_hat_22, Mu.col(i), Sigma.slice(i));
         } else {
           deltaSample = UpdateDelta(q, p, K, beta_hat, Mu.col(i), omega_hat_11, omega_hat_21, omega_hat_22, psiAdj, psi.col(i), Sigma.slice(i));
         }
         sigmaCube.slice(j) = UpdateSigma(eta, K, p, Lambda.col(i), delta.slice(i));
         lambdaOut.slice(i).row(j) = UpdateLambda(p, eta, Sigma.slice(i), xi);
         muOut.slice(i).row(j) = UpdateMean(q, p, K, SigmaMean, psi_hat, delta.slice(i), randomModel); //sampling fixed effects

       } else {
         if(p == q){
           deltaSample = UpdateDeltaInd(p, K, beta_hat, omega_hat_22, muOut.slice(i).row(j-1).t(), sigmaCube.slice(j-1));
         } else {
           rndMean = random_effect * muOut.slice(i).row(j-1).t();
           fixEff = fixed_effect * muOut.slice(i).row(j-1).t();

           deltaSample = UpdateDelta(q, p, K, beta_hat, rndMean, omega_hat_11, omega_hat_21, omega_hat_22, psiAdj, fixEff, sigmaCube.slice(j-1));
         }
         sigmaCube.slice(j) = UpdateSigma(eta, K, p, lambdaOut.slice(i).row(j-1).t(), deltaSample);
         lambdaOut.slice(i).row(j) = UpdateLambda(p, eta, sigmaCube.slice(j-1), xi);
         muOut.slice(i).row(j) = UpdateMean(q, p, K, SigmaMean, psi_hat, deltaSample, randomModel); //sampling fixed effects

       }
     }

    if(!randomModel){
       muOut.slice(i).row(j) = UpdateMean(q, p, K, SigmaMean, psi_hat, delta, randomModel); //sampling fixed effects
    }

     if(randomModel){
       for(int l = 0; l < K; ++l){
         deltaCube.slice(l).row(j) = deltaSample.col(l).t();
       }
     }

   }
   if(randomModel){
     deltaOut[i] = burnInRows(Nthin, NburnIn, deltaCube); //already burning in
     sigmaOut[i] = burnInSlices(Nthin, NburnIn, sigmaCube); //already burning in
   }

  }

  //////////////////////////////////////////////////////////////////////////////////

  if(randomModel){
   lambdaOut = burnInRows(Nthin, NburnIn, lambdaOut); //burning in the scale mixture of the cov matrix
  }

  muOut = burnInRows(Nthin, NburnIn, muOut); //burning in the fixed effects

  if(randomModel){//Mixed-effect model
    return(Rcpp::List::create(Rcpp::Named("delta") = deltaOut, //random effects
                              Rcpp::Named("sigma") = sigmaOut, //random-effect covariance matrix
                              Rcpp::Named("alpha") = lambdaOut, //mixture variable of cov-matrix prior
                              Rcpp::Named("mu") = muOut)); //fixed-effects
   } else {
     return(Rcpp::List::create(Rcpp::Named("mu") = muOut));
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
//' @param beta a 3-d array with dimension v x K x Nchain containing the starting values for the random effects parameters in the receiver model
//' @param beta_hat a matrix with dimensions v X K containing the MLEs of the random effects in the receiver model estimated from the independent networks
//' @param omega_hat a 3-d array with dimensions v x v x K containing the covariance matrices of the MLEs for the random effects in the receiver model estimated from the independent groups
//' @param omega_hat_11 covariance matrix decomposition of the fixed effects
//' @param omega_hat_21 covariance matrix decomposition of the covariance between fixed and random effects
//' @param omega_hat_22 covariance matrix decomposition of the random effects
//' @param MuB a matrix with dimensions v x Nchain containing the initial values of the random-effect means in the receiver model
//' @param Sigma a 3-d array with dimensions v x v x Nchain containing the initial values for the covariance matrix of the random-effects in the receiver model
//' @param eta an integer that we set equal to 2, so the covariances have an uniform prior, see Huang & Wand (2013)
//' @param LambdaB a matrix with dimensions Nchain x v containing initial values for the Lambda parameter in the receiver model
//' @param xi a scalar that can be an arbitrarily large number, see Gelman (2006)
//' @param SigmaPsi a matrix with dimensions u x K representing the standard errors of the MLEs for the fixed effects in the receiver model
//' @param psi a matrix with dimensions u x K representing the MLEs for the fixed effects in the receiver model independently estimated from every group
//' @param psi_hat a matrix with dimensions u x K representing the MLEs for the fixed effects in the receiver model independently estimated from every group
//' @param p an integer representing the number of random effects in the sender model
//' @param q an integer representing the number of fixed effects in the sender model
//' @param gamma a 3-d array with dimension p x K x Nchain containing the starting values for the random effects parameters in the sender
//' @param gamma_hat a matrix with dimensions p X K containing the MLEs of the random effects in the sender model estimated from the independent networks
//' @param zeta_hat a 3-d array with dimensions p x p x K containing the covariance matrices of the MLEs for the random effects in the sender model estimated from the independent groups
//' @param zeta_hat_11 covariance matrix decomposition of the fixed effects
//' @param zeta_hat_21 covariance matrix decomposition of the covariance between fixed and random effects
//' @param zeta_hat_22 covariance matrix decomposition of the random effects
//' @param MuG a matrix with dimensions p x Nchain containing the initial values of the random-effect means in the sender model
//' @param Zeta a 3-d array with dimensions p x p x Nchain containing the initial values for the covariance matrix of the random-effects in the sender model
//' @param LambdaG a matrix with dimensions Nchain x p containing initial values for the Lambda parameter in the sender model
//' @param SigmaPhi a matrix with dimensions q x K representing the standard errors of the MLEs for hte fixed effects in the sender model
//' @param phi a matrix with dimensions q x K representing the MLEs for the fixed effects in the sender model independently estimated from every group
//' @param phi_hat a matrix with dimensions q x K representing the MLEs for the fixed effects in the sender model independently estimated from every group
//' @param random_effect_rec matrix with 0's and 1's
//' @param fixed_effect_rec matrix with 0's and 1's
//' @param random_effect_snd matrix with 0's and 1's
//' @param fixed_effect_snd matrix with 0's and 1's
//' @param randomModelSnd a boolean, if equals to TRUE it states that the sender model contain only random effects
//' @param randomModelRec a boolean, if equals to TRUE it states that the receiver model contain only random effects
//'
//' @return a list containing two other lists with the MCMC sample for all parameters in the sender and receiver models
//'
// [[Rcpp::export]]
Rcpp::List samplerActor(const int Niter,
                        const int Nchain,
                        const int NburnIn,
                        const int Nthin,
                        const int v,
                        const int u,
                        const int K,
                        const arma::cube beta,
                        const arma::mat& beta_hat,
                        const arma::cube& omega_hat,
                        const arma::cube& omega_hat_11,
                        const arma::cube& omega_hat_21,
                        const arma::cube& omega_hat_22,
                        const arma::mat& MuB,
                        const arma::cube Sigma,
                        double eta,
                        const arma::mat& LambdaB,
                        const double& xi,
                        const arma::cube& SigmaPsi,
                        const arma::mat& psi,
                        const arma::mat& psi_hat,
                        const int p,
                        const int q,
                        const arma::cube gamma,
                        const arma::mat& gamma_hat,
                        const arma::cube& zeta_hat,
                        const arma::cube& zeta_hat_11,
                        const arma::cube& zeta_hat_21,
                        const arma::cube& zeta_hat_22,
                        const arma::mat& MuG,
                        const arma::cube Zeta,
                        const arma::mat& LambdaG,
                        const arma::cube& SigmaPhi,
                        const arma::mat& phi,
                        const arma::mat& phi_hat,
                        const arma::mat& random_effect_rec, //this will be 1's for random-effects means and 0 otherwise
                        const arma::mat& fixed_effect_rec,
                        const arma::mat& random_effect_snd, //this will be 1's for random-effects means and 0 otherwise
                        const arma::mat& fixed_effect_snd,
                        bool randomModelSnd,
                        bool randomModelRec){
  
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

  arma::cube muBOut(Niter, u, Nchain); //store grand mean samples

  Rcpp::List sigmaOut(Nchain); //store cov matrix samples

  arma::cube lambdaBOut(Niter, v, Nchain); //store scale mixture for the cov matrix

  //SENDER MODEL
  Rcpp::List gammaOut(Nchain); //store beta samples for every cluster

  arma::cube muGOut(Niter, q, Nchain); //store grand mean sample

  Rcpp::List zetaOut(Nchain); //store cov matrix samples

  arma::cube lambdaGOut(Niter, p, Nchain); //store scale mixture for the cov matrix


  /////////////////////////////////////////////////////////////////////////////

  //Now we get into the sampler algorithm

  for(int i = 0; i < Nchain; ++i){

   //RECEIVER MODEL
   arma::mat betaSample(v, K); //this an intermediate list to store the output of the posterior for the group-specific effects

   arma::cube betaCube(Niter, v, K); //This is a cube that stores beta samples for each cluster in its slices

   arma::cube sigmaCube(v, v, Niter); //Stores covariance matrix in each slice

   arma::vec rndMeanRec(p);
   arma::vec fixEffRec(q);
   arma::mat psiAdj;
   if(randomModelRec){
     psiAdj = fixed_effect_rec * psi_hat;
   } else {
     psiAdj = psi_hat;
   }

   //SENDER MODEL
   arma::mat gammaSample(p, K); //this an intermediate list to store the output of the posterior for the group-specific effects

   arma::cube gammaCube(Niter, p, K); //This is a cube that stores beta samples for each cluster in its slices

   arma::cube zetaCube(p, p, Niter); //Stores covariance matrix in each slice

   arma::vec rndMeanSnd(p);
   arma::vec fixEffSnd(q);
   arma::mat phiAdj;
   if(randomModelSnd){
     phiAdj = fixed_effect_snd * phi_hat;
   } else {
     phiAdj = phi_hat;
   }

   for(int j = 0; j < Niter; ++j){

     if(j == 0){ //using initial values

       if(randomModelRec){
         //RECEIVER MODEL
         if(u == v){
           betaSample = UpdateDeltaInd(v, K, beta_hat, omega_hat_22, MuB.col(i), Sigma.slice(i));
         } else {
           betaSample = UpdateDelta(u, v, K, beta_hat, MuB.col(i), omega_hat_11, omega_hat_21, omega_hat_22, psiAdj, psi.col(i), Sigma.slice(i));
         }
         muBOut.slice(i).row(j) = UpdateMean(u, v, K, SigmaPsi, psi_hat, beta.slice(i), randomModelRec);
         sigmaCube.slice(j) = UpdateSigma(eta, K, v, LambdaB.col(i), beta.slice(i));
         lambdaBOut.slice(i).row(j) = UpdateLambda(v, eta, Sigma.slice(i), xi);
       }
       if(randomModelSnd){
         //SENDER MODEL
         if(p == q){
           gammaSample = UpdateDeltaInd(p, K, gamma_hat, zeta_hat_22, MuG.col(i), Zeta.slice(i));
         } else {
           gammaSample = UpdateDelta(q, p, K, gamma_hat, MuG.col(i), zeta_hat_11, zeta_hat_21, zeta_hat_22, phiAdj, phi.col(i), Zeta.slice(i));
         }
         muGOut.slice(i).row(j) = UpdateMean(q, p, K, SigmaPhi, phi_hat, gamma.slice(i), randomModelSnd);
         zetaCube.slice(j) = UpdateSigma(eta, K, p, LambdaG.col(i), gamma.slice(i));
         lambdaGOut.slice(i).row(j) = UpdateLambda(p, eta, Zeta.slice(i), xi);
       }
     } else {
       if(randomModelRec){
           if(u == v){
           betaSample = UpdateDeltaInd(v, K, beta_hat, omega_hat_22, muBOut.slice(i).row(j-1).t(), sigmaCube.slice(j-1));
         } else {
           //separating the fixed-effect from the random-effect mean
           rndMeanRec = random_effect_rec * muBOut.slice(i).row(j-1).t();
           fixEffRec = fixed_effect_rec * muBOut.slice(i).row(j-1).t();
           //RECEIVER MODEL
           betaSample = UpdateDelta(u, v, K, beta_hat, rndMeanRec, omega_hat_11, omega_hat_21, omega_hat_22, psiAdj, fixEffRec, sigmaCube.slice(j-1));
         }
         muBOut.slice(i).row(j) = UpdateMean(u, v, K, SigmaPsi, psi_hat, betaSample, randomModelRec);
         sigmaCube.slice(j) = UpdateSigma(eta, K, v, lambdaBOut.slice(i).row(j-1).t(), betaSample);
         lambdaBOut.slice(i).row(j) = UpdateLambda(v, eta, sigmaCube.slice(j-1), xi);
       }
       if(randomModelSnd){
         if(p == q){
           gammaSample = UpdateDeltaInd(p, K, gamma_hat, zeta_hat_22, muGOut.slice(i).row(j-1).t(), zetaCube.slice(j-1));
         } else {
           //separating the fixed-effect from the random-effect mean
           rndMeanSnd = random_effect_snd * muGOut.slice(i).row(j-1).t();
           fixEffSnd = fixed_effect_snd * muGOut.slice(i).row(j-1).t();
           //SENDER MODEL
           gammaSample = UpdateDelta(q, p, K, gamma_hat, rndMeanSnd, zeta_hat_11, zeta_hat_21, zeta_hat_22, phiAdj, fixEffSnd, zetaCube.slice(j-1));
         }
         muGOut.slice(i).row(j) = UpdateMean(q, p, K, SigmaPhi, phi_hat, gammaSample, randomModelSnd);
         zetaCube.slice(j) = UpdateSigma(eta, K, p, lambdaGOut.slice(i).row(j-1).t(), gammaSample);
         lambdaGOut.slice(i).row(j) = UpdateLambda(p, eta, zetaCube.slice(j-1), xi);
       }
     }

     if(!randomModelRec){//RECEIVER MODEL
       muBOut.slice(i).row(j) = UpdateMean(u, v, K, SigmaPsi, psi_hat, beta, randomModelRec); //sampling fixed effects
     }

     if(!randomModelSnd){//SENDER MODEL
       muGOut.slice(i).row(j) = UpdateMean(q, p, K, SigmaPhi, phi_hat, gamma, randomModelSnd); //sampling fixed effects
     }

     if(randomModelRec){
       betaSample = betaSample.t(); //RECEIVER
       for(int l = 0; l < K; ++l){
         betaCube.slice(l).row(j) = betaSample.row(l); //RECEIVER
       }
     }

     if(randomModelSnd){
      gammaSample = gammaSample.t(); //SENDER
      for(int l = 0; l < K; ++l){
         gammaCube.slice(l).row(j) = gammaSample.row(l); //SENDER
      }
     }

   }

   if(randomModelRec){
     //RECEIVER MODEL
     betaOut[i] = burnInRows(Nthin, NburnIn, betaCube); //already burning in
     sigmaOut[i] = burnInSlices(Nthin, NburnIn, sigmaCube); //already burning in
   }

   if(randomModelSnd){
     //SENDER MODEL
     gammaOut[i] = burnInRows(Nthin, NburnIn, gammaCube); //already burning in
     zetaOut[i] = burnInSlices(Nthin, NburnIn, zetaCube); //already burning in
   }

 }

 //////////////////////////////////////////////////////////////////////////////////

 if(randomModelRec){//RECEIVER MODEL
   lambdaBOut = burnInRows(Nthin, NburnIn, lambdaBOut); //burning in the scale mixture of the cov matrix
 }

 if(randomModelSnd){//SENDER MODEL
   lambdaGOut = burnInRows(Nthin, NburnIn, lambdaGOut); //burning in the scale mixture of the cov matrix
 }
 //RECEIVER
 muBOut = burnInRows(Nthin, NburnIn, muBOut); //burning in the fixed effects

 //SENDER
 muGOut = burnInRows(Nthin, NburnIn, muGOut); //burning in the fixed effects

 ///OUTPUT OF THE FUNCTIONS

 Rcpp::List S1; //sender output
 Rcpp::List S2; //receiver output

 //Sender model contains fixed effects
 if(randomModelSnd){ //model contains
    S1 = Rcpp::List::create(Rcpp::Named("gamma") = gammaOut,
                            Rcpp::Named("mu") = muGOut,
                            Rcpp::Named("sigma") = zetaOut,
                            Rcpp::Named("alpha") = lambdaGOut);
  } else {
     S1 = Rcpp::List::create(Rcpp::Named("mu") = muGOut);
  }
 //Receiver model contains fixed effects
 if(randomModelRec){ //model contains
    S2 = Rcpp::List::create(Rcpp::Named("beta") = betaOut, //random effects
                            Rcpp::Named("mu") = muBOut, //random-effect means
                            Rcpp::Named("sigma") = sigmaOut, //random-effect covariance matrix
                            Rcpp::Named("alpha") = lambdaBOut); //mixture variable of cov-matrix prior
  } else {
    S2 = Rcpp::List::create(Rcpp::Named("mu") = muBOut);
  }
 return(Rcpp::List::create(Rcpp::Named("sender_model") = S1,
                           Rcpp::Named("receiver_model") = S2));
}
