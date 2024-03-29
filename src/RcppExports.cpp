// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// mat2DF
Rcpp::DataFrame mat2DF(const arma::mat& x);
RcppExport SEXP _remx_mat2DF(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(mat2DF(x));
    return rcpp_result_gen;
END_RCPP
}
// rehPoisson
Rcpp::DataFrame rehPoisson(Rcpp::List stats);
RcppExport SEXP _remx_rehPoisson(SEXP statsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type stats(statsSEXP);
    rcpp_result_gen = Rcpp::wrap(rehPoisson(stats));
    return rcpp_result_gen;
END_RCPP
}
// burnInRows
arma::cube burnInRows(double Nthin, double Nburn, const arma::cube& samples);
RcppExport SEXP _remx_burnInRows(SEXP NthinSEXP, SEXP NburnSEXP, SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Nthin(NthinSEXP);
    Rcpp::traits::input_parameter< double >::type Nburn(NburnSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(burnInRows(Nthin, Nburn, samples));
    return rcpp_result_gen;
END_RCPP
}
// burnInSlices
arma::cube burnInSlices(double Nthin, double Nburn, const arma::cube& samples);
RcppExport SEXP _remx_burnInSlices(SEXP NthinSEXP, SEXP NburnSEXP, SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type Nthin(NthinSEXP);
    Rcpp::traits::input_parameter< double >::type Nburn(NburnSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(burnInSlices(Nthin, Nburn, samples));
    return rcpp_result_gen;
END_RCPP
}
// rndMeans
arma::vec rndMeans(const arma::mat& x, const arma::vec& psi);
RcppExport SEXP _remx_rndMeans(SEXP xSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(rndMeans(x, psi));
    return rcpp_result_gen;
END_RCPP
}
// UpdateBeta
arma::mat UpdateBeta(int p, int K, const arma::mat& beta_hat, const arma::cube& omega_hat, const arma::vec& Mu, const arma::mat& Sigma);
RcppExport SEXP _remx_UpdateBeta(SEXP pSEXP, SEXP KSEXP, SEXP beta_hatSEXP, SEXP omega_hatSEXP, SEXP MuSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat(omega_hatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateBeta(p, K, beta_hat, omega_hat, Mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// UpdateOmega
arma::vec UpdateOmega(int p, int K, const arma::mat& beta_hat, const arma::cube& omega_hat, const arma::vec& Mu, const arma::mat& Sigma, const arma::vec& nu, const arma::mat& beta);
RcppExport SEXP _remx_UpdateOmega(SEXP pSEXP, SEXP KSEXP, SEXP beta_hatSEXP, SEXP omega_hatSEXP, SEXP MuSEXP, SEXP SigmaSEXP, SEXP nuSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat(omega_hatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateOmega(p, K, beta_hat, omega_hat, Mu, Sigma, nu, beta));
    return rcpp_result_gen;
END_RCPP
}
// UpdateMu
arma::mat UpdateMu(int p, int K, const arma::mat& TauMu, const arma::mat& Sigma, const arma::mat& beta);
RcppExport SEXP _remx_UpdateMu(SEXP pSEXP, SEXP KSEXP, SEXP TauMuSEXP, SEXP SigmaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type TauMu(TauMuSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateMu(p, K, TauMu, Sigma, beta));
    return rcpp_result_gen;
END_RCPP
}
// UpdateSigma
arma::mat UpdateSigma(int eta, int K, int p, const arma::vec& Lambda, const arma::mat& beta);
RcppExport SEXP _remx_UpdateSigma(SEXP etaSEXP, SEXP KSEXP, SEXP pSEXP, SEXP LambdaSEXP, SEXP betaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta(betaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateSigma(eta, K, p, Lambda, beta));
    return rcpp_result_gen;
END_RCPP
}
// UpdateLambda
arma::rowvec UpdateLambda(int p, int eta, const arma::mat& Sigma, double xi);
RcppExport SEXP _remx_UpdateLambda(SEXP pSEXP, SEXP etaSEXP, SEXP SigmaSEXP, SEXP xiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< double >::type xi(xiSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateLambda(p, eta, Sigma, xi));
    return rcpp_result_gen;
END_RCPP
}
// UpdatePsi
arma::rowvec UpdatePsi(int q, int K, const arma::vec& TauPsi, const arma::mat& SigmaPsi, const arma::mat& psi);
RcppExport SEXP _remx_UpdatePsi(SEXP qSEXP, SEXP KSEXP, SEXP TauPsiSEXP, SEXP SigmaPsiSEXP, SEXP psiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type TauPsi(TauPsiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type SigmaPsi(SigmaPsiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdatePsi(q, K, TauPsi, SigmaPsi, psi));
    return rcpp_result_gen;
END_RCPP
}
// UpdateMean
arma::rowvec UpdateMean(int q, int p, int K, const arma::cube& SigmaMean, const arma::mat& psi, const arma::mat& delta, bool randomModel);
RcppExport SEXP _remx_UpdateMean(SEXP qSEXP, SEXP pSEXP, SEXP KSEXP, SEXP SigmaMeanSEXP, SEXP psiSEXP, SEXP deltaSEXP, SEXP randomModelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type SigmaMean(SigmaMeanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type randomModel(randomModelSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateMean(q, p, K, SigmaMean, psi, delta, randomModel));
    return rcpp_result_gen;
END_RCPP
}
// UpdateDelta
arma::mat UpdateDelta(int q, int p, int K, const arma::mat& beta_hat, const arma::vec& Mu, const arma::cube& omega_hat_11, const arma::cube& omega_hat_21, const arma::cube& omega_hat_22, const arma::mat& psi_hat, const arma::vec& psi, const arma::mat& Sigma);
RcppExport SEXP _remx_UpdateDelta(SEXP qSEXP, SEXP pSEXP, SEXP KSEXP, SEXP beta_hatSEXP, SEXP MuSEXP, SEXP omega_hat_11SEXP, SEXP omega_hat_21SEXP, SEXP omega_hat_22SEXP, SEXP psi_hatSEXP, SEXP psiSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_11(omega_hat_11SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_21(omega_hat_21SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_22(omega_hat_22SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi_hat(psi_hatSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateDelta(q, p, K, beta_hat, Mu, omega_hat_11, omega_hat_21, omega_hat_22, psi_hat, psi, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// UpdateDeltaInd
arma::mat UpdateDeltaInd(int p, int K, const arma::mat& beta_hat, const arma::cube& omega_hat_22, const arma::vec& Mu, const arma::mat Sigma);
RcppExport SEXP _remx_UpdateDeltaInd(SEXP pSEXP, SEXP KSEXP, SEXP beta_hatSEXP, SEXP omega_hat_22SEXP, SEXP MuSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_22(omega_hat_22SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(UpdateDeltaInd(p, K, beta_hat, omega_hat_22, Mu, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// sampler
Rcpp::List sampler(const int Niter, const int Nchain, const int NburnIn, const int Nthin, const int p, const int q, const int K, const arma::cube delta, const arma::mat& beta_hat, const arma::cube& omega_hat, const arma::cube& omega_hat_11, const arma::cube& omega_hat_21, const arma::cube& omega_hat_22, const arma::mat& Mu, const arma::cube Sigma, double eta, const arma::mat& Lambda, const double& xi, const arma::cube& SigmaMean, const arma::mat& psi, const arma::mat& psi_hat, const arma::mat& random_effect, const arma::mat& fixed_effect, bool randomModel);
RcppExport SEXP _remx_sampler(SEXP NiterSEXP, SEXP NchainSEXP, SEXP NburnInSEXP, SEXP NthinSEXP, SEXP pSEXP, SEXP qSEXP, SEXP KSEXP, SEXP deltaSEXP, SEXP beta_hatSEXP, SEXP omega_hatSEXP, SEXP omega_hat_11SEXP, SEXP omega_hat_21SEXP, SEXP omega_hat_22SEXP, SEXP MuSEXP, SEXP SigmaSEXP, SEXP etaSEXP, SEXP LambdaSEXP, SEXP xiSEXP, SEXP SigmaMeanSEXP, SEXP psiSEXP, SEXP psi_hatSEXP, SEXP random_effectSEXP, SEXP fixed_effectSEXP, SEXP randomModelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const int >::type Nchain(NchainSEXP);
    Rcpp::traits::input_parameter< const int >::type NburnIn(NburnInSEXP);
    Rcpp::traits::input_parameter< const int >::type Nthin(NthinSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type delta(deltaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat(omega_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_11(omega_hat_11SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_21(omega_hat_21SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_22(omega_hat_22SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Mu(MuSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type Lambda(LambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type SigmaMean(SigmaMeanSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi_hat(psi_hatSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type random_effect(random_effectSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fixed_effect(fixed_effectSEXP);
    Rcpp::traits::input_parameter< bool >::type randomModel(randomModelSEXP);
    rcpp_result_gen = Rcpp::wrap(sampler(Niter, Nchain, NburnIn, Nthin, p, q, K, delta, beta_hat, omega_hat, omega_hat_11, omega_hat_21, omega_hat_22, Mu, Sigma, eta, Lambda, xi, SigmaMean, psi, psi_hat, random_effect, fixed_effect, randomModel));
    return rcpp_result_gen;
END_RCPP
}
// samplerActor
Rcpp::List samplerActor(const int Niter, const int Nchain, const int NburnIn, const int Nthin, const int v, const int u, const int K, const arma::cube beta, const arma::mat& beta_hat, const arma::cube& omega_hat, const arma::cube& omega_hat_11, const arma::cube& omega_hat_21, const arma::cube& omega_hat_22, const arma::mat& MuB, const arma::cube Sigma, double eta, const arma::mat& LambdaB, const double& xi, const arma::cube& SigmaPsi, const arma::mat& psi, const arma::mat& psi_hat, const int p, const int q, const arma::cube gamma, const arma::mat& gamma_hat, const arma::cube& zeta_hat, const arma::cube& zeta_hat_11, const arma::cube& zeta_hat_21, const arma::cube& zeta_hat_22, const arma::mat& MuG, const arma::cube Zeta, const arma::mat& LambdaG, const arma::cube& SigmaPhi, const arma::mat& phi, const arma::mat& phi_hat, const arma::mat& random_effect_rec, const arma::mat& fixed_effect_rec, const arma::mat& random_effect_snd, const arma::mat& fixed_effect_snd, bool randomModelSnd, bool randomModelRec);
RcppExport SEXP _remx_samplerActor(SEXP NiterSEXP, SEXP NchainSEXP, SEXP NburnInSEXP, SEXP NthinSEXP, SEXP vSEXP, SEXP uSEXP, SEXP KSEXP, SEXP betaSEXP, SEXP beta_hatSEXP, SEXP omega_hatSEXP, SEXP omega_hat_11SEXP, SEXP omega_hat_21SEXP, SEXP omega_hat_22SEXP, SEXP MuBSEXP, SEXP SigmaSEXP, SEXP etaSEXP, SEXP LambdaBSEXP, SEXP xiSEXP, SEXP SigmaPsiSEXP, SEXP psiSEXP, SEXP psi_hatSEXP, SEXP pSEXP, SEXP qSEXP, SEXP gammaSEXP, SEXP gamma_hatSEXP, SEXP zeta_hatSEXP, SEXP zeta_hat_11SEXP, SEXP zeta_hat_21SEXP, SEXP zeta_hat_22SEXP, SEXP MuGSEXP, SEXP ZetaSEXP, SEXP LambdaGSEXP, SEXP SigmaPhiSEXP, SEXP phiSEXP, SEXP phi_hatSEXP, SEXP random_effect_recSEXP, SEXP fixed_effect_recSEXP, SEXP random_effect_sndSEXP, SEXP fixed_effect_sndSEXP, SEXP randomModelSndSEXP, SEXP randomModelRecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type Niter(NiterSEXP);
    Rcpp::traits::input_parameter< const int >::type Nchain(NchainSEXP);
    Rcpp::traits::input_parameter< const int >::type NburnIn(NburnInSEXP);
    Rcpp::traits::input_parameter< const int >::type Nthin(NthinSEXP);
    Rcpp::traits::input_parameter< const int >::type v(vSEXP);
    Rcpp::traits::input_parameter< const int >::type u(uSEXP);
    Rcpp::traits::input_parameter< const int >::type K(KSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type beta_hat(beta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat(omega_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_11(omega_hat_11SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_21(omega_hat_21SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type omega_hat_22(omega_hat_22SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type MuB(MuBSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Sigma(SigmaSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LambdaB(LambdaBSEXP);
    Rcpp::traits::input_parameter< const double& >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type SigmaPsi(SigmaPsiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type psi_hat(psi_hatSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type q(qSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type gamma_hat(gamma_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type zeta_hat(zeta_hatSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type zeta_hat_11(zeta_hat_11SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type zeta_hat_21(zeta_hat_21SEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type zeta_hat_22(zeta_hat_22SEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type MuG(MuGSEXP);
    Rcpp::traits::input_parameter< const arma::cube >::type Zeta(ZetaSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type LambdaG(LambdaGSEXP);
    Rcpp::traits::input_parameter< const arma::cube& >::type SigmaPhi(SigmaPhiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type phi_hat(phi_hatSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type random_effect_rec(random_effect_recSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fixed_effect_rec(fixed_effect_recSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type random_effect_snd(random_effect_sndSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type fixed_effect_snd(fixed_effect_sndSEXP);
    Rcpp::traits::input_parameter< bool >::type randomModelSnd(randomModelSndSEXP);
    Rcpp::traits::input_parameter< bool >::type randomModelRec(randomModelRecSEXP);
    rcpp_result_gen = Rcpp::wrap(samplerActor(Niter, Nchain, NburnIn, Nthin, v, u, K, beta, beta_hat, omega_hat, omega_hat_11, omega_hat_21, omega_hat_22, MuB, Sigma, eta, LambdaB, xi, SigmaPsi, psi, psi_hat, p, q, gamma, gamma_hat, zeta_hat, zeta_hat_11, zeta_hat_21, zeta_hat_22, MuG, Zeta, LambdaG, SigmaPhi, phi, phi_hat, random_effect_rec, fixed_effect_rec, random_effect_snd, fixed_effect_snd, randomModelSnd, randomModelRec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_remx_mat2DF", (DL_FUNC) &_remx_mat2DF, 1},
    {"_remx_rehPoisson", (DL_FUNC) &_remx_rehPoisson, 1},
    {"_remx_burnInRows", (DL_FUNC) &_remx_burnInRows, 3},
    {"_remx_burnInSlices", (DL_FUNC) &_remx_burnInSlices, 3},
    {"_remx_rndMeans", (DL_FUNC) &_remx_rndMeans, 2},
    {"_remx_UpdateBeta", (DL_FUNC) &_remx_UpdateBeta, 6},
    {"_remx_UpdateOmega", (DL_FUNC) &_remx_UpdateOmega, 8},
    {"_remx_UpdateMu", (DL_FUNC) &_remx_UpdateMu, 5},
    {"_remx_UpdateSigma", (DL_FUNC) &_remx_UpdateSigma, 5},
    {"_remx_UpdateLambda", (DL_FUNC) &_remx_UpdateLambda, 4},
    {"_remx_UpdatePsi", (DL_FUNC) &_remx_UpdatePsi, 5},
    {"_remx_UpdateMean", (DL_FUNC) &_remx_UpdateMean, 7},
    {"_remx_UpdateDelta", (DL_FUNC) &_remx_UpdateDelta, 11},
    {"_remx_UpdateDeltaInd", (DL_FUNC) &_remx_UpdateDeltaInd, 6},
    {"_remx_sampler", (DL_FUNC) &_remx_sampler, 24},
    {"_remx_samplerActor", (DL_FUNC) &_remx_samplerActor, 41},
    {NULL, NULL, 0}
};

RcppExport void R_init_remx(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
