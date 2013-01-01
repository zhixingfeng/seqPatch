#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include "SpecialFuns.h"
#include "BayesianMixtureModel.h"


using namespace std;
using namespace Rcpp;
string SP_reverseSeq(const string &seq);

RcppExport SEXP R_API_BayesianMixtureModel (SEXP R_IPD, SEXP R_idx, SEXP R_theta_0, SEXP R_kappa_0, SEXP R_upsilon_0, SEXP R_tau2_0, 
				SEXP R_theta_1, SEXP R_kappa_1, SEXP R_upsilon_1, SEXP R_tau2_1, SEXP R_max_iter)
{
	BayesianMixtureModel_NC BayesianMixtureModelObj;
	BayesianMixtureModelObj.setHyperParametersNull(REAL(R_theta_0)[0], REAL(R_kappa_0)[0], REAL(R_upsilon_0)[0], REAL(R_tau2_0)[0]);
	BayesianMixtureModelObj.getMoleculeMeanIPD(REAL(R_IPD), REAL(R_idx), Rf_length(R_IPD), Rf_length(R_idx));
	BayesianMixtureModelObj.run(INTEGER(R_max_iter)[0]);		
	

	return Rcpp::List::create(Rcpp::Named("theta_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_theta_0_t_track()),
				Rcpp::Named("kappa_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_kappa_0_t_track()),
	 			Rcpp::Named("upsilon_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_upsilon_0_t_track()),
				Rcpp::Named("tau2_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_tau2_0_t_track()),
				Rcpp::Named("N_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_0_t_track()),
				Rcpp::Named("N_gamma_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_gamma_0_t_track()),

				Rcpp::Named("theta_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_theta_1_t_track()),
                                Rcpp::Named("kappa_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_kappa_1_t_track()),
                                Rcpp::Named("upsilon_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_upsilon_1_t_track()),
                                Rcpp::Named("tau2_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_tau2_1_t_track()),
                                Rcpp::Named("N_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_1_t_track()),
                                Rcpp::Named("N_gamma_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_gamma_1_t_track()),
		
				Rcpp::Named("gamma_0")=Rcpp::wrap(BayesianMixtureModelObj.get_gamma_0()),
				Rcpp::Named("gamma_1")=Rcpp::wrap(BayesianMixtureModelObj.get_gamma_1()));		

/*
	return Rcpp::List::create(Rcpp::Named("ipd_avg")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_avg()),
		Rcpp::Named("ipd_n")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_n()), 
		Rcpp::Named("ipd_var")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_var()));	
  
*/	
      	return R_NilValue;
}

RcppExport SEXP DetectModProp_NC(SEXP R_IPD, SEXP R_idx)
{

        return R_NilValue;
}



