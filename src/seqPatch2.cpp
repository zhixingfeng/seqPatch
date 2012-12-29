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
	BayesianMixtureModelObj.run();		
	

		
	return Rcpp::List::create(Rcpp::Named("ipd_avg")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_avg()),
		 Rcpp::Named("ipd_n")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_n() ));	
        return R_NilValue;
}

RcppExport SEXP DetectModProp_NC(SEXP R_IPD, SEXP R_idx)
{

        return R_NilValue;
}



