#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include "SpecialFuns.h"

using namespace std;
using namespace Rcpp;
string SP_reverseSeq(const string &seq);

inline int getMoleculeMeanIPD (double *ipd, double *idx, int len_ipd, int len_idx, vector<double> &ipd_avg, vector<double> &ipd_n)
{
	/*-----------Get mean IPD of each molecule--------------*/
        if (len_ipd != len_idx){ Rprintf("length of R_IPD and R_idx should be the same."); return 0;}
       	int len = len_ipd;
	pair<map<int,double>::iterator,bool> ret;
        map<int, double> ipd_sum_map;
        map<int, double> ipd_n_map;
        for (int i=0;i<len;i++){
                ret = ipd_sum_map.insert(pair<int, double>(int(idx[i]),ipd[i]));
                if (ret.second==true){
                        ipd_n_map.insert(pair<int, double>(int(idx[i]),1));
                }else{
                        ret.first->second += ipd[i];
                        ipd_n_map[int(idx[i])] += 1;
                }
        }

       	ipd_avg.clear();
        ipd_n.clear();
        map<int, double>::iterator it_avg = ipd_sum_map.begin();
        map<int, double>::iterator it_n = ipd_n_map.begin();
        while(it_avg!=ipd_sum_map.end()){
                ipd_avg.push_back(it_avg->second / it_n->second);
                ipd_n.push_back(it_n->second);
                it_avg++;
                it_n++;
        }
	return 0;	
}
vector<double> BayesianMixtureModel_core (vector<double> &ipd_avg, vector<double> &ipd_n, 
		double theta_0, double kappa_0, double upsilon_0, double tau2_0, 
		double theta_1, double kappa_1, double upsilon_1, double tau2_1)
{
	/*------------Fit Bayesian Gaussian Mixture Model----------------*/

        // set initial values
        vector<double> theta_0_t; vector<double> kappa_0_t;
        vector<double> upsilon_0_t; vector<double> tau2_0_t;
        vector<double> theta_1_t; vector<double> kappa_1_t;
        vector<double> upsilon_1_t; vector<double> tau2_1_t;
        vector<double> N_0_t; vector <double> N_gamma_0_t;
        vector<double> N_1_t; vector <double> N_gamma_1_t;

        theta_0_t.push_back(theta_0); kappa_0_t.push_back(kappa_0);
        upsilon_0_t.push_back(upsilon_0); tau2_0_t.push_back(tau2_0);

        theta_1_t.push_back(theta_1); kappa_1_t.push_back(kappa_1);
        upsilon_1_t.push_back(upsilon_1); tau2_1_t.push_back(tau2_1);

        N_0_t.push_back(0); N_gamma_0_t.push_back(0);
        N_1_t.push_back(0); N_gamma_1_t.push_back(0);
	
}
RcppExport SEXP BayesianMixtureModel (SEXP R_IPD, SEXP R_idx, SEXP R_theta_0, SEXP R_kappa_0, SEXP R_upsilon_0, SEXP R_tau2_0, 
				SEXP R_theta_1, SEXP R_kappa_1, SEXP R_upsilon_1, SEXP R_tau2_1, SEXP R_max_iter)
{
	/*-----------Get mean IPD of each molecule--------------*/
	vector<double> ipd_avg;
	vector<double> ipd_n;
	getMoleculeMeanIPD (REAL(R_IPD), REAL(R_idx), Rf_length(R_IPD), Rf_length(R_idx), ipd_avg, ipd_n);

	/*------------Fit Bayesian Gaussian Mixture Model----------------*/	
	// load hyper parameters 
	double theta_0 = REAL(R_theta_0)[0];	
	double kappa_0 = REAL(R_kappa_0)[0];
	double upsilon_0 = REAL(R_upsilon_0)[0];
	double tau2_0 = REAL(R_tau2_0)[0];

	double theta_1 = REAL(R_theta_1)[0];
        double kappa_1 = REAL(R_kappa_1)[0];
        double upsilon_1 = REAL(R_upsilon_1)[0];
        double tau2_1 = REAL(R_tau2_1)[0];
			

	// variational inference 
	BayesianMixtureModel_core(ipd_avg, ipd_n, theta_0, kappa_0, upsilon_0, tau2_0,
				theta_1, kappa_1, upsilon_1, tau2_1);	

		
	//return Rcpp::List::create(Rcpp::Named("ipd_avg")=Rcpp::wrap(ipd_avg), Rcpp::Named("ipd_n")=Rcpp::wrap(ipd_n));
	
        return R_NilValue;
}

RcppExport SEXP DetectModProp_NC(SEXP R_IPD, SEXP R_idx)
{
	for (int i=0;i<1000;i++){
		vector<double> ipd_avg;
        	vector<double> ipd_n;
        	getMoleculeMeanIPD (REAL(R_IPD), REAL(R_idx), Rf_length(R_IPD), Rf_length(R_idx), ipd_avg, ipd_n);
	}

        return R_NilValue;
}



