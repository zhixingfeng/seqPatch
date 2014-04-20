#ifndef EBMIXTURE_EM_NAIVE_H
#define EBMIXTURE_EM_NAIVE_H

#include <stdio.h>
#include <math.h>
#include "SpecialFuns.h"
#include "stl.h"
#include "math_utils.h"
#include "sampling.h"

#ifndef ERR
#define ERR 1e-10
#endif




class EBmixture_EM_naive
{
	public:
		EBmixture_EM_naive(){mu_0 = 0; sigma_0 = 1; n_mol = -1; }
		void setParameters(double _mu_0, double _sigma_0,int _max_iter)
		{
			mu_0 = _mu_0; 
			sigma_0 = _sigma_0;
			sigma_1 = _sigma_0;
			max_iter = _max_iter;
				
			mu_1 = _mu_0;
		}
		void clear()
		{
			gamma_0.clear(); gamma_1.clear();
			N_0.clear(); N_1.clear();
			prop.clear();
		}
		
		bool getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx);
		bool run();
		bool run_reinit();
		double get_mu_0(){return mu_0;}
		double get_sigma_0(){return sigma_0;}

	
		double f0(double x_avg, double x_var, double x_n);
		double f0_log(double x_avg, double x_var,  double x_n);
		double f1(double x_avg, double x_var, double x_n);
                double f1_log(double x_avg, double x_var,  double x_n);

		// get data 
		vector<double> get_ipd_avg(){ return ipd_avg; }
                vector<double> get_ipd_n(){ return ipd_n; }
                vector<double> get_ipd_var() { return ipd_var;}
	
		// get result
		int get_n_mol(){ return n_mol;}
		double get_avg_n () {return mean(ipd_n);}
		
		vector<double> get_gamma_0(){return gamma_0;}
		vector<double> get_gamma_1(){return gamma_1;}

		double get_N_0(){return N_0.back();}
		double get_N_1(){return N_1.back();}	
		double get_prop(){return prop.back();}
		
		double get_mu_1(){return mu_1;}

		vector<double> get_N_0_track(){return N_0;}
		vector<double> get_N_1_track(){return N_1;}
		vector<double> get_prop_track(){return prop;}	
		
		double get_n_iter(){return (double) prop.size()-1;}
		
		double get_log_likely_0(){return log_likely_0;}
		double get_log_likely_1(){return log_likely_1;}

		double get_BIC_0(){return BIC_0;}
		double get_BIC_1(){return BIC_1;}
	
		double get_AIC_0(){return AIC_0;}
		double get_AIC_1(){return AIC_1;}

	protected:
		// parameters
		double mu_0;
		double sigma_0;
		double sigma_1;
		int max_iter;		
			
		// data
		vector<double> ipd_avg;
		vector<double> ipd_n;
		vector<double> ipd_var;

		// results
		int n_mol;	 
		
		vector<double> prop;
		vector<double> gamma_0;
		vector<double> gamma_1;
		vector<double> N_0;
		vector<double> N_1;

		double mu_1;

		double log_likely_0;
		double log_likely_1;		
		
		double BIC_0;
		double BIC_1;
	
		double AIC_0;
		double AIC_1;
		int iter;

		vector<double> log_likely_1_all;
		vector<double> BIC_1_all;
		vector<double> AIC_1_all;
			
		vector<double> mu_1_all;
		vector<double> prop_all;
		vector<double> N_0_all;
		vector<double> N_1_all;
};




#endif


