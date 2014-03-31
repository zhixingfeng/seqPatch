#ifndef EBMIXTURE_EM_H
#define EBMIXTURE_EM_H

#include <stdio.h>
#include <math.h>
#include "SpecialFuns.h"
#include "stl.h"
#include "math_utils.h"
#include "sampling.h"

#ifndef ERR
#define ERR 1e-10
#endif




class EBmixture_EM
{
	public:
		EBmixture_EM(){mu_0 = 0; sigma_0 = 1; n_mol = -1; }
		void setParameters(double _mu_0, double _sigma_0, vector<double>  _f1_x, vector<double> _f1_y, int _max_iter)
		{
			mu_0 = _mu_0; 
			sigma_0 = _sigma_0;
			sigma_1 = _sigma_0;
			f1_x = _f1_x;
			f1_y = _f1_y;
			max_iter = _max_iter;
			for (int i = 0; i < (int)f1_x.size(); i++) 
				f1_x[i] = f1_x[i] + mu_0;
			s = f1_x[1] - f1_x[0];
		}
		void clear()
		{
			gamma_0.clear(); gamma_1.clear();
			N_0.clear(); N_1.clear();
			prop.clear();
		}
		
		bool getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx);
		vector<int> bin_search(double query, double *temp, int temp_len);		
		bool run(bool is_cal_mu_1 = true, double prop_min=0.2);

		double get_mu_0(){return mu_0;}
		double get_sigma_0(){return sigma_0;}
		vector<double> get_f1_x(){return f1_x;}
		vector<double> get_f1_y(){return f1_y;}

	
		double f0(double x_avg, double x_var, double x_n);
		double f0_log(double x_avg, double x_var,  double x_n);
		double f1(double x_avg, double x_var, double x_n, double _mu_1);
		double f1_log(double x_avg, double x_var, double x_n, double _mu_1);
		double f1_prior(double x);
		double f1_prior_log(double x);	

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
		double get_mu_d(){return mu_d;}

		vector<double> get_N_0_track(){return N_0;}
		vector<double> get_N_1_track(){return N_1;}
		vector<double> get_prop_track(){return prop;}	
		
		double get_n_iter(){return (double) prop.size()-1;}
		

	protected:
		// parameters
		double mu_0;
		double sigma_0;
		double sigma_1;
		vector<double> f1_x;
		vector<double> f1_y;
		int max_iter;		
		double s;	
			
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
		double mu_d;
			
		
		int iter;
};




#endif


