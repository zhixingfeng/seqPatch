#ifndef EBMIXTURE_POOLING_H
#define EBMIXTURE_POOLING_H

#include <stdio.h>
#include <math.h>
#include "SpecialFuns.h"
#include "stl.h"
#include "math_utils.h"

#ifndef ERR
#define ERR 1e-10
#endif


class EBmixture_pooling
{
	public:
		EBmixture_pooling(){mu_0 = 0; sigma_0 = 1; iter = 0; max_iter = 500;}
		void setParameters(double _mu_0, double _sigma_0, vector<double>  _f1_x, vector<double> _f1_y)
		{
			mu_0 = _mu_0; 
			sigma_0 = _sigma_0;
			f1_x = _f1_x;
			f1_y = _f1_y;
		}
		void setMaxIter(double _max_iter){max_iter = _max_iter;}
		void clear()
		{	
			rho_0.clear(); rho_1.clear();	
			gamma_0.clear(); gamma_1.clear();
			N_0.clear(); N_1.clear();
			prop.clear();
		}
		void setData(vector<double> z_score){z = z_score;}
		vector<int> bin_search(double query, double *temp, int temp_len);		
		bool run();

		double get_mu_0(){return mu_0;}
		double get_sigma_0(){return sigma_0;}
		vector<double> get_f1_x(){return f1_x;}
		vector<double> get_f1_y(){return f1_y;}

		vector<double> get_z(){return z;}	
	
		double f0(double x);
		double f0_log(double x);
		double f1(double x);
		double f1_log(double x);	

		vector<double> get_rho_0(){return rho_0;}
		vector<double> get_rho_1(){return rho_1;}
		vector<double> get_gamma_0(){return gamma_0;}
		vector<double> get_gamma_1(){return gamma_1;}
	
		double get_N_0(){return N_0.back();}
		double get_N_1(){return N_1.back();}	
		double get_prop(){return prop.back();}
	
		vector<double> get_N_0_track(){return N_0;}
		vector<double> get_N_1_track(){return N_1;}
		vector<double> get_prop_track(){return prop;}

	protected:
		// parameters
		double mu_0;
		double sigma_0;
		vector<double> f1_x;
		vector<double> f1_y;
		int max_iter;
			
		// data
		vector<double> z;

		// results
		vector<double> rho_0;
		vector<double> rho_1;
		vector<double> gamma_0;
		vector<double> gamma_1;
		vector<double> N_0;
		vector<double> N_1;
		vector<double> prop;	

		int iter;					
};




#endif


