#ifndef BAYESMODEL
#define BAYESMODEL

#include <stdio.h>
#include <math.h>
#include "SpecialFuns.h"
#include "stl.h"
#ifndef ERR
#define ERR 1e-10
#endif
class BayesianMixtureModel
{
	public:
		BayesianMixtureModel()
		{
			// default hyperparameters
			theta_0 = -2; kappa_0 = 0.0001; upsilon_0 = 0.0001; tau2_0 = 1;
			theta_1 = 2; kappa_1 = 0.0001; upsilon_1 = 0.0001; tau2_1 = 1;
			// lock 
			lock = true;
		}	

		// set hyperparameters
		bool setHyperParameters(double _theta_0, double _kappa_0, double _upsilon_0, double _tau2_0,
				double _theta_1, double _kappa_1, double _upsilon_1, double _tau2_1)
		{
			// load hyperparameters
			/*theta_0 = _theta_0; kappa_0 = _kappa_0 + 10000; upsilon_0 = _upsilon_0 + 10000; tau2_0 = _tau2_0;
			theta_1 = theta_0 ;
			upsilon_1 = _upsilon_0; tau2_1 = _tau2_0;*/
			
			theta_0 = _theta_0; kappa_0 = _kappa_0; upsilon_0 = _upsilon_0; tau2_0 = _tau2_0;
			theta_1 = _theta_1; kappa_1 = _kappa_1; upsilon_1 = _upsilon_1; tau2_1 = _tau2_1;
			// clear distribution parameters			
			theta_0_t.clear(); theta_1_t.clear();
			kappa_0_t.clear(); kappa_1_t.clear();
			
			upsilon_0_t.clear(); upsilon_1_t.clear();
			tau2_0_t.clear(); tau2_1_t.clear();
			
			N_0_t.clear(); N_gamma_0_t.clear();
			N_1_t.clear(); N_gamma_1_t.clear();

			// set initial values
			theta_0_t.push_back(theta_0); kappa_0_t.push_back(kappa_0);
        		upsilon_0_t.push_back(upsilon_0); tau2_0_t.push_back(tau2_0);

        		theta_1_t.push_back(theta_1); kappa_1_t.push_back(kappa_0);
        		upsilon_1_t.push_back(upsilon_0); tau2_1_t.push_back(tau2_0);

        		N_0_t.push_back(0); N_gamma_0_t.push_back(0);
        		N_1_t.push_back(0); N_gamma_1_t.push_back(0);

			// unlock 
			lock = false;
			return true;
		}
		
		// get data
		vector<double> get_ipd_avg(){return ipd_avg;}
		vector<double> get_ipd_n(){return ipd_n;}
		vector<double> get_ipd_var(){return ipd_var;}

		// get tracked distribution parameters
		vector<double> get_theta_0_t_track(){return theta_0_t;}
		vector<double> get_theta_1_t_track(){return theta_1_t;}
		
		vector<double> get_kappa_0_t_track(){return kappa_0_t;}
		vector<double> get_kappa_1_t_track(){return kappa_1_t;}

		vector<double> get_upsilon_0_t_track(){return upsilon_0_t;}
		vector<double> get_upsilon_1_t_track(){return upsilon_1_t;}

		vector<double> get_tau2_0_t_track(){return tau2_0_t;}
		vector<double> get_tau2_1_t_track(){return tau2_1_t;}

		vector<double> get_N_0_t_track(){return N_0_t;}
		vector<double> get_N_1_t_track(){return N_1_t;}

		vector<double> get_N_gamma_0_t_track(){return N_gamma_0_t;}
		vector<double> get_N_gamma_1_t_track(){return N_gamma_1_t;}

		// get estimated distribution parameters
		double get_theta_0_t(){return theta_0_t.back();}	
		double get_theta_1_t(){return theta_1_t.back();}

		double get_kappa_0_t(){return kappa_0_t.back();}
                double get_kappa_1_t(){return kappa_1_t.back();}

                double get_upsilon_0_t(){return upsilon_0_t.back();}
                double get_upsilon_1_t(){return upsilon_1_t.back();}

                double get_tau2_0_t(){return tau2_0_t.back();}
                double get_tau2_1_t(){return tau2_1_t.back();}

                double get_N_0_t(){return N_0_t.back();}
                double get_N_1_t(){return N_1_t.back();}

                double get_N_gamma_0_t(){return N_gamma_0_t.back();}
                double get_N_gamma_1_t(){return N_gamma_1_t.back();}

		vector<double> get_gamma_0(){return gamma_0;}
		vector<double> get_gamma_1(){return gamma_1;}
		
		double get_n_steps(){return  theta_0_t.size() - 1;}
	
		// get mean IPD of each molecule
		bool getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx);
	
	protected:
		// hyperprameters
		double theta_0;
		double kappa_0;
		double upsilon_0;
		double tau2_0;
                double theta_1;
		double kappa_1;
		double upsilon_1;
		double tau2_1;
			
		// track distribution parameters 
		vector<double> theta_0_t; vector<double> theta_1_t;
	        vector<double> kappa_0_t; vector<double> kappa_1_t;
		vector<double> upsilon_0_t; vector<double> upsilon_1_t;
       		vector<double> tau2_0_t; vector<double> tau2_1_t;
        	vector<double> N_0_t; vector <double> N_gamma_0_t;
        	vector<double> N_1_t; vector <double> N_gamma_1_t;

		// probability of each molecule to be modified or not
		vector<double> gamma_0;
		vector<double> gamma_1;

		// data (average IPD of each molecule and number of IPD of each molecule)
		vector<double> ipd_avg;
		vector<double> ipd_n;
		vector<double> ipd_var;
		// number of moelcules
		double T;
	
		// locker 
		bool lock;
};



class BayesianMixtureModel_NC: public BayesianMixtureModel
{
	public:
		bool run(int max_iter=5000);	
	
};

#endif
