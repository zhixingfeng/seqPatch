#include "BayesianMixtureModel.h"

bool BayesianMixtureModel::getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx)
{
	if (len_ipd != len_idx){ printf("length of R_IPD and R_idx should be the same."); return false;}
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
	
	if (ipd_avg.size()!=ipd_n.size()){printf("run error: size of ipd_avg and ipd_n should be the same."); return false;}
        int T = (int) ipd_avg.size();
	gamma_0.insert(gamma_0.end(), T, -1);
	gamma_1.insert(gamma_1.end(), T, -1);	
	return true;
}

bool BayesianMixtureModel_NC::run(int max_iter)
{
	if (lock==true){printf("run() is locked. rerun setHyperParametersNull()."); return false;}
	/*----------variational inference----------*/		
	for (int t=1;t<=max_iter;t++){
		// estimate q_delta (i.e. estimate gamma_0 and gamma_1)
		double E_log_1_p = digamma(N_0_t[t-1] + 1) - digamma(N_0_t[t-1] + N_1_t[t-1] + 2);
		double E_log_p = digamma(N_1_t[t-1] + 1) - digamma(N_0_t[t-1] + N_1_t[t-1] + 2);
		double E_log_sigma2_0 = log(upsilon_0_t[t-1]*tau2_0_t[t-1]/2) - digamma(upsilon_0_t[t-1]/2);
		double E_log_sigma2_1 = log(upsilon_1_t[t-1]*tau2_1_t[t-1]/2) - digamma(upsilon_1_t[t-1]/2);	
		
		double E_var_norm_0; double E_var_norm_1;
		double rho_0; double rho_1;

		double cur_N_0 = 0; double cur_N_1 = 0;
		double cur_N_gamma_0 = 0; double cur_N_gamma_1 = 0;
	
		double N_y_0_tilde = 0; double N_y_1_tilde = 0;
		double N_S2_0_tilde = 0; double N_S2_1_tilde = 0;
		
		for (int i=0;i<T;i++){
			E_var_norm_0 = 1/kappa_0_t[t-1] + pow(theta_0_t[t-1] - ipd_avg[i], 2)/tau2_0_t[t-1];
			E_var_norm_1 = 1/kappa_1_t[t-1] + pow(theta_1_t[t-1] - ipd_avg[i], 2)/tau2_1_t[t-1];
			rho_0 = exp(E_log_1_p - E_var_norm_0 - E_log_sigma2_0);
			rho_1 = exp(E_log_p - E_var_norm_1 - E_log_sigma2_1);
			
			// estimate gamma_0 and gamma_1
			gamma_0[i] = rho_0 / (rho_0 + rho_1);
			gamma_1[i] = 1 - gamma_0[i];
			
			// get N_0 and N_1
			cur_N_0 += gamma_0[i]; 
			cur_N_1 += gamma_1[i];

			// get N_gamma_0 and N_gamma_1
			double cur_N_gamma_0_i = gamma_0[i]*ipd_n[i];
			double cur_N_gamma_1_i = gamma_1[i]*ipd_n[i];
			cur_N_gamma_0 += cur_N_gamma_0_i;
			cur_N_gamma_1 += cur_N_gamma_1_i;			
		
			// get y_0_tilde and y_1_tilde
			double N_y_0_tilde_i = cur_N_gamma_0_i*ipd_avg[i];
			double N_y_1_tilde_i = cur_N_gamma_1_i*ipd_avg[i];
			N_y_0_tilde += N_y_0_tilde_i;
			N_y_1_tilde += N_y_1_tilde_i;
		
			// get S2_0_tilde and S2_1_tilde
			N_S2_0_tilde += N_y_0_tilde_i*ipd_avg[i];
			N_S2_1_tilde += N_y_1_tilde_i*ipd_avg[i];
		}
		if (cur_N_gamma_0 <= ERR)
			N_S2_0_tilde = ERR;
		else 
			N_S2_0_tilde = N_S2_0_tilde - N_y_0_tilde*N_y_0_tilde/N_y_0_tilde;

		if (cur_N_gamma_1 <= ERR)
                        N_S2_1_tilde = ERR;
                else
                        N_S2_1_tilde = N_S2_1_tilde - N_y_1_tilde*N_y_1_tilde/N_y_1_tilde;

		// estimate q_p (i.e. estimate N_0 and N_1)
		N_0_t.push_back(cur_N_0);
		N_1_t.push_back(cur_N_1);
		
		// estimate q_mu_sigma2 (i.e. estimate theta, kappa, upsilon, and tau2)
		N_gamma_0_t.push_back(cur_N_gamma_0);
		N_gamma_1_t.push_back(cur_N_gamma_1);

		double cur_kappa_0_t = kappa_0 + cur_N_gamma_0;
		double cur_theta_0_t = (N_y_0_tilde + kappa_0*theta_0)/cur_kappa_0_t;
		double cur_upsilon_0_t = upsilon_0 + cur_N_0;
		double cur_tau2_0_t;
		if (cur_N_gamma_0 <= ERR)
			cur_tau2_0_t = upsilon_0*tau2_0/cur_upsilon_0_t;
		else 
			cur_tau2_0_t = (N_S2_0_tilde + upsilon_0*tau2_0 + 
				kappa_0*cur_N_gamma_0*(N_y_0_tilde/cur_N_gamma_0 - theta_0)*(N_y_0_tilde/cur_N_gamma_0 - theta_0)/cur_kappa_0_t)/cur_upsilon_0_t;
	
		double cur_kappa_1_t = kappa_1 + cur_N_gamma_1;
                double cur_theta_1_t = (N_y_1_tilde + kappa_1*theta_1)/cur_kappa_1_t;
                double cur_upsilon_1_t = upsilon_1 + cur_N_1;
                double cur_tau2_1_t;
                if (cur_N_gamma_1 <= ERR)
                        cur_tau2_1_t = upsilon_1*tau2_1/cur_upsilon_1_t;
                else
                        cur_tau2_1_t = (N_S2_1_tilde + upsilon_1*tau2_1 +
                                kappa_1*cur_N_gamma_1*(N_y_1_tilde/cur_N_gamma_1 - theta_1)*(N_y_1_tilde/cur_N_gamma_1 - theta_1)/cur_kappa_1_t)/cur_upsilon_1_t;
		
		
	}

	// lock 
	lock = true;	
	return true;	
}


