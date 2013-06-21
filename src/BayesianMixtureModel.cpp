#include "BayesianMixtureModel.h"

bool BayesianMixtureModel::getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx)
{
	if (len_ipd != len_idx){ printf("length of R_IPD and R_idx should be the same.\n"); return false;}
	// get average IPD of each molecule 
	int len = len_ipd;
        pair<map<int,double>::iterator,bool> ret;
        map<int, double> ipd_sum_map;
        map<int, double> ipd_n_map;
	map<int, double> ipd_var_map;        
	for (int i=0;i<len;i++){
                ret = ipd_sum_map.insert(pair<int, double>(int(idx[i]),ipd[i]));
                if (ret.second==true){
                        ipd_n_map.insert(pair<int, double>(int(idx[i]),1));
			ipd_var_map.insert(pair<int, double>(int(idx[i]),ipd[i]*ipd[i]));
                }else{
                        ret.first->second += ipd[i];
                        ipd_n_map[int(idx[i])] += 1;
			ipd_var_map[int(idx[i])] += ipd[i]*ipd[i];
                }
        }

        ipd_avg.clear();
        ipd_n.clear();
	ipd_var.clear();
        map<int, double>::iterator it_avg = ipd_sum_map.begin();
        map<int, double>::iterator it_n = ipd_n_map.begin();
        map<int, double>::iterator it_var = ipd_var_map.begin();

        while(it_avg!=ipd_sum_map.end()){
                double cur_n = it_n->second;
		double cur_avg = it_avg->second / cur_n;
		ipd_avg.push_back(cur_avg);
                ipd_n.push_back(cur_n);
                ipd_var.push_back(it_var->second / cur_n - cur_avg*cur_avg);
		it_avg++;
                it_n++;
        	it_var++;
	}
	
	if ( !(ipd_avg.size()==ipd_n.size()&&ipd_n.size()==ipd_var.size()) ){
		printf("run error: size of ipd_avg, ipd_n and ipd_var should be the same.\n"); 
		return false;
	}
	
        T = (int) ipd_avg.size();
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
		double log_rho_0; double log_rho_1;

		double cur_N_0 = 0; double cur_N_1 = 0;
		double cur_N_gamma_0 = 0; double cur_N_gamma_1 = 0;
	
		double N_y_0_tilde = 0; double N_y_1_tilde = 0;
		double N_S2_0_tilde = 0; double N_S2_1_tilde = 0;
		double N_S2_0_bar = 0; double N_S2_1_bar = 0;	
	
	/*	printf("ipd_avg: ");
		for (int i=0;i<(int)ipd_avg.size();i++) printf("%lf ", ipd_avg[i]);
		printf("\n");
		
		printf("ipd_var: ");
                for (int i=0;i<(int)ipd_avg.size();i++) printf("%lf ", ipd_var[i]);
                printf("\n");
	*/	

		for (int i=0;i<T;i++){
			E_var_norm_0 = 1/kappa_0_t[t-1] + (pow(theta_0_t[t-1] - ipd_avg[i], 2) + ipd_var[i])/tau2_0_t[t-1];
			E_var_norm_1 = 1/kappa_1_t[t-1] + (pow(theta_1_t[t-1] - ipd_avg[i], 2) + ipd_var[i])/tau2_1_t[t-1];
			//printf("E_var_norm_0 : %lf, E_var_norm_1 : %lf\n", E_var_norm_0, E_var_norm_1);
			
			log_rho_0 = E_log_1_p - ipd_n[i]*E_var_norm_0 / 2 - ipd_n[i]*E_log_sigma2_0 / 2;
			log_rho_1 = E_log_p - ipd_n[i]*E_var_norm_1 / 2 - ipd_n[i]*E_log_sigma2_1 / 2;

			// estimate gamma_0 and gamma_1
			if (log_rho_1 - log_rho_0 >= 40)
				gamma_0[i] = 0;
			else 
				gamma_0[i] = 1 / (1 + exp(log_rho_1 - log_rho_0) );
			
			gamma_1[i] = 1 - gamma_0[i];
			//printf("[%lf; (%lf, %lf); (%lf, %lf); (%lf, %lf) ]", gamma_0[i], E_var_norm_0, E_var_norm_1, E_log_1_p, E_log_p ,
			//	E_log_sigma2_0, E_log_sigma2_1);
			//printf("(%lf, %lf, %lf) ", ipd_avg[i], ipd_var[i], gamma_0[i]);	
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
			
			// get S2_0_bar and S2_1_bar
			N_S2_0_bar += cur_N_gamma_0_i*ipd_var[i];
			N_S2_1_bar += cur_N_gamma_1_i*ipd_var[i];
			//printf("ipd_var : %lf, cur_N_gamma_0_i : %lf, cur_N_gamma_1_i : %lf \n", ipd_var[i], cur_N_gamma_0_i, cur_N_gamma_1_i);

		}
		//printf("\n");
		if (cur_N_gamma_0 <= ERR)
			N_S2_0_tilde = ERR;
		else 
			N_S2_0_tilde = N_S2_0_tilde - N_y_0_tilde*N_y_0_tilde/cur_N_gamma_0;

		if (cur_N_gamma_1 <= ERR)
                        N_S2_1_tilde = ERR;
                else
                        N_S2_1_tilde = N_S2_1_tilde - N_y_1_tilde*N_y_1_tilde/cur_N_gamma_1;

		// estimate q_p (i.e. estimate N_0 and N_1)
		N_0_t.push_back(cur_N_0);
		N_1_t.push_back(cur_N_1);
		
		// estimate q_mu_sigma2 (i.e. estimate theta, kappa, upsilon, and tau2)
		N_gamma_0_t.push_back(cur_N_gamma_0);
		N_gamma_1_t.push_back(cur_N_gamma_1);

		double cur_kappa_0_t = kappa_0 + cur_N_gamma_0;
		double cur_theta_0_t = (N_y_0_tilde + kappa_0*theta_0)/cur_kappa_0_t;
		double cur_upsilon_0_t = upsilon_0 + cur_N_gamma_0;
		double cur_tau2_0_t;
		if (cur_N_gamma_0 <= ERR)
			cur_tau2_0_t = tau2_0;
		else 
			cur_tau2_0_t = (N_S2_0_tilde + N_S2_0_bar + upsilon_0*tau2_0 + 
				kappa_0*cur_N_gamma_0*(N_y_0_tilde/cur_N_gamma_0 - theta_0)*(N_y_0_tilde/cur_N_gamma_0 - theta_0)/cur_kappa_0_t)/cur_upsilon_0_t;
	
		double cur_kappa_1_t = kappa_1 + cur_N_gamma_1;
                double cur_theta_1_t = (N_y_1_tilde + kappa_1*theta_1)/cur_kappa_1_t;
                double cur_upsilon_1_t = upsilon_1 + cur_N_gamma_1;
                double cur_tau2_1_t;
                if (cur_N_gamma_1 <= ERR)
                        cur_tau2_1_t = tau2_1;
                else
                        cur_tau2_1_t = (N_S2_1_tilde + N_S2_1_bar + upsilon_1*tau2_1 +
                                kappa_1*cur_N_gamma_1*(N_y_1_tilde/cur_N_gamma_1 - theta_1)*(N_y_1_tilde/cur_N_gamma_1 - theta_1)/cur_kappa_1_t)/cur_upsilon_1_t;
	
		//printf("N_S2_0_bar : %lf, N_S2_1_bar : %lf\n", N_S2_0_bar, N_S2_1_bar);	
		theta_0_t.push_back(cur_theta_0_t);	
		kappa_0_t.push_back(cur_kappa_0_t);
		upsilon_0_t.push_back(cur_upsilon_0_t);	
		tau2_0_t.push_back(cur_tau2_0_t);

		theta_1_t.push_back(cur_theta_1_t);
                kappa_1_t.push_back(cur_kappa_1_t);
                upsilon_1_t.push_back(cur_upsilon_1_t);
                tau2_1_t.push_back(cur_tau2_1_t);
		
		// check convergence
		double eps = 0.001;
		bool is_converge = ( fabs(cur_theta_0_t - theta_0_t[t-1]) < eps && fabs(cur_kappa_0_t - kappa_0_t[t-1])/cur_kappa_0_t < eps && 
					fabs(cur_upsilon_0_t - upsilon_0_t[t-1])/cur_upsilon_0_t < eps && fabs(cur_tau2_0_t - tau2_0_t[t-1])/cur_tau2_0_t < eps && 
					fabs(cur_theta_1_t - theta_1_t[t-1]) < eps && fabs(cur_kappa_1_t - kappa_1_t[t-1])/cur_kappa_1_t < eps &&
                                        fabs(cur_upsilon_1_t - upsilon_1_t[t-1])/cur_upsilon_1_t < eps && fabs(cur_tau2_1_t - tau2_1_t[t-1])/cur_tau2_1_t < eps);
		if (is_converge) break;
	}

	// lock 
	lock = true;	
	return true;	
}


