#include "EBmixture_EM_naive.h"

bool EBmixture_EM_naive::getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx)
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

        n_mol = (int) ipd_avg.size();
	return true;
}


double EBmixture_EM_naive::f0(double x_avg, double x_var, double x_n)
{
	return exp( -( x_n* ( (x_avg - mu_0)*(x_avg - mu_0) + x_var) ) / (2.0*sigma_0*sigma_0) ) / pow(sqrt(2.0*my_pi)*sigma_0, x_n);	

}
double EBmixture_EM_naive::f0_log(double x_avg, double x_var, double x_n)
{
	return -( x_n* ( (x_avg - mu_0)*(x_avg - mu_0) + x_var) ) / (2.0*sigma_0*sigma_0) - x_n*log(2.0*my_pi) / 2.0 - x_n*log(sigma_0);
}

double EBmixture_EM_naive::f1(double x_avg, double x_var, double x_n)
{
        return exp( -( x_n* ( (x_avg - mu_1)*(x_avg - mu_1) + x_var) ) / (2.0*sigma_1*sigma_1) ) / pow(sqrt(2.0*my_pi)*sigma_1, x_n);

}
double EBmixture_EM_naive::f1_log(double x_avg, double x_var, double x_n)
{
        return -( x_n* ( (x_avg - mu_1)*(x_avg - mu_1) + x_var) ) / (2.0*sigma_1*sigma_1) - x_n*log(2.0*my_pi) / 2.0 - x_n*log(sigma_1);
}


bool EBmixture_EM_naive::run()
{
	clear();
	
	// Initialize
	N_0.push_back(n_mol / 2.0);
	N_1.push_back(n_mol / 2.0);	
	for (int i=0; i<n_mol; i++){
		gamma_0.push_back(0.5);
		gamma_1.push_back(0.5);
	}
	prop.push_back(0.5);

	
	// iterate
	double eps = 1e-3;	
	iter = 1;
	while (iter <= max_iter){
		// update mu_1
		double numerator = 0;
		double denominator = 0;
		for (int i=0; i<n_mol; i++){
			double cur_coef = gamma_1[i] * ipd_n[i];
			numerator += cur_coef * ipd_avg[i];
			denominator += cur_coef;
		} 
		mu_1 = numerator / denominator;

		// update gamma_0 and gamma_1
		for (int i=0; i<n_mol; i++){
			double cur_f0_log = f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]);	
			double cur_lik_r = exp(f1_log(ipd_avg[i], ipd_var[i], ipd_n[i]) - f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]));
			gamma_0[i] = (1 - prop[iter-1]) / ( prop[iter-1] * cur_lik_r + 1 - prop[iter-1] );
			gamma_1[i] = 1 - gamma_0[i];
		}	

		// update prop
		double cur_N_0 = sum(gamma_0);
		double cur_N_1 = sum(gamma_1);
		N_0.push_back(cur_N_0);
		N_1.push_back(cur_N_1);
		prop.push_back(cur_N_1 / (cur_N_0 + cur_N_1) );	

		if (fabs(prop[iter] - prop[iter-1]) <= eps) break;
		iter ++;
	}	

	// calculate log likelihood of full model
	double p = this->get_prop();
	log_likely_1 = 0;
	for (int i=0; i<n_mol; i++)
		log_likely_1 += (1 - p) * f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]) + p * f1_log(ipd_avg[i], ipd_var[i], ipd_n[i]);
	
	// calculate log likelihood of null model
	
	log_likely_0 = 0;
	for (int i=0; i<n_mol; i++){
		log_likely_0 += f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]);
	}
	
	// BIC
	double sample_size = 0;
	for (int i=0; i<n_mol; i++){
		sample_size += ipd_n[i];
	}
	BIC_1 = -2 * log_likely_1 + 2 * log(sample_size);
	BIC_0 = -2 * log_likely_0 + log(sample_size);
	
	// AIC
	AIC_1 = -2 * log_likely_1 + 2*2;
	AIC_0 = -2 * log_likely_0 + 2;	

	return true;
}


