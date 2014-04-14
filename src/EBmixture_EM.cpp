#include "EBmixture_EM.h"

bool EBmixture_EM::getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx)
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

vector<int> EBmixture_EM::bin_search(double query, double *temp, int temp_len)
{
	vector<int> idx_pair(2,-1);
	if (temp_len<2) return idx_pair;
	temp_len--;	
	idx_pair[0] = 0;
	idx_pair[1] = temp_len;

	if (query<temp[0]){
		idx_pair[1] = 0;
		return idx_pair;
	}	
	if (query>temp[temp_len]){
		idx_pair[0] = temp_len;
                return idx_pair;
	}
	// assume temp is sorted (small to large)
	while (1){
		int cur_idx = idx_pair[0] + int((idx_pair[1] - idx_pair[0])/2.0);
	if (query <= temp[cur_idx]) 
			idx_pair[1] = cur_idx;
		else 
			idx_pair[0] = cur_idx;		

		if (idx_pair[1]-idx_pair[0]==1 || idx_pair[0]==temp_len || idx_pair[1]==0 )
			break;	

	}
	return idx_pair;
}

double EBmixture_EM::f0(double x_avg, double x_var, double x_n)
{
	return exp( -( x_n* ( (x_avg - mu_0)*(x_avg - mu_0) + x_var) ) / (2.0*sigma_0*sigma_0) ) / pow(sqrt(2.0*my_pi)*sigma_0, x_n);	

}
double EBmixture_EM::f0_log(double x_avg, double x_var, double x_n)
{
	return -( x_n* ( (x_avg - mu_0)*(x_avg - mu_0) + x_var) ) / (2.0*sigma_0*sigma_0) - x_n*log(2.0*my_pi) / 2.0 - x_n*log(sigma_0);
}

double EBmixture_EM::f1(double x_avg, double x_var, double x_n, double _mu_1)
{
        return exp( -( x_n* ( (x_avg - _mu_1)*(x_avg - _mu_1) + x_var) ) / (2.0*sigma_1*sigma_1) ) / pow(sqrt(2.0*my_pi)*sigma_1, x_n);

}
double EBmixture_EM::f1_log(double x_avg, double x_var, double x_n, double _mu_1)
{
        return -( x_n* ( (x_avg - _mu_1)*(x_avg - _mu_1) + x_var) ) / (2.0*sigma_1*sigma_1) - x_n*log(2.0*my_pi) / 2.0 - x_n*log(sigma_1);
}


double EBmixture_EM::f1_prior(double x)
{
	int len = f1_x.size();
	vector<int> idx_pair = bin_search(x, &f1_x[0], len);
	
	if (idx_pair[1]==0) return f1_y[0];
	if (idx_pair[0]==len-1) return f1_y[len-1];
        
	return f1_y[idx_pair[0]]+(x - f1_x[idx_pair[0]])* 
			(f1_y[idx_pair[1]] - f1_y[idx_pair[0]]) / (f1_x[idx_pair[1]] - f1_x[idx_pair[0]]);		

}

double EBmixture_EM::f1_prior_log(double x)
{
	double cur_f1 = f1_prior(x);
	if (cur_f1 < ERR) cur_f1 = ERR;
	return log(cur_f1);
}

bool EBmixture_EM::run(bool is_cal_mu_1, double prop_min)
{
	clear();
	mu_1 = sqrt(-1);
	mu_d = sqrt(-1);
	
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
		// update gamma_0 and gamma_1
		for (int i=0; i<n_mol; i++){
			double cur_f0_log = f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]);	
			if (ipd_avg[i] <= mu_0 && cur_f0_log <= -12){
                                gamma_0[i] = 1;
                                gamma_1[i] = 0;
                                continue;
                        }
                        if (ipd_avg[i] > mu_0 && cur_f0_log <= -12){
                                gamma_0[i] = 0;
                                gamma_1[i] = 1;
                                continue;
                        }


			double lik_r = 0;
			int cur_len = f1_x.size();
			for (int j=0; j<cur_len; j++){
				if (f1_y[j] < 1e-5)
					continue;
				double cur_lik_r_log = f1_log(ipd_avg[i], ipd_var[i], ipd_n[i],f1_x[j]) - cur_f0_log;	
				double cur_f1_y_log = log(f1_y[j]);
				if (cur_lik_r_log + cur_f1_y_log >= -9.21034)
					lik_r +=  exp(cur_lik_r_log + cur_f1_y_log);
			}
			lik_r *= s;
			gamma_0[i] = (1 - prop[iter-1]) / ( prop[iter-1]*lik_r + 1 - prop[iter-1] );
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
	
	// calculate mu_1
	
	if (is_cal_mu_1==true && this->get_prop() >= prop_min){
		mu_1 = 0;
		double C_norm = 0;
		for (int i = 0;i < n_mol; i++){
			mu_1 += gamma_1[i]*ipd_n[i]*ipd_avg[i];
			C_norm += gamma_1[i]*ipd_n[i];
		}		
		mu_1 = mu_1 / C_norm;
		mu_d = mu_1 - mu_0;
	}
	return true;
}


