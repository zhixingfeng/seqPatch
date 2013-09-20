#include "EBmixture.h"
bool EBmixture::getMoleculeMeanIPD(double *ipd, double *idx, int len_ipd, int len_idx)
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

vector<int> EBmixture::bin_search(double query, double *temp, int temp_len)
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
double EBmixture::f1_log_int(double x_avg, double x_var, double x_n)
{
	int len = f_mu_1.size();
	double rl = 0;
	for (int i=0; i<len; i++){
		double mu_1 = f1_x[i];
		rl += f_mu_1[i] * (x_avg - mu_1)*(x_avg - mu_1);
	}
	rl *= s;
	rl = - (rl + x_var) * x_n / (2.0*sigma_1*sigma_1) - x_n*log(sigma_1) - x_n*log(2.0*my_pi)/2.0;	
			
	return rl;
}

double EBmixture::f0(double x_avg, double x_var, double x_n)
{
	//return exp(-(x - mu_0)*(x - mu_0) / (2*sigma_0*sigma_0)) / (sqrt(2*my_pi)*sigma_0);	
	return exp( -( x_n* ( (x_avg - mu_0)*(x_avg - mu_0) + x_var) ) / (2.0*sigma_0*sigma_0) ) / pow(sqrt(2.0*my_pi)*sigma_0, x_n);	

}
double EBmixture::f0_log(double x_avg, double x_var, double x_n)
{
	return -( x_n* ( (x_avg - mu_0)*(x_avg - mu_0) + x_var) ) / (2.0*sigma_0*sigma_0) - x_n*log(2.0*my_pi) / 2.0 - x_n*log(sigma_0);
	//return -(x - mu_0)*(x - mu_0) / (2*sigma_0*sigma_0) - log(2*my_pi)/2 - log(sigma_0);
}

double EBmixture::f1(double x)
{
	int len = f1_x.size();
	vector<int> idx_pair = bin_search(x, &f1_x[0], len);
	
	if (idx_pair[1]==0) return f1_y[0];
	if (idx_pair[0]==len-1) return f1_y[len-1];
        
	return f1_y[idx_pair[0]]+(x - f1_x[idx_pair[0]])* 
			(f1_y[idx_pair[1]] - f1_y[idx_pair[0]]) / (f1_x[idx_pair[1]] - f1_x[idx_pair[0]]);		

}

double EBmixture::f1_log(double x)
{
	double cur_f1 = f1(x);
	if (cur_f1 < ERR) cur_f1 = ERR;
	return log(cur_f1);
}

bool EBmixture::run()
{
	clear();
	// Initialize
	f_mu_1 = f1_y;		
	N_0.push_back(n_mol / 2.0);
	N_1.push_back(n_mol / 2.0);	
	for (int i=0; i<n_mol; i++){
		gamma_0.push_back(-1);
		gamma_1.push_back(-1);
	}
	prop.push_back(-1);
	// iterate
	double eps = 1e-4;	
	iter = 1;
	while (iter <= max_iter){
		// update delta
		for (int i = 0;i < n_mol; i++){
			double f0_log_z = f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]);
                	double f1_log_z = f1_log_int(ipd_avg[i], ipd_var[i], ipd_n[i]);
			double gamma_0_core = f1_log_z + digamma(N_1[iter-1] + 1) - f0_log_z - digamma(N_0[iter-1] + 1);
			gamma_0[i] = 1 / (1 + exp(gamma_0_core)) ;
	                if (gamma_0_core <= -12)
				gamma_0[i] = 1;
			if (gamma_0_core >= 12)
				gamma_0[i] = 0;
			gamma_1[i] = 1 - gamma_0[i];
		}			

		// update N_0 and N_1
		double cur_N_0 = sum(gamma_0);
		double cur_N_1 = sum(gamma_1);
		N_0.push_back(cur_N_0);
		N_1.push_back(cur_N_1);
		prop.push_back(cur_N_1 / (cur_N_0 + cur_N_1) );
		
	
		// update f_mu_1;
		/*int s_len = f_mu_1.size(); 
		for (int i=0; i < s_len; i++){
			double cur_f_mu_1 = 0;
			double cur_mu_1 = f1_x[i];
			for (int j=0; j<n_mol; j++){
				cur_f_mu_1 += gamma_1[j]*ipd_n[j]*(ipd_avg[j] - cur_mu_1)*(ipd_avg[j] - cur_mu_1);
			}			
			
			cur_f_mu_1 = exp(- cur_f_mu_1/(2*sigma_1*sigma_1) + log(f1_y[i]) );
			if (f1_y[i]<=1e-12) cur_f_mu_1 = 1e-12;
			f_mu_1[i] = cur_f_mu_1;
		}
		double f_mu_1_C = 0;
		for (int i=0; i<s_len; i++){
			f_mu_1_C += f_mu_1[i];
		}		
		f_mu_1_C *= s;
	
		for (int i=0; i<s_len; i++){
                        f_mu_1[i] = f_mu_1[i] / f_mu_1_C;
                }*/

		if (fabs(prop[iter] - prop[iter-1]) <= eps) break;
		iter ++;
	}	

	return true;
}



