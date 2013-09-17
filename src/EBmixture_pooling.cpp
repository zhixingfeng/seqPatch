#include "EBmixture_pooling.h"

vector<int> EBmixture_pooling::bin_search(double query, double *temp, int temp_len)
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
double EBmixture_pooling::f0(double x)
{
	return exp(-(x - mu_0)*(x - mu_0) / (2*sigma_0*sigma_0)) / (sqrt(2*my_pi)*sigma_0);	
}
double EBmixture_pooling::f0_log(double x)
{
	return -(x - mu_0)*(x - mu_0) / (2*sigma_0*sigma_0) - log(2*my_pi)/2 - log(sigma_0);
}


double EBmixture_pooling::f1(double x)
{
	int len = f1_x.size();
	vector<int> idx_pair = bin_search(x, &f1_x[0], len);
	
	if (idx_pair[1]==0) return f1_y[0];
	if (idx_pair[0]==len-1) return f1_y[len-1];
        
	return f1_y[idx_pair[0]]+(x - f1_x[idx_pair[0]])* 
			(f1_y[idx_pair[1]] - f1_y[idx_pair[0]]) / (f1_x[idx_pair[1]] - f1_x[idx_pair[0]]);		

}

double EBmixture_pooling::f1_log(double x)
{
	double cur_f1 = f1(x);
	if (cur_f1 < ERR) cur_f1 = ERR;
	return log(cur_f1);
}

bool EBmixture_pooling::run()
{
	clear();	
	// Initialize gamma_0 (simply f0 / (f1 + f0) ), gamma_1, rho_0, rho_1, N_0, N_1, and prop
	N_0.push_back(0);
	N_1.push_back(0);	
	for (int i = 0; i < (int)z.size(); i++){
		double f0_log_z	= f0_log(z[i]);
		double f1_log_z = f1_log(z[i]);
		double cur_gamma_0 = 1 / (1 + exp(f1_log_z - f0_log_z) );
		if (z[i]<=-6)
			cur_gamma_0 = 1;
		if (z[i]>=6)
			cur_gamma_0 = 0;
		double cur_gamma_1 = 1 - cur_gamma_0;
			
		gamma_0.push_back(cur_gamma_0);
		gamma_1.push_back(cur_gamma_1);
		N_0[0] += cur_gamma_0;
		N_1[0] += cur_gamma_1;	
		
		rho_0.push_back(1);
		rho_1.push_back(1);
	}	
	prop.push_back(N_1[0] / (N_0[0] + N_1[0]));	
	
	// interate
	double eps = 1e-4;	
	iter = 0;
	while (iter < max_iter){
		for (int i = 0;i < (int)z.size(); i++){
			double f0_log_z = f0_log(z[i]);
                	double f1_log_z = f1_log(z[i]);
			gamma_0[i] = 1 / (1 + exp(f1_log_z + digamma(N_1[iter] + 1) - f0_log_z - digamma(N_0[iter] + 1)) );
	                if (z[i] <= -6)
				gamma_0[i] = 1;
			if (z[i] >= 6)
				gamma_0[i] = 0;
			gamma_1[i] = 1 - gamma_0[i];
		}			
		
		iter++;
		
		double cur_N_0 = sum(gamma_0);
		double cur_N_1 = sum(gamma_1);
		N_0.push_back(cur_N_0);
		N_1.push_back(cur_N_1);
		prop.push_back(cur_N_1 / (cur_N_0 + cur_N_1) );
		if (fabs(prop[iter] - prop[iter-1]) <= eps) break;
	}	
	return true;
}



