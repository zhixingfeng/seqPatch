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
	n_subreads.clear();
	n_subreads = ipd_n;
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
double EBmixture::f0(double x)
{
	return exp(-(x - mu_0)*(x - mu_0) / (2*sigma_0*sigma_0)) / (sqrt(2*my_pi)*sigma_0);	
}
double EBmixture::f0_log(double x)
{
	return -(x - mu_0)*(x - mu_0) / (2*sigma_0*sigma_0) - log(2*my_pi)/2 - log(sigma_0);
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
		
	return true;
}



