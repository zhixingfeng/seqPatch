#include "EBmixture.h"

bool operator<(const IPD_data& x, const IPD_data& y){ return x.avg < y.avg; }
bool operator<=(const IPD_data& x, const IPD_data& y){ return x.avg <= y.avg; }
bool operator>(const IPD_data& x, const IPD_data& y){ return x.avg > y.avg; }
bool operator>=(const IPD_data& x, const IPD_data& y){ return x.avg >= y.avg; }
bool operator==(const IPD_data& x, const IPD_data& y){ return x.avg == y.avg; }
bool operator!=(const IPD_data& x, const IPD_data& y){ return x.avg != y.avg; }




map<string, vector<double> >  EBmixture::subsample(double *ipd, double *idx, int len_ipd, int len_idx, double rate)
{
	map<string, vector<double> > rl;

	if (len_ipd != len_idx){ printf("length of R_IPD and R_idx should be the same.\n"); return rl;}
	
	ipd_map.clear();

	int len = len_ipd;
	for (int i=0; i<len; i++){
		ipd_map[int(idx[i])].push_back(ipd[i]);
	}

	map<int, vector<double> >::iterator it = ipd_map.begin();
	while(it!=ipd_map.end()){
		int M = int(rate*it->second.size());
		if (M < 1) M = 1;
		it->second = sample_SELECT(it->second, M);
		double cur_n = it->second.size();
		
		for (int i=0; i<cur_n; i++){
			rl["IPD"].push_back(it->second[i]);
			rl["idx"].push_back(double(it->first));	
		}
	
		it++;
	}
		
	return rl;
}

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

bool EBmixture::run(bool is_cal_mu_1, bool is_f1_varible, double prop_min)
{
	clear();
	mu_1 = sqrt(-1);
	mu_d = sqrt(-1);
	//ipd_sort();
	// Initialize
	f_mu_1 = f1_y;		
	N_0.push_back(n_mol / 2.0);
	N_1.push_back(n_mol / 2.0);	
	for (int i=0; i<n_mol; i++){
		gamma_0.push_back(-1);
		gamma_1.push_back(-1);
	}
	prop.push_back(-1);

	// integrate out mu1
	//vector<double> f1_log_z_all(n_mol, sqrt(-1));	
	//for (int i=0; i<n_mol; i++)
	//	f1_log_z_all[i] = f1_log_int(ipd_avg[i], ipd_var[i], ipd_n[i]);
	
	int s_len = f_mu_1.size();
	// iterate
	double eps = 1e-4;	
	iter = 1;
	while (iter <= max_iter){
		// update I
		
		//Rprintf("f0_log_z,f1_log_z,f1_log_z_norm, gamma_0_core:\n");

		for (int i = 0;i < n_mol; i++){
			double f0_log_z = f0_log(ipd_avg[i], ipd_var[i], ipd_n[i]);
                	double f1_log_z = f1_log_int(ipd_avg[i], ipd_var[i], ipd_n[i]);
			/*if (ipd_avg[i] <= mu_0 && (f1_log_z <= -12 || f0_log_z <= -12)){
				gamma_0[i] = 1;
				gamma_1[i] = 0;
				continue;
			}	
			if (ipd_avg[i] > mu_0 && f0_log_z <= -12){
                                gamma_0[i] = 0;
                                gamma_1[i] = 1;
				continue;
                        }*/
		
			double gamma_0_core = f1_log_z + digamma(N_1[iter-1] + 1) - f0_log_z - digamma(N_0[iter-1] + 1);
			gamma_0[i] = 1 / (1 + exp(gamma_0_core)) ;
	                if (gamma_0_core <= -12)
				gamma_0[i] = 1;
			if (gamma_0_core >= 12)
				gamma_0[i] = 0;
			gamma_1[i] = 1 - gamma_0[i];		
		}	
				
		//Rprintf("\n");
		// update N_0 and N_1
		double cur_N_0 = sum(gamma_0);
		double cur_N_1 = sum(gamma_1);
		N_0.push_back(cur_N_0);
		N_1.push_back(cur_N_1);
		prop.push_back(cur_N_1 / (cur_N_0 + cur_N_1) );
			
		// update f_mu_1;
		if (is_f1_varible==true){		
			// (pre-compute coefficient determined by data)
			vector<double> coefs_log_sum;
			double coefs_log_sum_max = -1e24 ;
			//Rprintf("coefs_log_sum : ");
			for (int i=0; i < s_len; i++){
				double cur_mu_1 = f1_x[i];
				double cur_coefs_log_sum = 0;
				for (int j=0; j<n_mol; j++){
					double log_dens_fun = -( ipd_n[j]* ( (ipd_avg[j] - cur_mu_1)*(ipd_avg[j] - cur_mu_1) + ipd_var[j]) ) / (2.0*sigma_1*sigma_1)
							 - ipd_n[j]*log(2.0*my_pi) / 2.0 - ipd_n[j]*log(sigma_1);
					cur_coefs_log_sum += gamma_1[j]*log_dens_fun;
					//double dens_fun = exp(gamma_1[j] * log_dens_fun); 
				}
				//Rprintf("%.4lf ", cur_coefs_log_sum);
				coefs_log_sum.push_back(cur_coefs_log_sum);
				if (coefs_log_sum_max < cur_coefs_log_sum + (1+cur_N_1)*log(f1_y[i]))
				coefs_log_sum_max = cur_coefs_log_sum + (1+cur_N_1)*log(f1_y[i]);
			}	
			//Rprintf("\n");
			//Rprintf("coefs_log_sum_max : %.4lf\n",coefs_log_sum_max);
		
			//Rprintf("f_mu_1_unnorm : ");
 			for (int i=0; i < s_len; i++){
				//double cur_f_mu_1_log =  coefs_log_sum[i] + (1+cur_N_1)*log(f1_y[i]) - coefs_log_sum_max;
				double cur_f_mu_1_log =  coefs_log_sum[i] + log(f1_y[i]) - coefs_log_sum_max;
				double cur_f_mu_1;
				if (cur_f_mu_1_log <= -12 | f1_y[i]<=1e-12)
					cur_f_mu_1 = 0;
				else 
					cur_f_mu_1 = exp(cur_f_mu_1_log);
				//if (f1_y[i]<=1e-12) cur_f_mu_1 = 1e-12;
                        	f_mu_1[i] = cur_f_mu_1;	
				//Rprintf("%.4lf ", f_mu_1[i]);
			}
			//Rprintf("\n");
			double f_mu_1_C = 0;
			for (int i=0; i<s_len; i++){
				f_mu_1_C += f_mu_1[i];
			}		
			f_mu_1_C *= s;
		
			for (int i=0; i<s_len; i++){
                        	f_mu_1[i] = f_mu_1[i] / f_mu_1_C;
			}
		}
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


map<string, vector<int> > EBmixture::buildGenomeIndex(string & genomeSeq, int seed_len)
{
        map<string, vector<int> > genomeIndex;
        for (int i=0; (int) i<genomeSeq.size()-seed_len+1;i++){
                string cur_context = genomeSeq.substr(i,seed_len);
                genomeIndex[cur_context].push_back(i);
        }
        return genomeIndex;
}

vector<int> EBmixture::findMaxContext(int cur_idx, string & genomeSeq, 
					map<string, vector<int> > &genomeIndex, int left_len, int right_len)
{
	int lim = 100;
	vector<int> max_index;
	int context_start = cur_idx - left_len;
	int context_end = cur_idx + right_len;
	if (context_start<0 | context_end>= (int) genomeSeq.size())
		return max_index;
	string context = genomeSeq.substr(context_start, context_end - context_start + 1);	
	
	vector<int> ref_index = genomeIndex[context];
	int max_left_len = 0;
	vector<int> left_ex(ref_index.size(),0);
	for (int i=0; i<(int) ref_index.size(); i++){
		left_ex[i] = 0;
		if (context_start == ref_index[i])
			continue;
		for (int j=1; j<=lim; j++){
			if (context_start - j<0 | ref_index[i] - j<0)
				break;	
			if (genomeSeq[context_start - j] != genomeSeq[ref_index[i] - j])
				break;	
			left_ex[i]++;	
		}
		if (left_ex[i]>max_left_len) 
			max_left_len = left_ex[i];
	}
	
	//if (max_left_len == 0)
	//	return max_index;
	for (int i=0; i<(int) ref_index.size(); i++){
		if(left_ex[i] == max_left_len && ref_index[i]!=context_start)
			max_index.push_back(ref_index[i] + left_len);
	}
	return max_index;	
}


void EBmixture::ipd_sort()
{
	vector<IPD_data> ipd_data;
	for (int i=0;(int)i<ipd_avg.size();i++){
		IPD_data cur_ipd_data; 
		cur_ipd_data.avg = ipd_avg[i];
		cur_ipd_data.var = ipd_var[i];
		cur_ipd_data.n = ipd_n[i];
		ipd_data.push_back(cur_ipd_data);
	}
	sort(ipd_data.begin(), ipd_data.end());
	for (int i=0;(int)i<ipd_avg.size();i++){
              	ipd_avg[i] = ipd_data[i].avg;
                ipd_var[i] = ipd_data[i].var;
                ipd_n[i] = ipd_data[i].n;
        }

}


