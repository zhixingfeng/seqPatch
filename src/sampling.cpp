#include "sampling.h"

vector<double> RNG_rbinom(int n, double size, double prob)
{
	GetRNGstate();
	vector<double> y;
	for (int i=0;i<n;i++) y.push_back(Rf_rbinom(size, prob));
	PutRNGstate();
	return y;
}

vector<double> RNG_gamma(int size, int shape, double rate)
{	
	vector<double> x;
	if (shape<1){
		Rprintf("shape should be greater than 1");
		return x;
	}
	GetRNGstate();
	for (int i=0;i<size;i++){
		x.push_back(0);
		for (int j=0;j<shape;j++){
			x[i] += exp_rand()/rate; 
		}
	}	
	PutRNGstate();	
	return x;
}

vector<double> RNG_normal(int size, double mu, double sigma)
{
        vector<double> x;
        GetRNGstate();
        for (int i=0;i<size;i++)
                x.push_back( norm_rand()*sigma +mu );
        PutRNGstate();
        return x;
}

vector<double> sample_SELECT (vector<double> &data, int M)
{
	if (M==0) M = (int) data.size();
	int n = data.size();
	vector<double> result(M,sqrt(-1));
	if (M>n) {Rprintf("Error: M should be less than length of data and >=0.\n"); return result; }
	GetRNGstate();
	for (int i=0;i<M;i++){
		int idx = int(Rf_runif(0,n-i));
		result[i] = data[idx];
		data[idx] = data[n-i-1];
	}
	PutRNGstate();
	return result;
}
vector<int> sample_SELECT (vector<int> &data, int M)
{
        if (M==0) M = (int) data.size();
        int n = data.size();
        vector<int> result(M,sqrt(-1));
        if (M>n) {Rprintf("Error: M should be less than length of data and >=0.\n"); return result; }
        GetRNGstate();
        for (int i=0;i<M;i++){
                int idx = int(Rf_runif(0,n-i));
                result[i] = data[idx];
                data[idx] = data[n-i-1];
        }
        PutRNGstate();
        return result;
}

vector<vector<double> > sample_group (vector<vector<double> > data_group)
{
	vector<double> data_pool;
	for (int i=0;i<(int)data_group.size();i++)
		for (int j=0;j<(int)data_group[i].size();j++)
			data_pool.push_back(data_group[i][j]);

	vector<double> data_pool_perm = sample_SELECT(data_pool);
	
	vector<vector<double> > data_group_perm;
	int idx=0;
	for (int i=0;i<(int)data_group.size();i++){
        	vector<double> cur_data_perm;
	        for (int j=0;j<(int)data_group[i].size();j++)
			cur_data_perm.push_back(data_pool_perm[idx+j]);
		data_group_perm.push_back(cur_data_perm);
		idx += (int) cur_data_perm.size();
	}
	return data_group_perm;	
}



