#ifndef POOLING_H
#define POOLING_H

#include "stl.h"
map<string, double> pooling(const vector<vector<double> >  & data)
{
	double mu = 0;
	double sigma = 0;
	double n = 0;
	for (int i=0 ;i < (int)data.size(); i++){
		for (int j=0;j < (int)data[i].size(); j++){
			mu += data[i][j];
			sigma += data[i][j]*data[i][j];
		}
		n += data[i].size();
	}			
	sigma = sigma/(n-1) - mu*mu/(n*n - n);
	mu = mu/n;	
	map<string, double> result;
	result["mu"] = mu;	
	result["sigma"] = sigma;
	result["sampleSize"] = n;
	return result;	

}

#endif


