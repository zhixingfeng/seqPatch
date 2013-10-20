#ifndef SAMPLING_H
#define SAMPLING_H  

#define sample_SELECT RNG_sample
#include <R.h>
#include <Rmath.h>
#include "stl.h"

vector<double> RNG_rbinom(int n, double size, double prob);

vector<double> RNG_gamma(int size, int shape, double rate);

vector<double> RNG_normal(int size, double mu, double sigma);

vector<double> sample_SELECT (vector<double> &data, int M = 0);

vector<int> sample_SELECT (vector<int> &data, int M = 0);

vector<vector<double> > sample_group (vector<vector<double> > data_group);


#endif

