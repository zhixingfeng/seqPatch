#ifndef SEQPATH2_HEADERS_H
#define SEQPATH2_HEADERS_H
#include "SpecialFuns.h"
#include "BayesianMixtureModel.h"

map<string, vector<double> > hieModelEB(vector<double> & ipd_ref_mean, vector<double> & ipd_ref_var,
                                vector<double> & ipd_ref_len, int max_iter=10);


#endif




