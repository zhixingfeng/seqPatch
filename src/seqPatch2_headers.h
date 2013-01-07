#ifndef SEQPATH2_HEADERS_H
#define SEQPATH2_HEADERS_H

#include <Rcpp.h>
#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include "SpecialFuns.h"
#include "BayesianMixtureModel.h"

map<string, vector<double> > hieModelEB(vector<double> & ipd_ref_mean, vector<double> & ipd_ref_var,
                                vector<double> & ipd_ref_len, int max_iter=10);

string SP_reverseSeq(const string &seq);

int getIdxByName (SEXP data, string query);


using namespace std;
using namespace Rcpp;

#endif




