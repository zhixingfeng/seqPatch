#include "R.h"
#include "Rinternals.h"
#include "stl.h"
#include "math_utils.h"
#include <math.h>
#include "SpecialFuns.h"
double get_t_stat(double mu_x, double s2_x, double n_x, double mu_y, double s2_y, double n_y)
{
	double S = sqrt(s2_x/n_x + s2_y/n_y);
        return (mu_x-mu_y)/S;	
}

inline vector<double> dnorm(const vector<double> &x, double mu, double sigma2)
{
	vector<double> y;
        for (int i=0;i<(int)x.size();i++)
                y.push_back( exp(-(x[i]-mu)*(x[i]-mu)/(2*sigma2))/sqrt(2*my_pi*sigma2)  );
        return y;

}

inline vector<double> ldnorm(const vector<double> &x, double mu, double sigma2)
{
	vector<double> y;
	for (int i=0;i<(int)x.size();i++)
		y.push_back(-(x[i]-mu)*(x[i]-mu)/(2*sigma2) - log(2*my_pi*sigma2)/2);
        return y;
}

inline vector<double> rbinom(int N, int size, int prob)
{
		
}

