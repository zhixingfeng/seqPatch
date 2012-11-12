#ifndef HIEMODEL_ALT_H
#define HIEMODEL_ALT_H
#include "hieModelEB.h"

double evalLogLikeliHood_alt(double mu0, double kappa0, double upsilon0, double sigma0,
                        vector<double> &mu_n, vector<double> &kappa_n, vector<double> &upsilon_n, vector<double> &sigma_n,
			 double upsilon_n_native, double sigma_n_native)
{
    	double k = mu_n.size();
    	// first term, mu0 and kappa0
    	double E_l_mu0_kappa0 = 0;
    	for (int i=0;i< k;i++){
        	E_l_mu0_kappa0 += 1/kappa_n[i] + (mu_n[i] - mu0)*(mu_n[i] - mu0) / sigma_n[i];
    	}
    	E_l_mu0_kappa0 = - kappa0*E_l_mu0_kappa0/2 + k*log(kappa0)/2;

    	// second term, upsilon0 and sigma0
    	double E_l_upsilon0_sigma0 = 0;
    	for (int i=0;i<k;i++){
        	E_l_upsilon0_sigma0 += log(upsilon_n[i]*sigma_n[i]/2) - digamma(upsilon_n[i]/2);
    	}	
    	E_l_upsilon0_sigma0 += log(upsilon_n_native*sigma_n_native/2) - digamma(upsilon_n_native/2);
    	E_l_upsilon0_sigma0 = (k+1)*upsilon0*log(upsilon0/2)/2 + (k+1)*upsilon0*(log(sigma0)-1)/2
        	                    - upsilon0* E_l_upsilon0_sigma0/2  - (k+1)* lgamma(upsilon0/2);
    	return E_l_mu0_kappa0 + E_l_upsilon0_sigma0;
}

vector<double> hyperParaEstimate_alt (double y_bar_native, double s2_native, double n_native, vector<double> &y_bar, vector<double> &s2,
				 vector <double> &n, double mu0, double kappa0, double upsilon0, double sigma0, int max_iter=10,
				 double min_change=0.02)
{
    	int k = (int) y_bar.size();

    	// get parameter of posterior distribution of mu and sigma given data and hyper parameters
    	vector<double> mu_n(k,0);
    	vector<double> kappa_n(k,0);
    	vector<double> upsilon_n(k,0);
    	vector<double> sigma_n(k,0);
    	for (int j=0;j<k;j++){
        	kappa_n[j] = kappa0 + n[j];
        	mu_n[j] = (kappa0*mu0+ n[j]*y_bar[j]) / kappa_n[j];
        	upsilon_n[j] = upsilon0 + n[j];
        	sigma_n[j] = (upsilon0*sigma0 + (n[j]-1)*s2[j] +
                      	kappa0*n[j]*(y_bar[j]-mu0)*(y_bar[j]-mu0)/kappa_n[j]) / upsilon_n[j];
    	}
	double upsilon_n_native = upsilon0 + n_native;
	double sigma_n_native = (upsilon0*sigma0 + (n_native-1)*s2_native)/upsilon_n_native;
		
    	// estimate hyper parameters by EM
    	for (int i=0; i< max_iter; i++){
        	// get parameter of posterior distribution of mu and sigma given data and hyper
        	for (int j=0;j<k;j++){
            		kappa_n[j] = kappa0 + n[j];
            		mu_n[j] = (kappa0*mu0+ n[j]*y_bar[j]) / kappa_n[j];
            		upsilon_n[j] = upsilon0 + n[j];
            		sigma_n[j] = (upsilon0*sigma0 + (n[j]-1)*s2[j] +
                         	 kappa0*n[j]*(y_bar[j]-mu0)*(y_bar[j]-mu0)/kappa_n[j]) / upsilon_n[j];
        	}
		upsilon_n_native = upsilon0 + n_native;
        	sigma_n_native = (upsilon0*sigma0 + (n_native-1)*s2_native)/upsilon_n_native;
        	
		// get original likelihood
                double L_og = evalLogLikeliHood_alt(mu0, kappa0, upsilon0, sigma0,
                                        mu_n, kappa_n, upsilon_n, sigma_n, upsilon_n_native, sigma_n_native);

		// estimate hyper parameters
        	double Z = 0;
        
		// mu0
        	mu0 = 0;
        	for (int j=0;j<k;j++){
            		mu0 += mu_n[j] / sigma_n[j];
            		Z += 1 / sigma_n[j];
        	}
        	mu0 = mu0 / Z;

        	// kappa0
        	kappa0 = 0;
        	for (int j=0;j<k;j++)
            		kappa0 += 1/kappa_n[j] + (mu_n[j] - mu0)*(mu_n[j] - mu0)/sigma_n[j];
        	kappa0 = k / kappa0;

        	// sigma0
        	sigma0 = 0;
        	for (int j=0;j<k;j++) sigma0 += 1/sigma_n[j];
		sigma0 += 1/sigma_n_native;
        	sigma0 = (k+1) / sigma0;

        	// upsilon0
       		double T = 0;
        	for (int j=0;j<k;j++) T += log(upsilon_n[j]*sigma_n[j]/2) - digamma(upsilon_n[j]/2);
        	T += log(upsilon_n_native*sigma_n_native/2) - digamma(upsilon_n_native/2);
		T = T / (k+1) - log(sigma0);
        	if (T <=0)
            		upsilon0 = 1e-4;
        	else{
            		upsilon0 = 2 / ( 3* (sqrt(1+4*T/3) - 1) );
        	}
        	if (upsilon0 <4.001) upsilon0 = 4.001;

        	// get likelihood
        	double L = evalLogLikeliHood_alt(mu0, kappa0, upsilon0, sigma0,
				 mu_n, kappa_n, upsilon_n, sigma_n,  upsilon_n_native, sigma_n_native);
        	//printf("\nEM step %d : \n", i+1);
        	//printf("mu0 = %lf\nkappa0 = %lf\nupsilon0 = %lf\nsigma0 = %lf\n",mu0, kappa0, upsilon0, sigma0);
        	//printf("log likelyhood change is %lf - %lf = %lf \n", L, L_og, L- L_og);
        	if ( L- L_og < min_change)
            		break;
    	}

    	vector<double> result;
    	result.push_back(mu0);
    	result.push_back(kappa0);
    	result.push_back(upsilon0);
    	result.push_back(sigma0);
    	return result;
}
map<string, vector<double> > hieModelEB_alt(vector<double> & data_native, vector<double> & ipd_ref_mean, vector<double> & ipd_ref_var, vector<double> & ipd_ref_len, int max_iter=10)
{
    	vector <double> sampleSize = ipd_ref_len;
	double sampleSize_native = data_native.size();

    	// get mean and variance for each group (sufficient statistics)
    	vector<double> s2 = ipd_ref_var;
    	vector<double> y_bar = ipd_ref_mean;
    	double s2_native = var(data_native);
	double y_bar_native = mean(data_native);	

    	/*------------------ set initial values of hyper parameters ------------------*/
    	// mu0
    	double mu0 = mean(y_bar);
	
    	// kappa0
    	double kappa0 = 0 ;
    	for (unsigned int i=0; i<y_bar.size(); i++){
        	kappa0 += (y_bar[i] - mu0)*(y_bar[i] - mu0) / s2[i];
    	}
    	kappa0 = y_bar.size() / kappa0;

    	// sigma0
    	double sigma0 = 0;
    	for (unsigned int i=0;i<s2.size();i++){
        	sigma0 += 1/s2[i];
    	}
	sigma0 += 1/s2_native;
    	sigma0 = (s2.size()+1)/sigma0;
    
	// upsilon0
    	double sigma_mean = (sum(s2) + s2_native)/(s2.size()+1) ;
    	double upsilon0 = 2*sigma_mean/(sigma_mean - sigma0);
    	if (upsilon0 < 4.001 ) upsilon0 = 4.001;

    	//printf("Initial values : \n");
    	//printf("mu0 = %lf\nkappa0 = %lf\nupsilon0 = %lf\nsigma0 = %lf\n",mu0, kappa0, upsilon0, sigma0);

    	/*---------------estimate hyper parameters----------------*/
    	vector<double> hyperPara = hyperParaEstimate_alt(y_bar_native, s2_native, sampleSize_native, 
						y_bar, s2, sampleSize, mu0, kappa0, upsilon0, sigma0, max_iter);
	
    	mu0 = hyperPara[0];
    	kappa0 = hyperPara[1];
    	upsilon0 = hyperPara[2];
    	sigma0 = hyperPara[3];

    	/*----------------get posterior mean of parameters ---------------------*/
    	int k = (int) y_bar.size();
    	vector<double> n = sampleSize;
	vector<double> mu_n(k,0);
    	vector<double> kappa_n(k,0);
    	vector<double> upsilon_n(k,0);
    	vector<double> sigma_n(k,0);
    	for (int i=0;i<k;i++){
        	kappa_n[i] = kappa0 + n[i];
        	mu_n[i] = (kappa0*mu0+ n[i]*y_bar[i]) / kappa_n[i];
        	upsilon_n[i] = upsilon0 + n[i];
        	sigma_n[i] = (upsilon0*sigma0 + (n[i]-1)*s2[i] +
                	      kappa0*n[i]*(y_bar[i]-mu0)*(y_bar[i]-mu0)/kappa_n[i]) / upsilon_n[i];
    	}

    	// get marginal likelihood given hyperparameters
    	map<string, vector<double> > result;
    	result["hyperPara"] = hyperPara;
    	result["mu"] = mu_n;
    	result["kappa"] = kappa_n;
    	result["upsilon"] = upsilon_n;
	result["sigma"] = sigma_n;
	return result;
}


#endif
