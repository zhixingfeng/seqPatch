BayesianMixtureModel <- function(IPD, idx, theta_0=0, kappa_0=0.001, upsilon_0=4.001, tau2_0=1000, 
			theta_1=theta_0+2, kappa_1=0.001, upsilon_1=4.001, tau2_1=1000, max.iter=100)
{
	.Call('BayesianMixtureModel',IPD, idx, theta_0, kappa_0, upsilon_0, tau2_0, 
		theta_1, kappa_1, upsilon_1, tau2_1)	

	
}






