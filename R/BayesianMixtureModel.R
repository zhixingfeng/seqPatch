BayesianMixtureModel <- function(IPD, idx, theta_0=0, kappa_0=0.001, upsilon_0=4.001, tau2_0=1000, 
			theta_1=theta_0+2, kappa_1=0.001, upsilon_1=4.001, tau2_1=1000, max.iter=10000)
{
	.Call('R_API_BayesianMixtureModel',as.numeric(IPD), as.numeric(idx), as.numeric(theta_0),
			 as.numeric(kappa_0), as.numeric(upsilon_0), as.numeric(tau2_0), 
			as.numeric(theta_1), as.numeric(kappa_1), as.numeric(upsilon_1),
			 as.numeric(tau2_1), as.integer(max.iter))	

	
}






