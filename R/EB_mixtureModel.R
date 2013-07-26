EB.mixture.model <- function(z.score, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	.Call('R_API_EBmixture', as.numeric(z.score), as.numeric(mu.0), as.numeric(sigma.0),
		as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

}


