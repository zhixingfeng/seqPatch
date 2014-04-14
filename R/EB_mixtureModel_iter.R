detectModPropEB.iter <- function(genomeF.native, genomeF.wga, f1.x, f1.y, is_f1_var=FALSE, max.iter = 100, is.EM=TRUE, d.eps=0.05, n.round=20)
{

	if (n.round < 2)
		stop('n.round should be >= 2')
	d.median <- numeric()
	for (i in 1:n.round){
		cat('iter ',i,':\n', sep='')
        	detect <- detectModPropEB(genomeF.native, genomeF.wga, f1.x, f1.y, is_f1_var=is_f1_var, max.iter=max.iter, is.EM=is.EM)
        	d <- detect$pos$chr$mu_d[!is.na(detect$pos$chr$mu_d)]
        	d.median[i] <- median(d,na.rm=TRUE)

        	f1.hist <- hist(d, breaks=50, plot=FALSE)
        	f1.x <- f1.hist$mids
        	f1.y <- f1.hist$density		
		if (i>1)
			if(abs(d.median[i] - d.median[i-1])<=d.eps)
				break
	}
	if (abs(d.median[i] - d.median[i-1]) > d.eps)
		warning('prior estimation does not converge, try to increase n.round')

	detect <- detectModPropEB(genomeF.native, genomeF.wga, f1.x, f1.y, is_f1_var=is_f1_var, max.iter=max.iter, is.EM=is.EM)
	
	rl <- list()
	rl$detect <- detect
	rl$f1.x <- f1.x
	rl$f1.y <- f1.y
	rl$d <- d
	rl$d.median <- d.median
	rl$n.iter <- i
	rl
}



