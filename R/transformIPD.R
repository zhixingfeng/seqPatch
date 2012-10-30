transformIPD <- function(alnsF, alpha = 0.005, lambda = 0.16)
{
	for (i in 1:length(alnsF)){
		alnsF[[i]]$IPD <- bc.transform(alnsF[[i]]$IPD+alpha, lambda)
		if (i==1000*floor(i/1000))cat('transformed ',i,' subreads\r')
	}
	cat('transformed ',length(alnsF),' subreads\n')
	alnsF
}


