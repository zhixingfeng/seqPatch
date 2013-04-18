normalizeBySubread <- function(alnsF, trim=0.1)
{
	for (i in 1:length(alnsF)){
		if (i==1000*floor(i/1000))
			cat('normalized ',i,' subreads\r')
		alnsF[[i]]$IPD <- alnsF[[i]]$IPD - mean(alnsF[[i]]$IPD, na.rm=TRUE, trim=trim)		
	}
	cat('normalized ',length(alnsF),' subreads\n')
	return(alnsF)
}


