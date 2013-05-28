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

normalizeBySubread.original <- function(alnsF, trim=0.1)
{
        for (i in 1:length(alnsF)){
                if (i==1000*floor(i/1000))
                        cat('normalized ',i,' subreads\r')
		alnsF[[i]]$IPD <- alnsF[[i]]$IPD / median(alnsF[[i]]$IPD,na.rm=TRUE)
		
		#IPD.log <- log(alnsF[[i]]$IPD + 0.0001)
                #alnsF[[i]]$IPD <- exp(IPD.log - mean(IPD.log, na.rm=TRUE, trim=trim)) - 0.0001
        }
        cat('normalized ',length(alnsF),' subreads\n')
        return(alnsF)
}


