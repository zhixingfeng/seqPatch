normalizeBySubRead.context <- function(alnsF, alnsIdx, left.len=6, right.len=1,len.cutoff=200, K=1, max.iter=1)
{
	# add 0.01 to IPD
	for (i in 1:length(alnsF)){
		if (i == 1000*floor(i/1000) )
			cat('modified ',i, ' subreads\r')
		alnsF[[i]]$IPD <- alnsF[[i]]$IPD + 0.01
	}
	cat('modified ',length(alnsF), ' subreads\n')
	
        context.effect.all <- getContextEffect(alnsF, left.len, right.len)				
	context.effect <- lapply(context.effect.all, function(x) {x<-filter.outlier(x); exp(mean(log(x))) })
	len <- sapply(context.effect.all, length)
	context.effect.all.trim <- context.effect.all[len>=len.cutoff]
	context.effect.trim <- context.effect[len>=len.cutoff]
	#context.effect.trim.mean <- exp(mean(log(as.numeric(context.effect.trim))))
	#context.effect.trim <- lapply(context.effect.trim, function(x, ratio){ x*ratio  } , ratio=K/context.effect.trim.mean )
	
	alnsF.norm <- correctContextEffect(alnsF, context.effect.trim, left.len, right.len)			
		
	for (i in 1:max.iter){
		# remove subread effect
		subread.effect <- sapply(alnsF.norm, function(x){ x<-filter.outlier(x$IPD[!is.na(x$IPD)]);exp(mean(log(x)))  } )
		subread.effect.mean <- exp(mean(log(subread.effect)))
		subread.effect <- sapply(subread.effect, function(x, ratio) x*ratio, ratio=K/subread.effect.mean )
		alnsF.norm <- alnsF
		for (j in 1:length(alnsF)){
			if (j == 1000*floor(j/1000) )
                        	cat('corrected ',j, ' subreads\r')
			alnsF.norm[[j]]$IPD <- alnsF[[j]]$IPD / subread.effect[j] 	
		}
		cat('corrected ',length(alnsF), ' subreads\n')

		# remove context effect
		context.effect.all <- getContextEffect(alnsF.norm, left.len, right.len)
        	context.effect <- lapply(context.effect.all, function(x) {x<-filter.outlier(x); exp(mean(log(x))) })
        	len <- sapply(context.effect.all, length)
        	context.effect.all.trim <- context.effect.all[len>=len.cutoff]
        	context.effect.trim <- context.effect[len>=len.cutoff]
		#context.effect.trim.mean <- exp(mean(log(as.numeric(context.effect.trim))))
        	#context.effect.trim <- lapply(context.effect.trim, function(x, ratio){ x*ratio  } ,ratio=K/context.effect.trim.mean )
                alnsF.norm <- correctContextEffect(alnsF, context.effect.trim, left.len, right.len)
	}	
	subread.effect <- sapply(alnsF.norm, function(x){ x<-filter.outlier(x$IPD[!is.na(x$IPD)]);exp(mean(log(x)))  } )
	subread.effect.mean <- exp(mean(log(subread.effect)))
        subread.effect <- sapply(subread.effect, function(x, ratio) x*ratio, ratio=K/subread.effect.mean )
	for (j in 1:length(alnsF)){
        	if (j == 1000*floor(j/1000) )
                	cat('corrected ',j, ' subreads\r')
                alnsF.norm[[j]]$IPD <- alnsF[[j]]$IPD / subread.effect[j]
        }
        cat('corrected ',length(alnsF), ' subreads\n')

	result <- list()
	result$alnsF <- alnsF.norm
	result$context.effect <- context.effect.trim	
	result$context.effect.all <- context.effect.all.trim
	result
}

normalizeBySubRead.context.pred <- function(alnsF, context.effect, left.len=6, right.len=1)
{
	if (nchar(names(context.effect)[1]) != left.len + right.len + 1){
		stop('context lendth should be equal to left.len + rigth.len + 1')
	}
	alnsF.norm <- correctContextEffect(alnsF, context.effect, left.len, right.len)
	subread.effect <- sapply(alnsF.norm, function(x){ x<-filter.outlier(x$IPD[!is.na(x$IPD)]);exp(mean(log(x)))  } )
        alnsF.norm <- alnsF
        for (j in 1:length(alnsF)){
        	if (j == 1000*floor(j/1000) )
                	cat('corrected ',j, ' subreads\r')
                alnsF.norm[[j]]$IPD <- alnsF[[j]]$IPD / subread.effect[j]
        }
        cat('corrected ',length(alnsF), ' subreads\n')
	alnsF.norm
}




 
