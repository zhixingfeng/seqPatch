mergeROC <- function(ROC.B, TPR.knots=seq(0.05,1,0.05))
{
	# spline interpolate
        model <- list()
        FDR.B <- list()
        for (i in 1:length(ROC.B)){
                model[[i]] <- list() 
                FDR.B[[i]] <- list()
                # fit spline for forward strand
                for (k in 1:length(ROC.B[[i]])){
                        TPR <- ROC.B[[i]][[k]]$TPR
                        FDR <- ROC.B[[i]][[k]]$FDR
                        chr <- names(ROC.B[[i]])[k]
			model[[i]][[chr]] <- approx(TPR,FDR, xout=TPR.knots, rule=2)
                        FDR.B[[i]][[chr]] <- model[[i]][[chr]]$y
                }
        }
	FDR.mean <-list()
        FDR.sd <-list()
        for (k in 1:length(FDR.B[[i]]) ){
                chr <- names(FDR.B[[i]])[k]
		FDR.mean[[chr]] <- 0
                FDR.sd[[chr]] <- 0
                for (i in 1:length(FDR.B)){
                        FDR.mean[[chr]] <- FDR.mean[[chr]] + FDR.B[[i]][[k]]
                        FDR.sd[[chr]] <- FDR.sd[[chr]] + FDR.B[[i]][[k]]^2
                }
                FDR.mean[[chr]] <- FDR.mean[[chr]] / length(FDR.B)
                FDR.sd[[chr]] <- sqrt(FDR.sd[[chr]]/length(FDR.B) - FDR.mean[[chr]]^2 + 1e-10)/sqrt(length(FDR.B))
                FDR.mean[[chr]][FDR.mean[[chr]]<0] <- 0
                FDR.mean[[chr]][FDR.mean[[chr]]>1] <- 1

        }
	list(FDR.B=FDR.B, FDR.mean=FDR.mean, FDR.sd= FDR.sd)
}




mergeROC.bak <- function(ROC.B, TPR.knots=seq(0.05,1,0.05))
{
	# spline interpolate
	model <- list()
	FDR.B <- list()
	for (i in 1:length(ROC.B)){
		model[[i]] <- list()
		FDR.B[[i]] <- list()
		# fit spline for forward strand
		for (k in 1:length(ROC.B[[i]]$pos)){
			TPR <- ROC.B[[i]]$pos[[k]]$TPR
			FDR <- ROC.B[[i]]$pos[[k]]$FDR
			model[[i]]$pos[[k]] <- approx(TPR,FDR, xout=TPR.knots, rule=2)
			FDR.B[[i]]$pos[[k]] <- model[[i]]$pos[[k]]$y
		}
		names(model[[i]]$pos) <- names(ROC.B[[i]]$pos)
		names(FDR.B[[i]]$pos) <- names(ROC.B[[i]]$pos)
		
		# fit spline for backward strand
                for (k in 1:length(ROC.B[[i]]$neg)){
                        TPR <- ROC.B[[i]]$neg[[k]]$TPR
                        FDR <- ROC.B[[i]]$neg[[k]]$FDR
			model[[i]]$neg[[k]] <- approx(TPR, FDR, xout=TPR.knots, rule=2)
			FDR.B[[i]]$neg[[k]] <- model[[i]]$neg[[k]]$y
                }
                names(model[[i]]$neg) <- names(ROC.B[[i]]$neg)
		names(FDR.B[[i]]$neg) <- names(ROC.B[[i]]$neg)
		
	}

	FDR.mean.pos <-list()
	FDR.sd.pos <-list()
	FDR.mean.neg <-list()
        FDR.sd.neg <-list()
	for (k in 1:length(ROC.B[[i]]$pos) ){
		FDR.mean.pos[[k]] <- 0
		FDR.sd.pos[[k]] <- 0
		for (i in 1:length(ROC.B)){
			FDR.mean.pos[[k]] <- FDR.mean.pos[[k]] + FDR.B[[i]]$pos[[k]]
			FDR.sd.pos[[k]] <- FDR.sd.pos[[k]] + FDR.B[[i]]$pos[[k]]^2	
		}
		FDR.mean.pos[[k]] <- FDR.mean.pos[[k]] / length(ROC.B)
		FDR.sd.pos[[k]] <- sqrt(FDR.sd.pos[[k]]/length(ROC.B) - FDR.mean.pos[[k]]^2 + 1e-10)/sqrt(length(ROC.B))
		FDR.mean.pos[[k]][FDR.mean.pos[[k]]<0] <- 0
	        FDR.mean.pos[[k]][FDR.mean.pos[[k]]>1] <- 1

	}

	for (k in 1:length(ROC.B[[i]]$neg) ){
                FDR.mean.neg[[k]] <- 0
                FDR.sd.neg[[k]] <- 0
                for (i in 1:length(ROC.B)){
                        FDR.mean.neg[[k]] <- FDR.mean.neg[[k]] + FDR.B[[i]]$neg[[k]]
                        FDR.sd.neg[[k]] <- FDR.sd.neg[[k]] + FDR.B[[i]]$neg[[k]]^2
                }
                FDR.mean.neg[[k]] <- FDR.mean.neg[[k]] / length(ROC.B)
                FDR.sd.neg[[k]] <- sqrt(FDR.sd.neg[[k]]/length(ROC.B) - FDR.mean.neg[[k]]^2 + 1e-10)/sqrt(length(ROC.B))
        	FDR.mean.neg[[k]][FDR.mean.neg[[k]]<0] <- 0
                FDR.mean.neg[[k]][FDR.mean.neg[[k]]>1] <- 1

	}
	names(FDR.mean.pos) <- names(ROC.B[[1]]$pos)
        names(FDR.sd.pos) <- names(ROC.B[[1]]$pos)
        names(FDR.mean.neg) <- names(ROC.B[[1]]$neg)
        names(FDR.mean.neg) <- names(ROC.B[[1]]$neg)
	
	result <- list()
	result$FDR.mean.pos <- FDR.mean.pos
        result$FDR.sd.pos <- FDR.sd.pos
        result$FDR.mean.neg <- FDR.mean.neg
        result$FDR.sd.neg <- FDR.sd.neg
	result$FDR.B <- FDR.B
	result
}



