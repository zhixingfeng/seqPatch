getROC.m6A.pin <- function(detection, idx.methy, genomeSeq,  method='effect.size', tol=20, idx.sel=NULL)
{
        ROC.pos <- list(); ROC.neg <- list()
        for (i in 1:length(detection$pos)){
                refname <- names(detection$pos)[i]
                if (method=='LR_log')
			t.stat <- detection$pos[[i]][[method]]
		else
			t.stat <- abs(detection$pos[[i]][[method]])
                if (is.null(idx.sel)) {
                        idx.sub <- which(!is.na(t.stat))
                }else {
                        idx.sub <- idx.sel$pos[[refname]]
                }
                # mask regions don't contain GATC
                cur.idx.methy.tol <- numeric(0)
                for (j in -tol:tol) cur.idx.methy.tol <- c(cur.idx.methy.tol, idx.methy$pos[[refname]]+j)
                cur.idx.methy.tol <- unique(cur.idx.methy.tol)
                idx.sub <- intersect(idx.sub, cur.idx.methy.tol)

                cutoff <- quantile(t.stat[intersect(idx.methy$pos[[refname]],idx.sub) ], prob=seq(0.95,0.05,-0.05))

                TPR <- seq(0.05,0.95,0.05)
                FDR <- numeric(length(TPR))
                idx.A <- gregexpr('T',genomeSeq$pos[[refname]])[[1]]
                idx.sub <- intersect(idx.sub, idx.A)
		idx.null <- setdiff(idx.sub, idx.methy$pos[[refname]])
               	idx.null <- setdiff(idx.null, idx.methy$pos[[refname]]+5) 
		for (j in 1:length(TPR)){
                        FDR[j] <- sum(t.stat[idx.null]>=cutoff[j]) / sum(t.stat[idx.sub]>=cutoff[j])
                }
		ROC.pos[[refname]]$TPR <- TPR
		ROC.pos[[refname]]$FDR <- FDR
        }

	for (i in 1:length(detection$neg)){
                refname <- names(detection$neg)[i]
		if (method=='LR_log')
                        t.stat <- detection$neg[[i]][[method]]
                else
                	t.stat <- abs(detection$neg[[i]][[method]])
                if (is.null(idx.sel)) {
                        idx.sub <- which(!is.na(t.stat))
                }else {
                        idx.sub <- idx.sel$neg[[refname]]
                }
                # mask regions don't contain GATC
                cur.idx.methy.tol <- numeric(0)
                for (j in -tol:tol) cur.idx.methy.tol <- c(cur.idx.methy.tol, idx.methy$neg[[refname]]+j)
                cur.idx.methy.tol <- unique(cur.idx.methy.tol)
                idx.sub <- intersect(idx.sub, cur.idx.methy.tol)

                cutoff <- quantile(t.stat[intersect(idx.methy$neg[[refname]],idx.sub) ], prob=seq(0.95,0.05,-0.05))

                TPR <- seq(0.05,0.95,0.05)
                FDR <- numeric(length(TPR))
                idx.A <- gregexpr('T',genomeSeq$neg[[refname]])[[1]]
                idx.sub <- intersect(idx.sub, idx.A)
		idx.null <- setdiff(idx.sub, idx.methy$neg[[refname]])
		idx.null <- setdiff(idx.null, idx.methy$pos[[refname]]-5)
                for (j in 1:length(TPR)){
                        FDR[j] <- sum(t.stat[idx.null]>=cutoff[j]) / sum(t.stat[idx.sub]>=cutoff[j])
                }
                ROC.neg[[refname]]$TPR <- TPR
                ROC.neg[[refname]]$FDR <- FDR
        }
	list(pos=ROC.pos, neg=ROC.neg)
}


getROC.meA <- function(detection, detection.null, idx.methy, idx.sel=NULL, method='effect.size')
{
        ROC.pos <- list(); ROC.neg <- list()
        for (i in 1:length(detection$pos)){
                refname <- names(detection$pos)[i]
		if (method=='LR_log'){
                        t.stat <- detection$pos[[i]][[method]]
			t.stat.null <- detection.null$pos[[i]][[method]]
		}else{
                	t.stat <- abs(detection$pos[[i]][[method]])
                	t.stat.null <- abs(detection.null$pos[[i]][[method]])
                }
		if (is.null(idx.sel)) {
                        idx.sub <- which(!is.na(t.stat) & !is.na(t.stat.null))
                }else {
                        idx.sub <- idx.sel$pos[[refname]]
                }
                cutoff <- quantile(t.stat[intersect(idx.methy$pos[[refname]],idx.sub) ], prob=seq(0.95,0.05,-0.05))
                TPR <- seq(0.05,0.95,0.05)
                FDR <- numeric(length(TPR))
                for (j in 1:length(TPR)){
                        FDR[j] <- sum(t.stat.null[idx.sub]>=cutoff[j]) / sum(t.stat[idx.sub]>=cutoff[j])
                }
                ROC.pos[[refname]]$TPR <- TPR
                ROC.pos[[refname]]$FDR <- FDR
        }

        for (i in 1:length(detection$neg)){
                refname <- names(detection$neg)[i]
		if (method=='LR_log'){
                        t.stat <- detection$neg[[i]][[method]]
                        t.stat.null <- detection.null$neg[[i]][[method]]
                }else{
                	t.stat <- abs(detection$neg[[i]][[method]])
                	t.stat.null <- abs(detection.null$neg[[i]][[method]])
                }
		if (is.null(idx.sel)){
                        idx.sub <- which(!is.na(t.stat) & !is.na(t.stat.null))
                }else{
                        idx.sub <- idx.sel$neg[[refname]]
                }
                cutoff <- quantile(t.stat[intersect(idx.methy$neg[[refname]],idx.sub) ], prob=seq(0.95,0.05,-0.05))
                TPR <- seq(0.05,0.95,0.05)
                FDR <- numeric(length(TPR))
                for (j in 1:length(TPR)){
                        FDR[j] <- sum(t.stat.null[idx.sub]>=cutoff[j]) / sum(t.stat[idx.sub]>=cutoff[j])
                }
                ROC.neg[[refname]]$TPR <- TPR
                ROC.neg[[refname]]$FDR <- FDR
        }
        list(pos=ROC.pos, neg=ROC.neg)
}

getROC.withNull <- function(detection, detection.null, idx.methy, tol=5)
{
	### get ROC for positive strand
	ROC <- list()
        for (i in 1:length(detection$pos)){
                chr <- names(detection$pos)[i]
                genome.start.pos <- detection$genome.start.pos[[chr]]
                if (is.null(genome.start.pos)){
                        cat('Error: can not find ',chr, ' in genome start position\n')
                        return (NULL)
                }
                t.stat <- abs(detection$pos[[i]]$t.stat)
		t.stat.null <- abs(detection.null$pos[[i]]$t.stat)
	
                ### set cutoff of t.stat and get TP and FP
                cur.idx.methy <- idx.methy$pos[[chr]]
                tmp <- getPvalueCutoff(t.stat, cur.idx.methy, tol)
                p.cutoff <- quantile(tmp$p.cutoff,prob=seq(0.95,0.05,-0.05),na.rm=TRUE)
                TP.total <- tmp$TP.total
                TPR <- seq(0.05,0.95,0.05)
                FD <- numeric(length(p.cutoff))
                FDR <- numeric(length(p.cutoff))
                detect.total <- numeric(length(p.cutoff))
                for (j in 1:length(FD)){
                        FD[j] <- sum(t.stat.null>=p.cutoff[j],na.rm=TRUE)
                        detect.total[j] <- sum(t.stat>=p.cutoff[j],na.rm=TRUE)
                        FDR[j] <- FD[j]/detect.total[j]
                }
                #ROC$pos[[i]]$TP <- TP
                ROC$pos[[i]]$FD <- FD
                ROC$pos[[i]]$TPR <- TPR
                ROC$pos[[i]]$FDR <- FDR
                ROC$pos[[i]]$TP.total <- TP.total
                ROC$pos[[i]]$detect.total <- detect.total


        }
        names(ROC$pos) <- names(detection$pos)

	### get ROC for negative strand
	for (i in 1:length(detection$neg)){
                chr <- names(detection$neg)[i]
                genome.start.neg <- detection$genome.start.neg[[chr]]
                if (is.null(genome.start.neg)){
                        cat('Error: can not find ',chr, ' in genome start position\n')
                        return (NULL)
                }
                t.stat <- abs(detection$neg[[i]]$t.stat)
                t.stat.null <- abs(detection.null$neg[[i]]$t.stat)

                ### set cutoff of t.stat and get TP and FP
                cur.idx.methy <- idx.methy$neg[[chr]]
                tmp <- getPvalueCutoff(t.stat, cur.idx.methy, tol)
                p.cutoff <- quantile(tmp$p.cutoff,prob=seq(0.95,0.05,-0.05),na.rm=TRUE)
                TP.total <- tmp$TP.total
                TPR <- seq(0.05,0.95,0.05)
                FD <- numeric(length(p.cutoff))
                FDR <- numeric(length(p.cutoff))
                detect.total <- numeric(length(p.cutoff))
                for (j in 1:length(FD)){
                        FD[j] <- sum(t.stat.null>=p.cutoff[j],na.rm=TRUE)
                        detect.total[j] <- sum(t.stat>=p.cutoff[j],na.rm=TRUE)
                        FDR[j] <- FD[j]/detect.total[j]
                }
                ROC$neg[[i]]$FD <- FD
                ROC$neg[[i]]$TPR <- TPR
                ROC$neg[[i]]$FDR <- FDR
                ROC$neg[[i]]$TP.total <- TP.total
                ROC$neg[[i]]$detect.total <- detect.total


        }
        names(ROC$neg) <- names(detection$neg)
	ROC

}

getROC <- function(detection, idx.methy, tol=5, method='effect.size', idx.sel=NULL)
{	
	if (any(names(detection$pos)!=names(detection$neg))) stop('names of pos and neg of detection should be the same\n')
	ROC <- list()
	for (i in 1:length(detection$pos)){
		chr <- names(detection$pos)[i]
		genome.start.pos <- detection$genome.start.pos[[i]]
		genome.start.neg <- detection$genome.start.neg[[i]]
		if (method=='LR_log'){
			t.stat.pos <- detection$pos[[i]][[method]]
                        t.stat.neg <- detection$neg[[i]][[method]]
		}else{
			t.stat.pos <- abs(detection$pos[[i]][[method]])
			t.stat.neg <- abs(detection$neg[[i]][[method]])
                }
		if (is.null(t.stat.pos) | is.null(t.stat.neg)) stop('incorrect method.\n')
                if (!is.null(idx.sel)) {
			idx.sel$pos[[chr]] <- idx.sel$pos[[chr]][idx.sel$pos[[chr]]>=genome.start.pos]
			idx.sel$neg[[chr]] <- idx.sel$neg[[chr]][idx.sel$neg[[chr]]>=genome.start.neg]
			t.stat.pos[-(idx.sel$pos[[chr]]-genome.start.pos+1)] <- NaN; 
			t.stat.neg[-(idx.sel$neg[[chr]]-genome.start.neg+1)] <- NaN
		}
			
		cur.idx.methy.pos <- idx.methy$pos[[chr]] - genome.start.pos + 1
		cur.idx.methy.neg <- idx.methy$neg[[chr]] - genome.start.neg + 1

                cur.idx.methy.tol.pos <- numeric(0)
                cur.idx.methy.tol.neg <- numeric(0)
		for (j in -tol:tol) { cur.idx.methy.tol.pos<- c(cur.idx.methy.tol.pos, cur.idx.methy.pos+j)  }
                for (j in -tol:tol) { cur.idx.methy.tol.neg<- c(cur.idx.methy.tol.neg, cur.idx.methy.neg+j)  }
		cur.idx.methy.tol.pos <- unique(cur.idx.methy.tol.pos)
                cur.idx.methy.tol.neg <- unique(cur.idx.methy.tol.neg)
	
		tmp <- getPvalueCutoff(t.stat.pos, cur.idx.methy.pos, tol)
                p.cutoff.pos <- rev(tmp$p.cutoff)
                TP.total.pos <- tmp$TP.total

		tmp <- getPvalueCutoff(t.stat.neg, cur.idx.methy.neg, tol)
                p.cutoff.neg <- rev(tmp$p.cutoff)
                TP.total.neg <- tmp$TP.total
		
		p.cutoff <- sort(c(p.cutoff.pos, p.cutoff.neg), decreasing=TRUE)
		TP.total <- TP.total.pos + TP.total.neg
	
		TP <- 1:length(p.cutoff)
                TPR <- TP/TP.total
                FD <- numeric(length(p.cutoff))
                FDR <- numeric(length(p.cutoff))
                detect.total <- numeric(length(p.cutoff))

		for (j in 1:length(FD)){
                        FD[j] <- sum(t.stat.pos[-cur.idx.methy.tol.pos]>=p.cutoff[j],na.rm=TRUE) + 
				sum(t.stat.neg[-cur.idx.methy.tol.neg]>=p.cutoff[j],na.rm=TRUE)
                        detect.total[j] <- sum(t.stat.pos>=p.cutoff[j],na.rm=TRUE) + 
					sum(t.stat.neg>=p.cutoff[j],na.rm=TRUE) 
                        FDR[j] <- FD[j]/detect.total[j]
                }
                ROC[[chr]]$TP <- TP
                ROC[[chr]]$FD <- FD
                ROC[[chr]]$TPR <- TPR
                ROC[[chr]]$FDR <- FDR
                ROC[[chr]]$TP.total <- TP.total
                ROC[[chr]]$detect.total <- detect.total
	}		
	ROC
	
}

getROC.bak <- function(detection, idx.methy, tol=5, method='effect.size', idx.sel=NULL)
{
	### get ROC for positive strand
	ROC <- list()
	for (i in 1:length(detection$pos)){
		chr <- names(detection$pos)[i]
		genome.start.pos <- detection$genome.start.pos[[i]]
		t.stat <- abs(detection$pos[[i]][[method]])
		if (is.null(t.stat)) stop('incorrect method.\n')	
		if (!is.null(idx.sel)) t.stat[idx.sel$pos[[chr]]] <- NaN
		### set cutoff of t.stat and get TP and FP
		cur.idx.methy <- idx.methy$pos[[chr]]
		cur.idx.methy.tol <- numeric(0)
		for (j in -tol:tol) { cur.idx.methy.tol<- c(cur.idx.methy.tol,cur.idx.methy+j)  }
		cur.idx.methy.tol <- unique(cur.idx.methy.tol)
		tmp <- getPvalueCutoff(t.stat, cur.idx.methy, tol)
		p.cutoff <- rev(tmp$p.cutoff)
		TP.total <- tmp$TP.total
		TP <- 1:length(p.cutoff)
		TPR <- TP/TP.total
		FD <- numeric(length(p.cutoff))
		FDR <- numeric(length(p.cutoff))
		detect.total <- numeric(length(p.cutoff))
		for (j in 1:length(FD)){	
			FD[j] <- sum(t.stat[-cur.idx.methy.tol]>=p.cutoff[j],na.rm=TRUE)
			detect.total[j] <- sum(t.stat>=p.cutoff[j],na.rm=TRUE)
			FDR[j] <- FD[j]/detect.total[j]	
		}
		ROC$pos[[i]]$TP <- TP
                ROC$pos[[i]]$FD <- FD
                ROC$pos[[i]]$TPR <- TPR
                ROC$pos[[i]]$FDR <- FDR
                ROC$pos[[i]]$TP.total <- TP.total
                ROC$pos[[i]]$detect.total <- detect.total	
	}
	names(ROC$pos) <- names(detection$pos)

	### get ROC for negative strand
        for (i in 1:length(detection$neg)){
                chr <- names(detection$neg)[i]
                genome.start.neg <- detection$genome.start.neg[[i]]
                t.stat <- abs(detection$neg[[i]][[method]])
                if (is.null(t.stat)) stop('incorrect method.\n')
		if (!is.null(idx.sel)) t.stat[idx.sel$neg[[chr]]] <- NaN
                ### set cutoff of t.stat and get TP and FP
		cur.idx.methy <- idx.methy$neg[[chr]]
		cur.idx.methy.tol <- numeric(0)
                for (j in -tol:tol) { cur.idx.methy.tol<- c(cur.idx.methy.tol,cur.idx.methy+j)  }
                cur.idx.methy.tol <- unique(cur.idx.methy.tol)
		tmp <- getPvalueCutoff(t.stat, cur.idx.methy, tol)
		p.cutoff <- rev(tmp$p.cutoff)
		TP.total <- tmp$TP.total
		TP <- 1:length(p.cutoff)
                TPR <- TP/TP.total
                FD <- numeric(length(p.cutoff))
                FDR <- numeric(length(p.cutoff))
		detect.total <- numeric(length(p.cutoff))
                for (j in 1:length(FD)){
                        FD[j] <- sum(t.stat[-cur.idx.methy.tol]>=p.cutoff[j],na.rm=TRUE)
			detect.total[j] <- sum(t.stat>=p.cutoff[j],na.rm=TRUE)
                        FDR[j] <- FD[j]/detect.total[j]
                }

                ROC$neg[[i]]$TP <- TP
                ROC$neg[[i]]$FD <- FD
                ROC$neg[[i]]$TPR <- TPR
                ROC$neg[[i]]$FDR <- FDR
		ROC$neg[[i]]$TP.total <- TP.total
                ROC$neg[[i]]$detect.total <- detect.total

        }
        names(ROC$neg) <- names(detection$neg)
	
	ROC
}


getPvalueCutoff <- function(t.stat, cur.idx.methy, tol)
{
	p.cutoff <- sapply(cur.idx.methy, function(x,t.stat,tol) if (is.na(t.stat[x])) {return (NaN)} else{return(max(t.stat[(x-tol):(x+tol)],na.rm=TRUE ))} , t.stat=t.stat, tol=tol  )
	result <- list()
	result$p.cutoff <- sort(unique(p.cutoff[!is.na(p.cutoff)]))
	result$TP.total <-  length(result$p.cutoff)
	result	
}


getIdxByMotif <- function(genomeSeq, motif, shift=0)
{
	idx.methy <- list()
	genome.seq.pos <- as.list(genomeSeq$pos)
	genome.seq.neg <- as.list(genomeSeq$neg)
	for (i in 1:length(genome.seq.pos)){
		idx.methy$pos[[i]] <- gregexpr(motif, genome.seq.pos[[i]])[[1]] + shift
	}	
	names(idx.methy$pos) <- names(genome.seq.pos)
	
	for (i in 1:length(genome.seq.neg)){
                idx.methy$neg[[i]] <- gregexpr(reverseSeq(motif), genome.seq.neg[[i]])[[1]] + getLengthOfRegularExpr(motif) - 1 - shift
        }
	names(idx.methy$neg) <- names(genome.seq.neg)	
	idx.methy
}


