write.z.score.to.CSV <- function(z.score, dir, prefix='csv_z_score')
{
	if (any(names(z.score$pos)!=names(z.score$genome.start.pos)) |
		any(names(z.score$neg)!=names(z.score$genome.start.neg)))
		stop('incorrect reference genome names.')	
	if (!file.exists(dir)) 
		dir.create(dir, recursive=TRUE)	
	
	# postive strand 
	for (i in 1:length(z.score$pos)){
		cur.ref <- names(z.score$pos[i])
		z <- z.score$pos[[i]]
		ref.pos <- z.score$genome.start.pos[[cur.ref]] - 1 + 1:length(z)
		if (substr(dir, nchar(dir), nchar(dir))=='/'){
                        file.name <- paste(dir, prefix,'_',cur.ref,'_pos.csv',sep='')
                }else{
                        file.name <- paste(dir,'/',prefix,'_',cur.ref,'_pos.csv',sep='')
                }
		cat('convert forward strand: ', cur.ref, '\n',sep='')
		.Call('R_API_write_z_score_to_CSV', z, ref.pos, file.name)
	}

	# negative strand 
	for (i in 1:length(z.score$neg)){
                cur.ref <- names(z.score$neg[i])
                z <- z.score$neg[[i]]
                ref.neg <- z.score$genome.start.neg[[cur.ref]] - 1 + 1:length(z)
                if (substr(dir, nchar(dir), nchar(dir))=='/'){
                        file.name <- paste(dir, prefix,'_',cur.ref,'_neg.csv',sep='')
                }else{
                        file.name <- paste(dir,'/',prefix,'_',cur.ref,'_neg.csv',sep='')
                }
                cat('convert backward strand: ', cur.ref, '\n',sep='')
                .Call('R_API_write_z_score_to_CSV', z, ref.neg, file.name)
        }

}

write.genomeFtoCSV <- function(genomeF, dir, prefix='csv_data')
{
	if (any(names(genomeF$features$ipd_pos)!=names(genomeF$genome.start.pos)) | 
		any(names(genomeF$features$ipd_neg)!=names(genomeF$genome.start.neg)) | 
		any(names(genomeF$features$moleculeID_pos)!=names(genomeF$genome.start.pos)) | 
		any(names(genomeF$features$moleculeID_neg)!=names(genomeF$genome.start.neg)))
			stop('names of genomeF are not correct.')
	# positive strand
	cat('convet positive strand\n')
	for (i in 1:length(genomeF$features$ipd_pos)){
		cur.ref <- names(genomeF$features$ipd_pos[i])
		ipd <- genomeF$features$ipd_pos[[i]]
		mol.ID <- genomeF$features$moleculeID_pos[[i]]
		if (length(ipd)!=length(mol.ID)) stop('length of IPD and moleculeID are not the same.')
		ref.pos <- genomeF$genome.start.pos[[i]] - 1 + 1:length(ipd)
		if (substr(dir, nchar(dir), nchar(dir))=='/'){
			file.name <- paste(dir, prefix,'_',cur.ref,'_pos.csv',sep='')
		}else{
			file.name <- paste(dir,'/',prefix,'_',cur.ref,'_pos.csv',sep='')
		}
		cat('convert ', cur.ref, '\n',sep='')
		.Call('R_API_writeGenomeFtoCSV', ipd, mol.ID, ref.pos, file.name)	
	}
	
	# negative strand
	cat('convet negative strand\n')
        for (i in 1:length(genomeF$features$ipd_neg)){
                cur.ref <- names(genomeF$features$ipd_neg[i])
                ipd <- genomeF$features$ipd_neg[[i]]
                mol.ID <- genomeF$features$moleculeID_neg[[i]]
                if (length(ipd)!=length(mol.ID)) stop('length of IPD and moleculeID are not the same.')
                ref.pos <- genomeF$genome.start.neg[[i]] - 1 + 1:length(ipd)
                if (substr(dir, nchar(dir), nchar(dir))=='/'){
                        file.name <- paste(dir, prefix,'_',cur.ref,'_neg.csv',sep='')
                }else{
                        file.name <- paste(dir,'/',prefix,'_',cur.ref,'_neg.csv',sep='')
                }
		cat('convert ', cur.ref, '\n',sep='')
                .Call('R_API_writeGenomeFtoCSV', ipd, mol.ID, ref.pos, file.name)
        }

}

plotLR_log_density <- function(LR, LR.flank, file.name, xlim=NULL)
{
        den.motif <- density(LR$LR_log, na.rm=TRUE)
        den.flank <- density(LR.flank$LR_log, na.rm=TRUE)
        png(file.name, type='cairo')
                if (is.null(xlim)) xlim = range(den.motif$x,den.flank$x)
                plot(den.motif, col='red', xlim=xlim, main='',
                        ylim=range(den.motif$y, den.flank$y), xlab='log likelihood ratio')
                lines(den.flank, col='blue')
        dev.off()
}

map_genome_locus <- function(idx, map, value, val.name)
{
        if ( any(sort(names(idx$pos))!=sort(names(value$pos))) )
                stop('inconsistent names.')
        rl <- list(); rl$pos <- list(); rl$neg <- list()

        # forward strand
        for (i in 1:length(idx$pos) ){
                if (length(idx$pos[[i]]) == 0 ) next
                cur.name <- names(idx$pos[i])
                cur.idx <- idx$pos[[i]]
                cur.value <- value$pos[[i]]
                cur.map <- map[[cur.name]]
                if (is.null(cur.map)) next

                rl$pos[[cur.name]] <- list()
                rl$pos[[cur.name]]$idx <- numeric(0)
                rl$pos[[cur.name]]$hit.name <- character(0)
                rl$pos[[cur.name]]$hit.idx <- numeric(0)
                rl$pos[[cur.name]]$strand <- numeric(0)
                rl$pos[[cur.name]][[val.name]] <- numeric(0)

                if (length(cur.idx)!=length(cur.value)) stop('inconsistent idx and value.')
                for (j in 1:length(cur.idx)){
                        hit.ref <- which(cur.idx[j] == cur.map[,1])
                        if (length(hit.ref)==0) next
                        rl$pos[[cur.name]]$idx <- c(rl$pos[[cur.name]]$idx, cur.idx[j])
                        rl$pos[[cur.name]]$hit.name <- c(rl$pos[[cur.name]]$hit.name, as.character(cur.map[hit.ref,2]) )
                        rl$pos[[cur.name]]$hit.idx <- c(rl$pos[[cur.name]]$hit.idx, as.numeric(cur.map[hit.ref,3]) )
                        rl$pos[[cur.name]]$strand <- c(rl$pos[[cur.name]]$strand, as.character(cur.map[hit.ref,4]))
                        rl$pos[[cur.name]][[val.name]] <- c(rl$pos[[cur.name]][[val.name]], cur.value[j])
                }
		if (length(rl$pos[[cur.name]]$idx) == 0) rl$pos[[cur.name]] <- NULL
        }

	# backward strand
        for (i in 1:length(idx$neg) ){
                if (length(idx$neg[[i]]) == 0 ) next
                cur.name <- names(idx$neg[i])
                cur.idx <- idx$neg[[i]]
                cur.value <- value$neg[[i]]
                cur.map <- map[[cur.name]]
                if (is.null(cur.map)) next

                rl$neg[[cur.name]] <- list()
                rl$neg[[cur.name]]$idx <- numeric(0)
                rl$neg[[cur.name]]$hit.name <- character(0)
                rl$neg[[cur.name]]$hit.idx <- numeric(0)
                rl$neg[[cur.name]]$strand <- numeric(0)
                rl$neg[[cur.name]][[val.name]] <- numeric(0)

                if (length(cur.idx)!=length(cur.value)) stop('inconsistent idx and value.')
                for (j in 1:length(cur.idx)){
                        hit.ref <- which(cur.idx[j] == cur.map[,1])
                        if (length(hit.ref)==0) next
                        rl$neg[[cur.name]]$idx <- c(rl$neg[[cur.name]]$idx, cur.idx[j])
                        rl$neg[[cur.name]]$hit.name <- c(rl$neg[[cur.name]]$hit.name, as.character(cur.map[hit.ref,2]) )
                        rl$neg[[cur.name]]$hit.idx <- c(rl$neg[[cur.name]]$hit.idx, as.numeric(cur.map[hit.ref,3]) )
                        rl$neg[[cur.name]]$strand <- c(rl$neg[[cur.name]]$strand, as.character(cur.map[hit.ref,4]))
                        rl$neg[[cur.name]][[val.name]] <- c(rl$neg[[cur.name]][[val.name]], cur.value[j])
                }
		if (length(rl$neg[[cur.name]]$idx) == 0) rl$neg[[cur.name]] <- NULL

        }

	rl
}




mergeLR_log <- function(detection, idx.list)
{
	LR_log_pos <- numeric(0)
	LR_log_neg <- numeric(0)

	if ( any(names(idx.list$pos)!=names(idx.list$neg)) ) stop('inconsistent chromosome.')

	for (i in 1:length(idx.list$pos)){
        	cur.ref <- names(idx.list$pos[i])

        	# forward
        	idx <- idx.list$pos[[i]] - detection$genome.start.pos[[cur.ref]] + 1
        	idx <- idx[idx>=1]
        	cur.LR <- detection$pos[[cur.ref]]$LR_log[idx]
        	LR_log_pos <- c(LR_log_pos, cur.LR)

        	# backward
        	idx <- idx.list$neg[[i]] - detection$genome.start.neg[[cur.ref]] + 1
        	idx <- idx[idx>=1]
        	cur.LR <- detection$neg[[cur.ref]]$LR_log[idx]
        	LR_log_neg <- c(LR_log_neg, cur.LR)
	}
	rl <- list()
	rl$LR_log_pos <- LR_log_pos
	rl$LR_log_neg <- LR_log_neg
	rl$LR_log <- c(LR_log_pos, LR_log_neg)
	rl 
	
}

getLengthOfRegularExpr <- function(word)
{
	left <- gregexpr('[[]',word)[[1]]
	right <- gregexpr('[]]', word)[[1]]

	if (sum(left!=-1)==0 & sum(right!=-1)==0) return (nchar(word))	
	if (sum(left!=-1)!=sum(right!=-1)) stop('number of [ and ] are not the same')
	if (any(left >= right)) stop('[ should be in left of ]')	
	
	len <- sum(right - left - 2)

	nchar(word) - length(left) - length(right) - len
}

getStatByIdx <- function(detection, idx.sel, stat='LR_log')
{
	# forward strand
	rl <- list()
	rl$pos <- list()
	for (i in 1:length(idx.sel$pos)){
		cur.ref <- names(idx.sel$pos[i])
		cur.idx <- idx.sel$pos[[i]] - detection$genome.start.pos[[cur.ref]] + 1
		cur.idx <- cur.idx[cur.idx>=1]
		rl$pos[[cur.ref]] <- detection$pos[[cur.ref]][[stat]][cur.idx] 
	}
	
	# backward strand 	
	rl$neg <- list()
        for (i in 1:length(idx.sel$neg)){
                cur.ref <- names(idx.sel$neg[i])
		cur.idx <- idx.sel$neg[[i]] - detection$genome.start.neg[[cur.ref]] + 1
                cur.idx <- cur.idx[cur.idx>=1]
                rl$neg[[cur.ref]] <- detection$neg[[cur.ref]][[stat]][cur.idx]
        }
	rl 
}

getRegionByCutoff <- function(detection, cutoff, idx.sel=NULL, stat='LR_log', method='less')
{
	idx <- list()
	idx$pos <- list()
	idx$neg <- list()

	# forward strand
	for (i in 1:length(detection$pos)){
		cur.ref <- names(detection$pos[i])
		if (method=='less'){
			idx$pos[[cur.ref]] <- detection$genome.start.pos[[i]] - 1 + which(detection$pos[[i]][[stat]]<=cutoff)
		}else{ 
			idx$pos[[cur.ref]] <- detection$genome.start.pos[[i]] - 1 + which(detection$pos[[i]][[stat]]>cutoff)
		}
		if (!is.null(idx.sel)){
			 idx$pos[[cur.ref]] <- intersect(idx$pos[[cur.ref]], idx.sel$pos[[cur.ref]] )
		}
	}				
	# backward strand	
	for (i in 1:length(detection$neg)){
                cur.ref <- names(detection$neg[i])
                if (method=='less'){
                        idx$neg[[cur.ref]] <- detection$genome.start.neg[[i]] - 1 + which(detection$neg[[i]][[stat]]<=cutoff)
                }else{
                        idx$neg[[cur.ref]] <- detection$genome.start.neg[[i]] - 1 + which(detection$neg[[i]][[stat]]>cutoff)
                }
                if (!is.null(idx.sel)){
                         idx$neg[[cur.ref]] <- intersect(idx$neg[[cur.ref]], idx.sel$neg[[cur.ref]] )
                }
        }
	idx
}



getR2OfContextEffect <- function(context.effect.mean)
{
	mean.all <- sum(sapply(context.effect.mean, sum)) / sum(sapply(context.effect.mean, length))
	var.all <- sum(sapply(context.effect.mean, function(x,m) sum((x-m)^2), m=mean.all))	
	var.err <- sum(sapply(context.effect.mean, function(x) sum((x-mean(x))^2)))	
	1 - var.err/var.all
}

combine.genomeSeq <- function(genomeSeq.1, genomeSeq.2)
{
	genomeSeq <- genomeSeq.1
	genomeSeq$pos <- c(genomeSeq.1$pos, genomeSeq.2$pos)	
        genomeSeq$neg <- c(genomeSeq.1$neg, genomeSeq.2$neg)
	genomeSeq
}

combine.genomeF <- function(genomeF.1, genomeF.2)
{
	genomeF <- genomeF.1
	genomeF$features$ipd_pos <- c(genomeF.1$features$ipd_pos, genomeF.2$features$ipd_pos)	
        genomeF$features$ipd_neg <- c(genomeF.1$features$ipd_neg, genomeF.2$features$ipd_neg)
	genomeF$genome.start.pos <- c(genomeF.1$genome.start.pos, genomeF.2$genome.start.pos)
        genomeF$genome.start.neg <- c(genomeF.1$genome.start.neg, genomeF.2$genome.start.neg)
	genomeF
}


outlier.rm <- function(x,k=2.25)
{
	tmp <- quantile(x,prob=c(0.25,0.75))
	inter.len <- tmp[2] - tmp[1]
	upper.bd <- tmp[2] + k*inter.len
	lower.bd <- tmp[1] - k*inter.len
	x[x <=upper.bd]
}
bc.transform <- function(x,lambda)
{
	(x^lambda-1)/lambda
}

refresh <- function()
{
	dyn.unload(paste(system.file(package = "seqPatch"), "/libs/seqPatch.so",sep = ""))
	dyn.load(paste(system.file(package = "seqPatch"), "/libs/seqPatch.so",sep = ""))
}
reverseSeq <- function(seq)
{
        .Call('reverseSeq',seq)
}


trim.context.effect <- function(context.effect.bypos, cvg.cutoff = 20, context.ex=NULL)
{
        ce.len <- sapply(context.effect.bypos, function(x) sapply(x,length))
	if (is.null(context.ex))
        	return (mapply( function(x,len,cutoff){ x[len>=cutoff]  } ,context.effect.bypos, ce.len, cutoff = cvg.cutoff ))
	else{
		if (any(names(context.effect.bypos)!=names(context.ex)))
			stop('imcompatible context.effect.bypos and context.ex\n')
		for (i in 1:length(context.effect.bypos)){
			if (i==1000*floor(i/1000))cat('filtered ',i,' contexts\r')
			context.effect.bypos[[i]] <- context.effect.bypos[[i]][ce.len[[i]]>=cvg.cutoff]
			context.ex[[i]] <- context.ex[[i]][ce.len[[i]]>=cvg.cutoff]
		}
		cat('filtered ', length(context.effect.bypos), ' contexts\n')
		return (list(context.effect=context.effect.bypos,context.ex=context.ex))
	}
	
}


filter.outlier <- function(x, k=3.5, tail='right')
{
	if (tail!='right' & tail!='left' & tail!='both') 
		stop('tail should be \'right\', \'left\' or \'both\'')
	.Call('filter_outlier_wrap',as.numeric(x), as.numeric(k), tail)
}

filter.outlier.byGenomeF <- function(ipd)
{
	.Call('filter_outlier_by_genomeF',ipd)
}

logTransGenomeF <- function(genomeF)
{
	for (i in 1:length(genomeF$features$ipd_pos))
                genomeF$features$ipd_pos[[i]] <- sapply(genomeF$features$ipd_pos[[i]], function(x) log(x+0.01) )
        for (i in 1:length(genomeF.native$features$ipd_neg))
                genomeF$features$ipd_neg[[i]] <- sapply(genomeF$features$ipd_neg[[i]], function(x) log(x+0.01) )
	genomeF
}


CollapseContextEffectByPos <- function(context.effect.bypos, cvg.cutoff = 20, is.log=TRUE)
{
	rl <- .Call('collapseContextEffectByPos', context.effect.bypos, as.integer(cvg.cutoff))	
	rl <- as.data.frame(rl)
	names(rl)[1] <- 'ipd.log'
	names(rl)[2:ncol(rl)] <- paste('c',1:(ncol(rl)-1),sep='')
	for (i in 2:ncol(rl)) rl[,i] <- as.factor(rl[,i])
	rl
}

get.subread.R2 <- function(data, is.log = TRUE)
{
	if (is.log==TRUE) for (i in 1:length(data)) data[[i]]$IPD <- log(data[[i]]$IPD+0.01)
	
	### get total variance 
	S <- sapply(data, function(x) sum(x$IPD,na.rm=TRUE))
	len <- sapply(data,nrow)
	S.mean <- sum(S)/sum(len)
	
	T <- sum(sapply(data, function(x,S.mean) sum((x$IPD-S.mean)^2,na.rm=TRUE), S.mean=S.mean))

	### get variance within group
	W <- sum(sapply(data, function(x) sum((x$IPD - mean(x$IPD,na.rm=TRUE))^2 ,na.rm=TRUE)))	
	(T - W)/T	
}


plot.cmp.context.effect <- function(cmp.context.effect, outfile='~/cmp_contexteffect.png')
{
        data <- as.data.frame(cbind(cmp.context.effect$x, cmp.context.effect$y))
        names(data) <- c('x','y')
        model <- lm(y~x,data=data)

        png(file=outfile,type='cairo')
        min.sig<- min(cmp.context.effect$x, cmp.context.effect$y)
        max.sig<- max(cmp.context.effect$x, cmp.context.effect$y)
        plot(cmp.context.effect$x, cmp.context.effect$y, pch=19,cex=0.2, col='blue', xlim=c(min.sig, max.sig), ylim=c(min.sig, max.sig), xlab='context effect E. Coli WGA', ylab='context effect MP WGA ', main=paste('PCC=',cor(cmp.context.effect$x,cmp.context.effect$y),sep=''))
        abline(model$coef[1], model$coef[2],col='red')
        abline(0,1,col='green')
        dev.off()


}




getOverDisp <- function(x)
{
        bd <- quantile(x,prob=c(0.01,0.99) )
        y <- x[x>bd[1] & x<bd[2]]
        return ( sd(y)/mean(abs(y)))

}

getSd.filter <- function(x, is.log=FALSE)
{
	y <- filter.outlier(x, method='high')
	if (is.log==TRUE){
		return(sd(log(y+0.01)))	
	}else{
		return(sd(y))
	}
}

sub.sample.reads <- function(alnsF, alnsIdx , ratio)
{
        idx.sel <- sample(1:length(alnsF), size=floor(length(alnsF)*ratio+0.5))
        alnsF.sel <- alnsF[idx.sel]
        alnsIdx.sel <- alnsIdx[idx.sel, ]
        result <- list()
        result$alnsF <- alnsF.sel
        result$alnsIdx <- alnsIdx.sel
        result
}


