match.nmol <- function(n.mol, n.mol.null, prop.null, cvg.wga, n.mol.min = 5, cvg.wga.cutoff = 25)
{
	if(!(length(n.mol.null)==length(prop.null) & length(prop.null)==length(cvg.wga))){
		stop('lenght of n.mol.null, prop.null and cvg.wga should be the same.')
	}
	idx <- which(n.mol.null>=n.mol.min)
	n.mol <- as.integer(n.mol[n.mol>=n.mol.min])
	n.mol.null <- as.integer(n.mol.null[idx])
	cvg.wga <- as.integer(cvg.wga[idx])
	prop.null <- prop.null[idx]
	
	idx.wga <- which(cvg.wga>=cvg.wga.cutoff)
	n.mol.null <- n.mol.null[idx.wga]
	prop.null <- prop.null[idx.wga]

	
	n.mol.dens <- hist(n.mol, plot=FALSE, breaks=0.5:(max(n.mol,na.rm=TRUE)+0.5))
	n.mol.null.dens <- hist(n.mol.null, plot=FALSE, breaks=0.5:(max(n.mol.null, na.rm=TRUE)+0.5))
	
	bd <- min(max(n.mol,na.rm=TRUE), max(n.mol.null, na.rm=TRUE))
	#bd <- quantile(n.mol.null, probs=1)	
	#mids <- 1:bd
	
	counts <- c(n.mol.dens$counts[n.mol.dens$mids<=bd], sum(n.mol.dens$counts[n.mol.dens$mids>bd]))
	counts.null <- c(n.mol.null.dens$counts[n.mol.null.dens$mids<=bd], sum(n.mol.null.dens$counts[n.mol.null.dens$mids>bd]))
	dens <- counts / sum(counts)
	dens.null <- counts.null / sum(counts.null)

	p.sample <- dens / dens.null	
	p.sample[dens.null==0] <- 0		
	C <- max(p.sample) + 0.0001
	
	p.sample <- p.sample / C
	
	N <- round(counts.null*p.sample)
	prop.sample <- list()
	n.mol.sample <- list()
	for (i in 1:bd){
		cur.idx <- which(n.mol.null==i)
		if(N[i]<1){
			prop.sample[[i]] <- numeric(0)
			n.mol.sample[[i]] <- numeric(0)
		}else{
			prop.sample[[i]] <- sample(prop.null[cur.idx], N[i]) 
			n.mol.sample[[i]] <- sample(n.mol.null[cur.idx], N[i])	
		}
	}
	
	#cur.idx <- which(n.mol.null>bd)
        #if(N[bd+1]<1){
        #	prop.sample[[bd+1]] <- numeric(0)
        #        n.mol.sample[[bd+1]] <- numeric(0)
        #}else{
       # 	prop.sample[[bd+1]] <- sample(prop.null[cur.idx], N[bd+1])
       #         n.mol.sample[[bd+1]] <- sample(n.mol.null[cur.idx], N[bd+1])
       # }
	
	unlist(prop.sample)
}

getDetectedLoci.prop <- function(detection, prop.cutoff=0.35, n.mol.min=5, cvg.wga.min=20)
{
	detected.loci <- list()
	detected.loci$pos <- list()
	for (i in 1:length(detection$pos)){
		cur.ref <- names(detection$pos[i])
		n.mol <- detection$pos[[i]]$N_0 + detection$pos[[i]]$N_1
		cvg.wga <- detection$pos[[i]]$cvg_wga
		detected.loci$pos[[cur.ref]] <- which(detection$pos[[i]]$prop >= prop.cutoff &
							n.mol >= n.mol.min & cvg.wga >= cvg.wga.min)
	}

	detected.loci$neg <- list()
        for (i in 1:length(detection$neg)){
                cur.ref <- names(detection$neg[i])
                n.mol <- detection$neg[[i]]$N_0 + detection$neg[[i]]$N_1
                cvg.wga <- detection$neg[[i]]$cvg_wga
                detected.loci$neg[[cur.ref]] <- which(detection$neg[[i]]$prop >= prop.cutoff &
                                                        n.mol >= n.mol.min & cvg.wga >= cvg.wga.min)
        }
	
	detected.loci	
}



get.detected.sites <- function(detection, left.cutoff, right.cutoff, stat='z.stat',
			 genomeSeq=NULL, motif=NULL, shift=NULL)
{
	if (!(is.null(genomeSeq)&is.null(motif)&is.null(shift) | !is.null(genomeSeq)&!is.null(motif)&!is.null(shift)))
			stop('genomeSeq, motif and shift should be all NULL or not NULL')
	if (!is.null(motif)){
		if(length(motif)!=length(shift))
			stop('inconsistent motif and shift')
	}

	# get significant genome-wide sites 
	detected.sites <- list()
	detected.sites$pos <- list()
	detected.sites$neg <- list()
	for (i in 1:length(detection$pos)){
		cur.ref <- names(detection$pos[i])
		detected.sites$pos[[cur.ref]] <- which(detection$pos[[i]][[stat]]<=left.cutoff|
						detection$pos[[i]][[stat]]>=right.cutoff) + 
						detection$genome.start.pos[[cur.ref]] - 1

		detected.sites$neg[[cur.ref]] <- which(detection$neg[[i]][[stat]]<=left.cutoff|
                                                detection$neg[[i]][[stat]]>=right.cutoff) +
                                                detection$genome.start.neg[[cur.ref]] - 1
	}
	
	if (!is.null(motif)){
		# get all motif sites
		idx.motif.merge <- list()
		idx.motif.merge$pos <- list()
		idx.motif.merge$neg <- list()
		for (i in 1:length(motif)){
			idx.motif <- getIdxByMotif(genomeSeq, motif[i], shift[i])
			#cat(length(idx.motif$pos[[1]]),' : ', length(idx.motif$neg[[1]]),'\n')
			for (j in 1:length(idx.motif$pos)){
				cur.ref <- names(idx.motif$pos[j])
				idx.motif.merge$pos[[cur.ref]] <- c(idx.motif.merge$pos[[cur.ref]], idx.motif$pos[[cur.ref]])
				idx.motif.merge$neg[[cur.ref]] <- c(idx.motif.merge$neg[[cur.ref]], idx.motif$neg[[cur.ref]])
			}	
		}
		# get overlap of significant sites 
		for (i in 1:length(detected.sites$pos)){
			cur.ref <- names(detected.sites$pos[i])
			detected.sites$pos[[cur.ref]] <- intersect(detected.sites$pos[[cur.ref]],
						idx.motif.merge$pos[[cur.ref]])
			detected.sites$neg[[cur.ref]] <- intersect(detected.sites$neg[[cur.ref]],
                                                idx.motif.merge$neg[[cur.ref]])
		}
	}	
	return(detected.sites)	
}

get.z.stat.from.locfdr <- function(locfdr.output, FDR.cutoff=0.1)
{
	rl <- list()
	idx.left <- which(locfdr.output$mat[,3]<=FDR.cutoff)
	idx.sel <- idx.left[which.max(locfdr.output$mat[idx.left,3])]
	rl$left.x <- locfdr.output$mat[idx.sel,1]
	rl$left.FDR <- locfdr.output$mat[idx.sel,3]
	
	idx.right <- which(locfdr.output$mat[,4]<=FDR.cutoff)
        idx.sel <- idx.right[which.max(locfdr.output$mat[idx.right,4])]
	rl$right.x <- locfdr.output$mat[idx.sel,1]
	rl$right.FDR <- locfdr.output$mat[idx.sel,4]	

	rl
}


getFDRcutoff.genome.wide <- function(detection, file=NULL, nulltype=0, df=20)
{
	# merge z.stat of detection
	z.stat <- numeric(0)
	for (i in 1:length(detection$pos)){
		cur.z.stat <- detection$pos[[i]]$z.stat
		z.stat <- c(z.stat, cur.z.stat)
	}
	for (i in 1:length(detection$neg)){
                cur.z.stat <- detection$neg[[i]]$z.stat
                z.stat <- c(z.stat, cur.z.stat)
        }
	z.stat <- z.stat[!is.na(z.stat)]
	
	# get FDR report 
	if (!is.null(file)){
		pdf(file=file)
		tmp <- locfdr(z.stat, mlests=c(0,1), plot=1, nulltype=1, df=df)
		dev.off()
	}
	rl <- locfdr(z.stat, mlests=c(0,1), plot=0, nulltype=nulltype, df=df)
	rl
}

mergeStatByMotif <- function(stat.motif)
{
	# stat.motif is the output of getStatByMotif
	stat.pool <- numeric(0)
	for (i in 1:length(stat.motif)){
		for (j in 1:length(stat.motif[[i]]$pos)){
			stat.pool <- c(stat.pool, stat.motif[[i]]$pos[[j]])
		}
		for (j in 1:length(stat.motif[[i]]$neg)){
                        stat.pool <- c(stat.pool, stat.motif[[i]]$neg[[j]])
                }
	}
	stat.pool
}

getStatByMotif <- function(detection, genomeSeq, stat, motif, shift)
{
	if (length(motif)!=length(shift)) stop('unmatched motif and shift')
	stat.motif <- list()
	idx.motif <- list()
	for (i in 1:length(motif)){
		idx.motif[[i]] <- getIdxByMotif(genomeSeq, motif[i], shift[i])
		stat.motif[[i]] <- getStatByIdx(detection, idx.motif[[i]], stat=stat)
	}
	names(stat.motif) <- motif
	stat.motif
}

getIdxByMotif.wrap <- function(genomeSeq, motif, shift)
{
	if (length(motif)!=length(shift)) stop('unmatched motif and shift')
	idx.motif <- list()
        for (i in 1:length(motif)){
                idx.motif[[i]] <- getIdxByMotif(genomeSeq, motif[i], shift[i])
        }
        names(idx.motif) <- motif
	idx.motif
}

getPvalue.from.t.stat <- function(detection, is.return.full = FALSE)
{
	p.value <- list()
	p.value$genome.start.pos <- detection$genome.start.pos
	p.value$genome.start.neg <- detection$genome.start.neg
	p.value$pos <- list()
	p.value$neg <- list()

	# forward strand
	for (i in 1:length(detection$pos)){
		cur.ref <- names(detection$pos[i])
		t.stat <- detection$pos[[i]]$t.stat
		mu.native <- detection$pos[[i]]$mu_native
		mu.ctrl <- detection$pos[[i]]$mu_ctrl
		s2.native <- detection$pos[[i]]$sigma2_native
		s2.ctrl <- detection$pos[[i]]$sigma2_ctrl
		N.native <- detection$pos[[i]]$sampleSize_native
		N.ctrl <- detection$pos[[i]]$sampleSize_ctrl
		
		# calculate degree of freedom according to 
		# http://en.wikipedia.org/wiki/Welch's_t_test

		v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
		p.value$pos[[cur.ref]] <- 2*pt(-abs(t.stat),df=v)	
		if (is.return.full==TRUE) 
			detection$pos[[i]]$p.value <- p.value$pos[[cur.ref]]
	}

	# backward strand
        for (i in 1:length(detection$neg)){
                cur.ref <- names(detection$neg[i])
                t.stat <- detection$neg[[i]]$t.stat
                mu.native <- detection$neg[[i]]$mu_native
                mu.ctrl <- detection$neg[[i]]$mu_ctrl
                s2.native <- detection$neg[[i]]$sigma2_native
                s2.ctrl <- detection$neg[[i]]$sigma2_ctrl
                N.native <- detection$neg[[i]]$sampleSize_native
                N.ctrl <- detection$neg[[i]]$sampleSize_ctrl

                # calculate degree of freedom according to
                # http://en.wikipedia.org/wiki/Welch's_t_test

                v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
                p.value$neg[[cur.ref]] <- 2*pt(-abs(t.stat),df=v)
		if (is.return.full==TRUE)
                        detection$neg[[i]]$p.value <- p.value$neg[[cur.ref]]

        }
	if (is.return.full==TRUE)
		return(detection)
	else
		return(p.value)

}

getZvalue.from.t.stat <- function(detection, is.return.full = FALSE)
{
        z.stat <- list()
        z.stat$genome.start.pos <- detection$genome.start.pos
        z.stat$genome.start.neg <- detection$genome.start.neg
        z.stat$pos <- list()
        z.stat$neg <- list()

        # forward strand
        for (i in 1:length(detection$pos)){
                cur.ref <- names(detection$pos[i])
                t.stat <- detection$pos[[i]]$t.stat
                mu.native <- detection$pos[[i]]$mu_native
                mu.ctrl <- detection$pos[[i]]$mu_ctrl
                s2.native <- detection$pos[[i]]$sigma2_native
                s2.ctrl <- detection$pos[[i]]$sigma2_ctrl
                N.native <- detection$pos[[i]]$sampleSize_native
                N.ctrl <- detection$pos[[i]]$sampleSize_ctrl

                # calculate degree of freedom according to
                # http://en.wikipedia.org/wiki/Welch's_t_test

                v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
                z.stat$pos[[cur.ref]] <- qnorm(pt(t.stat,df=v, log.p=TRUE)  ,log.p=TRUE)
        	if (is.return.full == TRUE)
			detection$pos[[i]]$z.stat <- z.stat$pos[[cur.ref]]
	}

        # backward strand
        for (i in 1:length(detection$neg)){
                cur.ref <- names(detection$neg[i])
                t.stat <- detection$neg[[i]]$t.stat
                mu.native <- detection$neg[[i]]$mu_native
                mu.ctrl <- detection$neg[[i]]$mu_ctrl
                s2.native <- detection$neg[[i]]$sigma2_native
                s2.ctrl <- detection$neg[[i]]$sigma2_ctrl
                N.native <- detection$neg[[i]]$sampleSize_native
                N.ctrl <- detection$neg[[i]]$sampleSize_ctrl

                # calculate degree of freedom according to
                # http://en.wikipedia.org/wiki/Welch's_t_test

                v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
                z.stat$neg[[cur.ref]] <- qnorm(pt(t.stat,df=v, log.p=TRUE)  ,log.p=TRUE)
		if (is.return.full == TRUE)
                        detection$neg[[i]]$z.stat <- z.stat$neg[[cur.ref]]

        }
	if (is.return.full == TRUE)
		return(detection)
	else
        	return(z.stat)
}

simOverlap <- function(total.1, sampled.1, total.2, sampled.2, B=1000)
{
	x <- numeric(B)
	for (i in 1:B){
		idx.1 <- sample(1:total.1, size = sampled.1)	
		idx.2 <- sample(1:total.2, size = sampled.2)
		x[i] <- length(intersect(idx.1, idx.2))
	}
	x
}

