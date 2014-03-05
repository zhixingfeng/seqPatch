checkHemiModification <- function(genomeF.native, genomeF.ctrl=NULL, idx.pos, idx.neg, dir=NULL, is.simple=FALSE, format='pdf')
{
	if (length(idx.pos) != length(idx.neg)) 
		stop('length of idx.pos and idx.neg should be the same.')	
	cor.all <- numeric(length(idx.pos))
	for (i in 1:length(idx.pos)){
		ipd.pos <- split(genomeF.native$features$ipd_pos$Streptococcus_pneumoniae_ST556[[idx.pos[i]]],
                	genomeF.native$features$moleculeID_pos$Streptococcus_pneumoniae_ST556[[idx.pos[i]]])
		ipd.neg <- split(genomeF.native$features$ipd_neg$Streptococcus_pneumoniae_ST556[[idx.neg[i]]],
                	genomeF.native$features$moleculeID_neg$Streptococcus_pneumoniae_ST556[[idx.neg[i]]])
		mol.id.pos <- names(ipd.pos)
		mol.id.neg <- names(ipd.neg)
		mol.id.consensus <- intersect(mol.id.pos, mol.id.neg)

		ipd.pos.avg <- sapply(ipd.pos[mol.id.consensus], mean)
		ipd.neg.avg <- sapply(ipd.neg[mol.id.consensus], mean)

		if (length(ipd.pos.avg) < 3 | length(ipd.neg.avg) < 3) next
		cor.all[i] <- cor(ipd.pos.avg, ipd.neg.avg)
		if (is.null(dir))
			next
		if (!file.exists(dir))	
			dir.create(dir, recursive = TRUE)
		if (format=='pdf'){
			file <- paste(dir,'/','hemi_mod_pos_',idx.pos[i],'_neg_',idx.neg[i],'.pdf',sep='')
			pdf(file)
		}else{
			file <- paste(dir,'/','hemi_mod_pos_',idx.pos[i],'_neg_',idx.neg[i],'.png',sep='')
			png(file)
		}
		if (is.simple == TRUE){
			plot(ipd.pos.avg, ipd.neg.avg, xlim=range(ipd.pos.avg, ipd.neg.avg), ylim=range(ipd.pos.avg, ipd.neg.avg),
				xlab='IPD forward strand', ylab='IPD backward strand', main=paste('PCC =', round(cor.all[i],2)))
		}else{
			ipd.ctrl.pos <- genomeF.ctrl$features$ipd_pos$Streptococcus_pneumoniae_ST556[[idx.pos[i]]]
			ipd.ctrl.neg <- genomeF.ctrl$features$ipd_neg$Streptococcus_pneumoniae_ST556[[idx.neg[i]]]
			if (length(ipd.ctrl.pos) < 3 | length(ipd.ctrl.neg) < 3) next
			den.ctrl.pos <- density(ipd.ctrl.pos)
			den.ctrl.neg <- density(ipd.ctrl.neg)
			den.pos <- density(ipd.pos.avg)
			den.neg <- density(ipd.neg.avg)
			#par(mfrow=c(2,2))
			layout(matrix(c(1,0,2,3), 2, 2, byrow = TRUE))
			plot(den.pos, xlim=range(ipd.pos.avg, ipd.neg.avg, ipd.ctrl.pos, ipd.ctrl.neg),
				ylim=range(den.pos$y, den.ctrl.pos$y),
				 main='', xlab='IPD forward strand', col='red')
			lines(den.ctrl.pos, col='blue', lty=22)	
			#plot(1,type='n')
			plot(ipd.pos.avg, ipd.neg.avg, xlim=range(ipd.pos.avg, ipd.neg.avg, ipd.ctrl.pos, ipd.ctrl.neg),
				ylim=range(ipd.pos.avg, ipd.neg.avg, ipd.ctrl.pos, ipd.ctrl.neg),
				xlab='IPD forward strand', ylab='IPD backward strand', main=paste('PCC =', round(cor.all[i],2)))
			plot(den.neg$y, den.neg$x, ylim=range(ipd.pos.avg, ipd.neg.avg, ipd.ctrl.pos, ipd.ctrl.neg),
				xlim=range(den.neg$y, den.ctrl.neg$y),
				type='l', xlab ='Density',ylab='IPD backward strand', col='red')	
			lines(den.ctrl.neg$y, den.ctrl.neg$x, col='blue', lty=22)
		}
		dev.off()	
	}
	if (is.null(dir))
		return(cor.all)
}



plotAvgIPD.density <- function(IPD, mol.id, IPD.ctrl=NULL, mol.id.ctrl=NULL, loci=NULL,
			 ref.genome=NULL, strand=NULL, context.hyperPara=NULL, out.file, tl='density', left.len =6, right.len = 1)
{
	if (length(IPD) != length(mol.id)) 
		stop('inconsistant IPD and mol.id.')
	IPD.list <- split(IPD, mol.id)
	IPD.avg <- sapply(IPD.list, mean)
	dens <- density(IPD.avg)

	if (!is.null(IPD.ctrl) & !is.null(mol.id.ctrl)){
		if (length(IPD.ctrl) != length(mol.id.ctrl))
	                stop('inconsistant IPD.ctrl and mol.id.ctrl.')
		IPD.list.ctrl <- split(IPD.ctrl, mol.id.ctrl)
        	IPD.avg.ctrl <- sapply(IPD.list.ctrl, mean)	
	}
	if (length(IPD.avg)<3 | length(IPD.avg.ctrl)<3) return (FALSE)
	
	if (is.null(loci)){
		if (!is.null(IPD.ctrl) & !is.null(mol.id.ctrl)) {
                        dens.ctrl <- density(IPD.avg.ctrl)
                        png(file=out.file, type='cairo')
                        plot(dens, col='red', xlab='average IPD', xlim=range(dens$x, dens.ctrl$x),
                                ylim = range(dens$y, dens.ctrl$y), main=tl)
                        lines(dens.ctrl, col='green')
                        dev.off()

                }else{
                        png(file=out.file, type='cairo')
                        plot(dens, col='red', xlab='average IPD', xlim=range(dens$x),ylim = range(dens$y),
                                main=tl)
                        dev.off()
                }
	
	}else{	
		if (strand  == 0){
			cur.context <- substr(ref.genome, loci - left.len, loci + right.len)
		}else{
		cur.context <- reverseSeq(substr(ref.genome, loci - right.len, loci + left.len))
		}
		cur.hyperPara <- context.hyperPara[[cur.context]]
		
		x <- seq(cur.hyperPara$theta - 5*sqrt(cur.hyperPara$tau2), 
			cur.hyperPara$theta + 5*sqrt(cur.hyperPara$tau2), length.out=500)
        	y <- dnorm(x, cur.hyperPara$theta, sqrt(cur.hyperPara$tau2) )
	
		if (!is.null(IPD.ctrl) & !is.null(mol.id.ctrl))	{
			dens.ctrl <- density(IPD.avg.ctrl)
			png(file=out.file, type='cairo')
                	plot(dens, col='red', xlab='average IPD', xlim=range(x,dens$x, dens.ctrl$x),
				ylim = range(y, dens$y, dens.ctrl$y), main=tl)
                	lines(dens.ctrl, col='green')
			lines(x,y, col = 'blue', lty=2)
			dev.off()	

		}else{
			png(file=out.file, type='cairo')
			plot(dens, col='red', xlab='average IPD', xlim=range(x,dens$x),ylim = range(y, dens$y),
				main=tl)
			lines(x,y, col = 'blue', lty=2)	
			dev.off()
		}	
	}	
	return (TRUE)
}




