hieModel.NC <- function(genomeF.native, genomeF.ref, genomeSeq.native, genomeSeq.ref, min.cvg=5, min.pos=4, left.len=6, right.len=1)
{


}

buildGenomeIndex <- function(genomeSeq, seed.len=8)
{
	genomeSeqIndex <- list()
	
	## build genome index 
	for (i in 1:length(genomeSeq$pos)){
		genomeSeqIndex$pos[[i]] <- .Call('buildGenomeIndexR',genomeSeq$pos[[i]],as.integer(seed.len) )
		genomeSeqIndex$neg[[i]] <- .Call('buildGenomeIndexR',genomeSeq$neg[[i]],as.integer(seed.len) )
	}
	names(genomeSeqIndex$pos) <- names(genomeSeq$pos)
	names(genomeSeqIndex$neg) <- names(genomeSeq$neg)
	genomeSeqIndex
}

buildControl.NC <- function(genomeSeq, genomeSeq.ref, genomeRef.index, genomeF.ref, min.cvg = 3, left.len = 6, right.len = 1 )
{
	
	if (left.len + right.len + 1 != nchar(names(genomeRef.index$pos[[1]])[1])) stop('left.len + right.len + 1 is not equal to context length\n')
	control.data <- list()
	for (i in 1:length(genomeSeq$pos) ){
		control.data$pos[[i]] <- .Call('buildControl_NC', genomeSeq$pos[[i]],genomeSeq.ref$pos, genomeSeq.ref$neg, 
					genomeF.ref$features$ipd_pos, genomeF.ref$features$ipd_neg, genomeF.ref$genome.start.pos,
					 genomeF.ref$genome.start.neg, genomeRef.index$pos, genomeRef.index$neg, as.integer(0),
					 as.integer(min.cvg), as.integer(left.len), as.integer(right.len))

		control.data$neg[[i]] <- .Call('buildControl_NC', genomeSeq$neg[[i]],genomeSeq.ref$pos, genomeSeq.ref$neg,
					genomeF.ref$features$ipd_pos, genomeF.ref$features$ipd_neg, genomeF.ref$genome.start.pos,
					genomeF.ref$genome.start.neg, genomeRef.index$pos, genomeRef.index$neg, as.integer(1),
					as.integer(min.cvg), as.integer(left.len), as.integer(right.len))
	}
        names(control.data$pos) <- names(genomeSeq$pos)
        names(control.data$neg) <- names(genomeSeq$neg)
	control.data
}

