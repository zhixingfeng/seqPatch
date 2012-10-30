getGenomeSeq <- function(fastafile)
{
	genomeSeq <- .Call('getGenomeSeq', fastafile)
	genomeSeq$pos <- as.list(genomeSeq$pos)	
        genomeSeq$neg <- as.list(genomeSeq$neg)
	genomeSeq
}	

writeGenomeSeq <- function(seq.obj, outfile)
{
	if (file.exists(outfile)==TRUE) is.rm <- file.remove(outfile)
	rl <- mapply(function(tag,sequence,outfile){cat('>',tag,'\n',sequence,'\n',file=outfile,append=TRUE,sep='')} , 
		 names(seq.obj),seq.obj,MoreArgs=list(outfile=outfile)  )
}


