loadCmpH5toGenomeF <- function(cmpH5.file, out.dir, n.chunk=10, normalization.method='bySubread', is.use.CCS=TRUE, is.return=FALSE, mapQV.cutoff=255, reagent='C2')
{
	
	if (!file.exists(out.dir))
                dir.create(out.dir, recursive=TRUE)

	mapQV.cutoff <- mapQV.cutoff - 1

	cmpH5 <- PacBioCmpH5(cmpH5.file)
	alnsIdx <- alnIndex(cmpH5)
	
	### split chunks
	n.subread <- nrow(alnsIdx)
	chunk.size <- floor(n.subread/ n.chunk)
	idx.list <- list()
	for (i in 1:(n.chunk-1)){
		idx.list[[i]] <- ((i-1)*chunk.size+1):(i*chunk.size)
	}
	idx.list[[n.chunk]] <- ((n.chunk-1)*chunk.size+1):n.subread
		
			
	### load data form each chunk
	for (i in 1:n.chunk ){
		alnsF.cur <- getAlignmentsWithFeatures(cmpH5, idx=idx.list[[i]], fxs=list(IPD=getIPD))
		alnsIdx.cur <- alnsIdx[idx.list[[i]], ]
		
		# remove low mapQV subreads
		idx.sel <- which (alnsIdx.cur$mapQV>=mapQV.cutoff)
		alnsF.cur <- alnsF.cur[idx.sel]	
		alnsIdx.cur <- alnsIdx.cur[idx.sel,]

		if (reagent=='C2')
			alnsF.cur <- transformIPD(alnsF.cur, alpha = 0.005, lambda = 0.16)
		if (normalization.method=='bySubread')
			alnsF.cur <- normalizeBySubread(alnsF.cur)
		genomeF.cur <- getFeaturesAlongGenome(alnsF.cur, alnsIdx.cur, is.use.CCS)
		file.out <- paste(out.dir,'/genomeF.chunk.',i,'.Rdata',sep='')
		save(genomeF.cur, file=file.out)
		rm(alnsF.cur, alnsIdx.cur);gc()
		cat('load chunk ',i,'\n')
	}
	cat('load chunk', n.chunk,'\n')
		
	### merge chunks
	print('start to merge chunks')
	rm(genomeF.cur)
	
	genomeF.list <- list()
	for (i in 1:n.chunk){
		file.in <- paste(out.dir,'/genomeF.chunk.',i,'.Rdata',sep='')	
		load (file.in)
		genomeF.list[[i]] <- genomeF.cur
		rm(genomeF.cur)
		cat('load chunk ',i,'\r')
	}
	cat('finished loading ', n.chunk, ' chunks\n')
	genomeF	<- mergeGenomeF(genomeF.list)	
	save(genomeF, file= paste(out.dir,'/genomeF.Rdata',sep=''))
	print('finished')
	rm (genomeF.list);gc()
	if (is.return==TRUE)
		genomeF
}



