getZscoreSingleMolecule.CC <- function(genomeF.native, genomeF.ctrl)
{
	z.score <- list()
	z.score$pos <- list()
	z.score$neg <- list()
	if (any(names(genomeF.native$features$ipd_pos) != names(genomeF.native$features$ipd_neg) ) ){
		stop('incompatible genomeF.native')
	}
		
	# forward strand 
	cat('forward strand\n')
	for (i in 1:length(genomeF.native$features$ipd_pos)){
		cur.ref <- names(genomeF.native$features$ipd_pos[i])

		start.native <- genomeF.native$genome.start.pos[[cur.ref]]
		start.ctrl <- genomeF.ctrl$genome.start.pos[[cur.ref]]

		ipd.native <- genomeF.native$features$ipd_pos[[cur.ref]]
		ipd.ctrl <- genomeF.ctrl$features$ipd_pos[[cur.ref]]

		mol.id.native <- genomeF.native$features$moleculeID_pos[[cur.ref]]	
		mol.id.ctrl <- genomeF.ctrl$features$moleculeID_pos[[cur.ref]]
			
		cat(cur.ref, '\n')
		z.score$pos[[cur.ref]] <- .Call('R_API_getZscoreSingleMolecule_CC', ipd.native, ipd.ctrl,
			 mol.id.native, mol.id.ctrl, as.integer(start.native), as.integer(start.ctrl))
	}		
	
	# backward strand
	cat('backward strand\n')
	for (i in 1:length(genomeF.native$features$ipd_neg)){
                cur.ref <- names(genomeF.native$features$ipd_neg[i])

                start.native <- genomeF.native$genome.start.neg[[cur.ref]]
                start.ctrl <- genomeF.ctrl$genome.start.neg[[cur.ref]]

                ipd.native <- genomeF.native$features$ipd_neg[[cur.ref]]
                ipd.ctrl <- genomeF.ctrl$features$ipd_neg[[cur.ref]]

                mol.id.native <- genomeF.native$features$moleculeID_neg[[cur.ref]]
                mol.id.ctrl <- genomeF.ctrl$features$moleculeID_neg[[cur.ref]]
		
		cat(cur.ref, '\n')
                z.score$neg[[cur.ref]] <- .Call('R_API_getZscoreSingleMolecule_CC', ipd.native, ipd.ctrl,
                         mol.id.native, mol.id.ctrl, as.integer(start.native), as.integer(start.ctrl))
        }
 
	z.score
}

getZscoreSingleMolecule.CC.one.site <- function(x, x.idx, y)
{
	x.group <- split(x,x.idx)
	sapply(x.group, function(x,y) t.test(x,y,var.equal=TRUE)$stat, y=y)

}




