combinePvalues <- function(detection, left.len=1, right.len=6, method = 'L2',is.log=FALSE)
{

	if (length(detection$pos) != length(detection$neg)){
		cat('number of chr should be equal for forward and backward strand.\n')
		return(NULL)
	}
	
	for (i in 1:length(detection$pos)){
		x2.pos <- rep(NaN, length(detection$pos[[i]]$p.value))
		x2.neg <- rep(NaN, length(detection$neg[[i]]$p.value))
		if (is.log==TRUE){
			cur.p.pos <- log10(detection$pos[[i]]$p.value)
			cur.p.neg <- log10(detection$neg[[i]]$p.value)
		}else{
			cur.p.pos <- detection$pos[[i]]$p.value
                        cur.p.neg <- detection$neg[[i]]$p.value
		}
		cur.p.pos[is.na(cur.p.pos)] <- median(cur.p.pos[!is.na(cur.p.pos)])
		cur.p.neg[is.na(cur.p.neg)] <- median(cur.p.neg[!is.na(cur.p.neg)])
		
		### forward strand 
		for (j in 1:length(cur.p.pos)){
			if (j - left.len < 1 | j + right.len > length(cur.p.pos) )
				next
			p.win <- cur.p.pos[(j - left.len):(j + right.len)]
			if (method=='L1'){
				x2.pos[j] <-  sum( p.win[!is.na(p.win)] )
			}else{
				if (method=='L2'){
					x2.pos[j] <- - sum( p.win[!is.na(p.win)]^2 )
				}else{
					stop('method should be \'L1\' or \'L2\'\n')
				}
			}		
		}
		
		### backward strand
		for (j in 1:length(cur.p.neg)){
                        if (j - right.len < 1 | j + left.len > length(cur.p.neg) )
                                next
                        p.win <- cur.p.neg[(j - right.len):(j + left.len)]
			if (method=='L1'){
                                x2.neg[j] <-  sum( p.win[!is.na(p.win)] )
                        }else{
                                if (method=='L2'){
                                        x2.neg[j] <- - sum( p.win[!is.na(p.win)]^2 )
                                }else{
                                        stop('method should be \'L1\' or \'L2\'\n')
                                }
                        }
		}

		detection$pos[[i]]$p.value <- x2.pos
		detection$neg[[i]]$p.value <- x2.neg
	}	
	detection
}



