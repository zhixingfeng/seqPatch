enum_seq <- function(motif.len)
{
	motif <- character(4^motif.len)

        for (i in 1:4^motif.len){
                id <- rep(0,motif.len)
                dec.num <- i-1
                for (j in 1:motif.len){
                        id[j] <- dec.num - 4*floor(dec.num/4)
                        dec.num <- floor(dec.num/4)
                }
                motif.tmp <- rep('A', motif.len)
                motif.tmp[id == 1] <- 'C'
                motif.tmp[id == 2] <- 'G'
                motif.tmp[id == 3] <- 'T'
                motif[i] <- paste( motif.tmp, collapse = '')
        }
	motif
}


