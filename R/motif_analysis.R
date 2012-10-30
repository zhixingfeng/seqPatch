bootstrap.motif.freq.bg <- function(motif, genomeSeq, idx.sel, strand, tol, B, n.context)
{
	motif.len <- nchar(motif[1])
	.Call('bootstrap_bg_motif', motif, genomeSeq, as.integer(idx.sel - 1),
	        as.integer(strand), as.integer(tol), as.integer(motif.len), as.integer(B), as.integer(n.context) )

}
getSigfreq <- function(motif, context.sig)
{
        motif.len <- nchar(motif[1])
        sig.count <- list()
        sig.freq <- list()
        for (i in 1:length(motif)){
                cur.idx <- gregexpr(motif[i], context.sig)
                is.find <- sapply(cur.idx, function(x) all(x!=-1))
                sig.count[[ motif[i] ]] <- sum(is.find)
                sig.freq[[ motif[i] ]] <- sig.count[[ motif[i] ]] / length(context.sig)
                cat(i,'\r')
        }
        cat(length(motif),'\n')
        list(sig.count=sig.count, sig.freq=sig.freq, total = length(context.sig))
}




