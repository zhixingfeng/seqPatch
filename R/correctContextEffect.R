correctContextEffect <- function(alnsF, context.effect, left.len, right.len)
{
	if (left.len + right.len + 1 != nchar(names(context.effect)[1])){
		cat('error: left.len + right.len + 1 should be equal to context length\n')
		return(NULL)	
	}
	assign('alnsF.bak',alnsF)
	alnsF.bak[[1]]$PulseWidth <- 1*alnsF.bak[[1]]$PulseWidth
	tmp <- .Call('correctContextEffect', alnsF.bak, sapply(context.effect, function(x) x), names(context.effect), as.integer(left.len), as.integer(right.len))
	return (alnsF.bak)
}



