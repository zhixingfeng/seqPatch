mixContextEffect.bak <- function(context.effect.bypos.1, context.effect.bypos.2)
{
	if (nchar(names(context.effect.bypos.1)[1]) != nchar(names(context.effect.bypos.2)[1]))
		stop('context length should be the same!\n')
	context.len <- nchar(names(context.effect.bypos.1)[1])
	context.enum <- enum_seq(context.len)

	context.effect.bypos.mix <- lapply(context.enum, function(x,context.effect.1, context.effect.2)
	 			c(context.effect.1[[x]], context.effect.2[[x]]),
				context.effect.1=context.effect.bypos.1, context.effect.2=context.effect.bypos.2)
	names(context.effect.bypos.mix) <- context.enum
	context.effect.bypos.mix
}

mixContextEffect <- function(context.effect.bypos.1, context.effect.bypos.2)
{
        if (nchar(names(context.effect.bypos.1)[1]) != nchar(names(context.effect.bypos.2)[1]))
                stop('context length should be the same!\n')
	.Call('mixContextEffect',context.effect.bypos.1, context.effect.bypos.2)	

}




