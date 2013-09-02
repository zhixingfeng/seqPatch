detectModPropEB <- function(z.score, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	rl <- list()
	rl$pos <- list()
	rl$neg <- list()
	# forward strand
	cat('forward strand:\n')
	for (i in 1:length(z.score$pos)){
		cur.ref <- names(z.score$pos[i])
		cat(cur.ref,':\n', sep='')
		rl$pos[[cur.ref]] <- .Call('R_API_detectModProp_EB', z.score$pos[[i]], as.numeric(mu.0), as.numeric(sigma.0),
			                as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

	}

	# backward strand
	cat('backward strand:\n')
	for (i in 1:length(z.score$neg)){
                cur.ref <- names(z.score$neg[i])
		cat(cur.ref,':\n', sep='')
                rl$neg[[cur.ref]] <- .Call('R_API_detectModProp_EB', z.score$neg[[i]], as.numeric(mu.0), as.numeric(sigma.0),
                                        as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

        }
	
	rl$genome.start.pos <- z.score$genome.start.pos
	rl$genome.start.neg <- z.score$genome.start.neg
	rl

}

EB.mixture.model <- function(z.score, mu.0 = 0, sigma.0 = 1, f1.x, f1.y, max.iter = 100)
{
	.Call('R_API_EBmixture', as.numeric(z.score), as.numeric(mu.0), as.numeric(sigma.0),
		as.numeric(f1.x), as.numeric(f1.y), as.integer(max.iter))

}


