getPvalue.from.t.stat <- function(detection)
{
	p.value <- list()
	p.value$genome.start.pos <- detection$genome.start.pos
	p.value$genome.start.neg <- detection$genome.start.neg
	p.value$pos <- list()
	p.value$neg <- list()

	# forward strand
	for (i in 1:length(detection$pos)){
		cur.ref <- names(detection$pos[i])
		t.stat <- detection$pos[[i]]$t.stat
		mu.native <- detection$pos[[i]]$mu_native
		mu.ctrl <- detection$pos[[i]]$mu_ctrl
		s2.native <- detection$pos[[i]]$sigma2_native
		s2.ctrl <- detection$pos[[i]]$sigma2_ctrl
		N.native <- detection$pos[[i]]$sampleSize_native
		N.ctrl <- detection$pos[[i]]$sampleSize_ctrl
		
		# calculate degree of freedom according to 
		# http://en.wikipedia.org/wiki/Welch's_t_test

		v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
		p.value$pos[[cur.ref]] <- 2*pt(-abs(t.stat),df=v)	
	}

	# backward strand
        for (i in 1:length(detection$neg)){
                cur.ref <- names(detection$neg[i])
                t.stat <- detection$neg[[i]]$t.stat
                mu.native <- detection$neg[[i]]$mu_native
                mu.ctrl <- detection$neg[[i]]$mu_ctrl
                s2.native <- detection$neg[[i]]$sigma2_native
                s2.ctrl <- detection$neg[[i]]$sigma2_ctrl
                N.native <- detection$neg[[i]]$sampleSize_native
                N.ctrl <- detection$neg[[i]]$sampleSize_ctrl

                # calculate degree of freedom according to
                # http://en.wikipedia.org/wiki/Welch's_t_test

                v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
                p.value$neg[[cur.ref]] <- 2*pt(-abs(t.stat),df=v)
        }
	p.value
}

getZvalue.from.t.stat <- function(detection)
{
        z.stat <- list()
        z.stat$genome.start.pos <- detection$genome.start.pos
        z.stat$genome.start.neg <- detection$genome.start.neg
        z.stat$pos <- list()
        z.stat$neg <- list()

        # forward strand
        for (i in 1:length(detection$pos)){
                cur.ref <- names(detection$pos[i])
                t.stat <- detection$pos[[i]]$t.stat
                mu.native <- detection$pos[[i]]$mu_native
                mu.ctrl <- detection$pos[[i]]$mu_ctrl
                s2.native <- detection$pos[[i]]$sigma2_native
                s2.ctrl <- detection$pos[[i]]$sigma2_ctrl
                N.native <- detection$pos[[i]]$sampleSize_native
                N.ctrl <- detection$pos[[i]]$sampleSize_ctrl

                # calculate degree of freedom according to
                # http://en.wikipedia.org/wiki/Welch's_t_test

                v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
                z.stat$pos[[cur.ref]] <- qnorm(pt(t.stat,df=v, log.p=TRUE)  ,log.p=TRUE)
        }

        # backward strand
        for (i in 1:length(detection$neg)){
                cur.ref <- names(detection$neg[i])
                t.stat <- detection$neg[[i]]$t.stat
                mu.native <- detection$neg[[i]]$mu_native
                mu.ctrl <- detection$neg[[i]]$mu_ctrl
                s2.native <- detection$neg[[i]]$sigma2_native
                s2.ctrl <- detection$neg[[i]]$sigma2_ctrl
                N.native <- detection$neg[[i]]$sampleSize_native
                N.ctrl <- detection$neg[[i]]$sampleSize_ctrl

                # calculate degree of freedom according to
                # http://en.wikipedia.org/wiki/Welch's_t_test

                v <- (s2.native/N.native + s2.ctrl/N.ctrl)^2 / (s2.native^2/(N.native^2*(N.native-1)) + s2.ctrl^2/(N.ctrl^2*(N.ctrl-1)))
                z.stat$neg[[cur.ref]] <- qnorm(pt(t.stat,df=v, log.p=TRUE)  ,log.p=TRUE)
        }
        z.stat
}

