DLcurve.plot.all <- function (mcmc.list = NULL, sim.dir = NULL, 
					output.dir = file.path(getwd(), "DLcurves"), 
					output.type="png",
					burnin = NULL, verbose = FALSE,  ...) {
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	if(is.null(mcmc.list)) mcmc.list <- get.tfr.mcmc(sim.dir=sim.dir, verbose=verbose, burnin=burnin)
	mcmc.list <- get.mcmc.list(mcmc.list)
	meta <- mcmc.list[[1]]$meta
	dl.countries <- meta$id_DL
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'

    for (country in dl.countries) {
        country.obj <- get.country.object(country, meta, index=TRUE)
        if (verbose) 
            cat("Creating DL graph for", country.obj$name, '(', country.obj$code, ')\n')
        do.call(output.type, list(file.path(output.dir, 
										paste('DLplot_c', country.obj$code, '.', postfix, sep=''))))
        DLcurve.plot(mcmc.list = mcmc.list, country = country.obj$code, 
            burnin = burnin, ...)
        dev.off()
    }
    if (verbose) 
        cat("\nDL plots stored into", output.dir, "\n")
}

stop.if.country.not.DL <- function(country.obj, meta) {
	if (!is.element(country.obj$index, meta$id_DL))
    	stop('Country ', country.obj$name, ' not estimated because no decline observed.')
}

DLcurve.plot <- function (mcmc.list, country, burnin = NULL, pi = 80, tfr.max = 10, 
    nr.curves = NULL, ylim = NULL, xlab = "TFR (reversed)", ylab = "TFR decrement", 
    main = NULL, ...
    ) 
{	
	if(class(mcmc.list) == 'bayesTFR.prediction') {
		if(!is.null(burnin) && burnin != mcmc.list$burnin)
			warning('Prediction was generated with different burnin. Burnin set to ', mcmc.list$burnin)
		burnin <- 0 # because burnin was already cut of the traces
	}
	if(is.null(burnin)) burnin <- 0
    mcmc.list <- get.mcmc.list(mcmc.list)
    meta <- mcmc.list[[1]]$meta
    country <- get.country.object(country, meta)
    stop.if.country.not.DL(country, meta)
    tfr_plot <- seq(0, tfr.max, 0.1)
    dlc <- c()
    U.var <- paste("U_c", country$code, sep = "")
    d.var <- paste("d_c", country$code, sep = "")
    Triangle_c4.var <- paste("Triangle_c4_c", country$code, sep = "")
    gamma.vars <- paste("gamma_", 1:3, "_c", country$code, sep = "")
    nr.curves.from.mc <- if (!is.null(nr.curves)) floor(max(nr.curves, 2000)/length(mcmc.list))
    						else NULL
    for (mcmc in mcmc.list) {
    	th.burnin <- get.thinned.burnin(mcmc,burnin)
    	thincurves.mc <- get.thinning.index(nr.points=nr.curves.from.mc, 
            all.points=mcmc$length - th.burnin)
        traces <- load.tfr.parameter.traces.cs(mcmc, country$code, 
        						burnin=th.burnin, 
								thinning.index=thincurves.mc$index)
        theta <- (traces[, U.var] - traces[, Triangle_c4.var] ) * 
            exp(traces[, gamma.vars])/apply(exp(traces[,gamma.vars]), 1, sum)
        theta <- cbind(theta, traces[, Triangle_c4.var], 
            traces[, d.var])
        dlc <- rbind(dlc, t(apply(theta, 1, DLcurve, tfr = tfr_plot, 
            p1 = mcmc$meta$dl.p1, p2 = mcmc$meta$dl.p2)))
    }
    miny <- min(dlc)
    maxy <- max(dlc)
    decr <- -diff(meta$tfr_matrix[1:meta$T_end, country$index])
    thincurves <- get.thinning.index(nr.curves, dim(dlc)[1])
    ltype <- "l"
    if (thincurves$nr.points == 0) {
        ltype <- "n"
        thincurves$index <- 1
    }
    if (is.null(main)) main <- country$name
    if (is.null(ylim)) ylim <- c(miny, maxy)
    plot(dlc[thincurves$index[1], ] ~ tfr_plot, col = "grey", 
        type = "n", xlim = c(max(tfr_plot), min(tfr_plot)), 
        ylim = ylim, ylab = ylab, xlab = xlab, main = main, ...
        )
    if (thincurves$nr.points > 1) {
        for (i in 2:thincurves$nr.points) {
            lines(dlc[thincurves$index[i], ] ~ tfr_plot, col = "grey")
        }
    }
    dl50 <- apply(dlc, 2, quantile, 0.5)
    lines(dl50 ~ tfr_plot, col = "red", lwd = 2)
    lty <- 2:(length(pi) + 1)
    for (i in 1:length(pi)) {
        al <- (1 - pi[i]/100)/2
        dlpi <- apply(dlc, 2, quantile, c(al, 1 - al))
        lines(dlpi[1, ] ~ tfr_plot, col = "red", lty = lty[i], 
            lwd = 2)
        lines(dlpi[2, ] ~ tfr_plot, col = "red", lty = lty[i], 
            lwd = 2)
    }
    points(decr ~ meta$tfr_matrix[1:(meta$T_end - 
        1), country$index], pch = 19)
    legend("topright", legend = c("median", paste("PI", pi)), 
        lty = c(1, lty), bty = "n", col = "red")
}


tfr.trajectories.table <- function(tfr.pred, country, pi=c(80, 95), half.child.variant=TRUE) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	country <- get.country.object(country, tfr.pred$mcmc.set$meta)
	l <- tfr.pred$nr.projections
	data.matrix <- get.data.imputed(tfr.pred)
	x1 <- as.integer(rownames(data.matrix))
	x2 <- seq(max(x1)+5, by=5, length=l)
	tfr <- as.matrix(data.matrix[,country$index], ncol=1)
	rownames(tfr) <- x1
	pred.table <- matrix(NA, ncol=2*length(pi)+1, nrow=l)
	pred.table[,1] <- get.median.from.prediction(tfr.pred, country$index, country$code)[2:(l+1)]
	colnames(pred.table) <- c('median', rep(NA,ncol(pred.table)-1))
	trajectories <- get.trajectories(tfr.pred, country$code)
	idx <- 2
	for (i in 1:length(pi)) {
		cqp <- get.traj.quantiles(tfr.pred, country$index, country$code, trajectories$trajectories, pi[i])
		al <- (1-pi[i]/100)/2
		if (!is.null(cqp)) {
			pred.table[,idx:(idx+1)] <- t(cqp[,2:(l+1)])
		} else{
			pred.table[,idx:(idx+1)] <- matrix(NA, nrow=l, ncol=2)
		}
		colnames(pred.table)[idx:(idx+1)] <- c(al, 1-al)
		idx <- idx+2
	}
	rownames(pred.table) <- x2
	cn <- colnames(pred.table)[2:ncol(pred.table)]
	pred.table[,2:ncol(pred.table)] <- pred.table[,cn[order(cn)]]
	colnames(pred.table)[2:ncol(pred.table)] <- cn[order(cn)]
	if(half.child.variant) {
		up.low <- get.half.child.variant(median=c(0, pred.table[,1]))
		pred.table <- cbind(pred.table, t(up.low[2:1,2:ncol(up.low)]))
		colnames(pred.table)[(ncol(pred.table)-1):ncol(pred.table)] <- c('-0.5child', '+0.5child')
	}
	return(rbind(cbind(tfr, matrix(NA, nrow=nrow(tfr), ncol=ncol(pred.table)-1)), pred.table))
}

get.trajectories <- function(tfr.pred, country, nr.traj=NULL) {
	traj.file <- file.path(tfr.pred$output.dir, paste('traj_country', country, '.rda', sep=''))
	if (file.exists(traj.file)) {
		load(traj.file)
		thintraj <- get.thinning.index(nr.traj, dim(trajectories)[2]) 
		if (thintraj$nr.points == 0) return(list(trajectories=NULL))
		traj.idx <- thintraj$index
	} else {
		trajectories <- NULL
		traj.idx <- NULL
	}
	if(!is.null(trajectories)) {
		shift <- get.tfr.shift(country, tfr.pred)
	 	if(!is.null(shift)) trajectories <- trajectories + shift
	 	rownames(trajectories) <- get.prediction.years(tfr.pred$mcmc.set$meta, dim(trajectories)[1])
	 }
	return(list(trajectories=trajectories, index=traj.idx))
}

get.median.from.prediction <- function(tfr.pred, country.index, country.code) {
	median <- tfr.pred$quantiles[country.index, '0.5',]
	shift <- get.tfr.shift(country.code, tfr.pred)
	if(!is.null(shift)) median <- median + shift
	return(median)
}
	
get.traj.quantiles <- function(tfr.pred, country.index, country.code, trajectories=NULL, pi=80) {
	al <- (1-pi/100)/2
	quantile.values <- as.numeric(dimnames(tfr.pred$quantiles)[[2]])
	alidx<-round(quantile.values,6)==round(al,6)
	cqp <- NULL
	if (any(alidx)) { # pre-saved quantiles
		alidx2 <- round(quantile.values,6)==round(1-al,6)
		cqp <- rbind(tfr.pred$quantiles[country.index, alidx,], 
							tfr.pred$quantiles[country.index, alidx2,])
	} else { # non-standard quantiles
		reload <- FALSE
		if (is.null(trajectories)) {
			if(tfr.pred$nr.traj > 0) reload <- TRUE
		} else { 
			if (dim(trajectories)[2] < 2000 & tfr.pred$nr.traj > dim(trajectories)[2]) reload <- TRUE
		}
		if(reload) {
			#load 2000 trajectories maximum for computing quantiles
			traj.reload <- get.trajectories(tfr.pred, tfr.pred$mcmc.set$meta$regions$country_code[country.index], 2000)
			trajectories <- traj.reload$trajectories
		}
		if (!is.null(trajectories)) {
			cqp <- apply(trajectories, 1, 
						quantile, c(al, 1-al), na.rm = TRUE)
		}
	}
	shift <- get.tfr.shift(country.code, tfr.pred)
	if(!is.null(shift)) cqp <- cqp + matrix(shift, nrow=nrow(cqp), ncol=ncol(cqp), byrow=TRUE)
	return(cqp)
}
	
tfr.trajectories.plot.all <- function(tfr.pred, 
									output.dir=file.path(getwd(), 'TFRtrajectories'),
									output.type="png", verbose=FALSE, ...) {
	# plots TFR trajectories for all countries
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.countries <- country.names(tfr.pred$mcmc.set$meta)
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	for (country in all.countries) {
		country.obj <- get.country.object(country, tfr.pred$mcmc.set$meta)
		if(verbose)
			cat('Creating TFR graph for', country, '(', country.obj$code, ')\n')

		do.call(output.type, list(file.path(output.dir, 
										paste('TFRplot_c', country.obj$code, '.', postfix, sep=''))))
		tfr.trajectories.plot(tfr.pred, country=country.obj$code, ...)
		dev.off()
	}
	if(verbose)
		cat('\nTrajectory plots stored into', output.dir, '\n')
}

get.half.child.variant <- function(median) {
	l <- length(median)
	increment <- c(0, 0.25, 0.4, 0.5)
	upper <- lower <- c()
	for (i in 1:l) {
		upper <- c(upper, median[i]+increment[min(i,4)])
		lower <- c(lower, median[i]-increment[min(i,4)])
	}
	return(rbind(lower, upper))	
}

tfr.trajectories.plot <- function(tfr.pred, country, pi=c(80, 95), 
								  half.child.variant=TRUE, nr.traj=NULL,
								  xlim=NULL, ylim=NULL, type='b', 
								  xlab='Year', ylab='TFR', main=NULL, ...
								  ) {
	if (missing(country)) {
		stop('Argument "country" must be given.')
	}
	country <- get.country.object(country, tfr.pred$mcmc.set$meta)
	tfr_observed <- tfr.pred$mcmc.set$meta$tfr_matrix_observed
	T_end_c <- tfr.pred$mcmc.set$meta$T_end_c
	tfr_matrix_reconstructed <- get.tfr.reconstructed(tfr.pred$tfr_matrix_reconstructed, tfr.pred$mcmc.set$meta)
	x1 <- as.integer(rownames(tfr_matrix_reconstructed))
	x2 <- as.numeric(dimnames(tfr.pred$quantiles)[[3]])
	#x2 <- seq(max(x1), by=5, length=tfr.pred$nr.projections+1)
	lpart1 <- T_end_c[country$index]
	y1.part1 <- tfr_observed[1:T_end_c[country$index],country$index]
	y1.part2 <- NULL
	lpart2 <- tfr.pred$mcmc.set$meta$T_end - T_end_c[country$index]
	if (lpart2 > 0) 
		y1.part2 <- tfr_matrix_reconstructed[(T_end_c[country$index]+1):nrow(tfr_matrix_reconstructed),country$index]

	trajectories <- get.trajectories(tfr.pred, country$code, nr.traj)
	if(is.null(xlim)) xlim <- c(min(x1,x2), max(x1,x2))
	if(is.null(ylim)) ylim <- c(0, max(trajectories$trajectories, y1.part1, y1.part2, na.rm=TRUE))
	if(is.null(main)) main <- country$name
	# plot historical data: observed
	plot(x1[1:lpart1], y1.part1, type=type, xlim=xlim, ylim=ylim, ylab=ylab, xlab=xlab, main=main, ...
					)
	if(lpart2 > 0) {
		lines(x1[(lpart1+1): length(x1)], y1.part2, pch=2, type='b', col='green')
		lines(x1[lpart1:(lpart1+1)], c(y1.part1[lpart1], y1.part2[1]), col='green') # connection between the two parts
	}
	
	# plot trajectories
	if(!is.null(trajectories$trajectories)) { 
		for (i in 1:length(trajectories$index)) {
			lines(x2, trajectories$trajectories[,trajectories$index[i]], type='l', col='gray')
		}
	}
	# plot median
	tfr.median <- get.median.from.prediction(tfr.pred, country$index, country$code)
	lines(x2, tfr.median, type='l', col='red', lwd=2) 
	# plot given CIs
	lty <- 2:(length(pi)+1)
	for (i in 1:length(pi)) {
		cqp <- get.traj.quantiles(tfr.pred, country$index, country$code, trajectories$trajectories, pi[i])
		if (!is.null(cqp)) {
			lines(x2, cqp[1,], type='l', col='red', lty=lty[i], lwd=2)
			lines(x2, cqp[2,], type='l', col='red', lty=lty[i], lwd=2)
		}
	}
	legend <- c('median', paste('PI', pi))
	col <- rep('red', length(lty)+1)
	if (half.child.variant) {
		lty <- c(lty, max(lty)+1)
		llty <- length(lty)
		up.low <- get.half.child.variant(median=tfr.median)
		lines(x2, up.low[1,], type='l', col='blue', lty=lty[llty])
		lines(x2, up.low[2,], type='l', col='blue', lty=lty[llty])
		legend <- c(legend, '+/- 0.5 child')
		col <- c(col, 'blue')
	}
	legend <- c(legend, 'observed TFR')
	col <- c(col, 'black')
	lty <- c(lty, 1)
	pch <- c(rep(-1, length(legend)-1), 1)
	if(lpart2 > 0) {
		legend <- c(legend, 'imputed TFR')
		col <- c(col, 'green')
		lty <- c(lty, 1)
		pch <- c(pch, 2)
	}
	legend('bottomleft', legend=legend, lty=c(1,lty), bty='n', col=col, pch=pch)
	#abline(h=1, lty=3)
	#abline(h=1.5, lty=3)
	#abline(h=2.1, lty=3)
}

extract.plot.args <- function(...) {
	# split '...' into plot arguments and the rest
	all.plot.args <- names(formals(plot.default))
	args <- list(...)
	which.plot.args <- pmatch(names(args), all.plot.args)
	is.fun.arg <- is.na(which.plot.args)
	return(list(plot.args=args[!is.fun.arg], other.args=args[is.fun.arg]))
}
	
do.plot.tfr.partraces <- function(mcmc.list, func, par.names, main.postfix='', chain.ids=NULL, 
									nr.points=NULL, dev.ncol=5, ...) {
	mcmc.list <- get.mcmc.list(mcmc.list)
	if (is.null(chain.ids)) {
		nr.chains <- length(mcmc.list)
		chain.ids <- 1:nr.chains
	} else {
		nr.chains <- length(chain.ids)
		}
	pars <- list()
	iter <- rep(NA, nr.chains)
	mclen <- rep(0, nr.chains)
	
	# split '...' into function arguments and plot arguments
	split.args <- extract.plot.args(...)
	fun.args <- split.args$other.args
	plot.args <- split.args$plot.args
	thin <- fun.args$thin
	fun.args$thin <- NULL
	if(is.null(fun.args$burnin)) fun.args$burnin <- 0
	orig.burnin <- fun.args$burnin
	i <- 1
	for(chain in chain.ids) {
		mcmc <- mcmc.list[[chain]]
		if (!is.null(thin) || mcmc$thin > 1) {
			consolidated.burn.thin <- burn.and.thin(mcmc, orig.burnin, 
										if (is.null(thin)) mcmc$thin else thin)
			fun.args$burnin <- consolidated.burn.thin$burnin
			if (!is.null(consolidated.burn.thin$index)) fun.args$thinning.index <- consolidated.burn.thin$index
			else {
				if(!is.null(thin)) {
					thin <- max(thin, mcmc$thin)
					fun.args$thinning.index <- unique(round(seq(1,mcmc$length-consolidated.burn.thin$burnin, 
												by=thin/mcmc$thin)))
				} else  {
						fun.args$thinning.index <- seq(1,mcmc$length-consolidated.burn.thin$burnin)
						mclen[i] <- length(fun.args$thinning.index)
					}
			}
			
		}
		pars[[i]] <- eval(do.call(func, c(list(mcmc, par.names=par.names), fun.args)))
		pars[[i]] <- filter.traces(pars[[i]], par.names)
		iter[i] <- mcmc$finished.iter
		if (i==1) {
			par.names.l <- length(colnames(pars[[1]]))
			maxy <- rep(NA, par.names.l)
			miny <- rep(NA, par.names.l)
		}
		ipara<-1
		for (para in colnames(pars[[i]])) {
			maxy[ipara] <- max(maxy[ipara], pars[[i]][,para], na.rm=TRUE)
			miny[ipara] <- min(miny[ipara], pars[[i]][,para], na.rm=TRUE)
			ipara <- ipara+1
		}
		i <- i+1
	}

	if (par.names.l < dev.ncol) {
		ncols <- par.names.l
		nrows <- 1
	} else {
		ncols <- dev.ncol
		nrows <- ceiling(par.names.l/dev.ncol)
	}
	par(mfrow=c(nrows,ncols))
	col <- 1:nr.chains
	maxx<-max(iter)
	if(is.null(plot.args$xlim)) plot.args$xlim <- c(1+orig.burnin,maxx)
	if(is.null(plot.args$xlab)) plot.args$xlab <- 'iterations'
	if(is.null(plot.args$ylab)) plot.args$ylab <- ''
	ylim <- plot.args$ylim
	ipara <- 1
	for (para in colnames(pars[[1]])) {
		#mx <- length(pars[[1]][,para])
		if (is.null(thin)) {
			maxmclen <- max(mclen)
			xindex <- if(maxmclen > 0) seq(1+orig.burnin, maxx, length=maxmclen)
						else (1+orig.burnin):maxx 
		} else xindex <- seq(1+orig.burnin, maxx, by=thin)
		thinpoints <- get.thinning.index(nr.points, length(xindex))
		ltype='l'
		if (thinpoints$nr.points > 0) {
			plot.args$ylim <- if(is.null(ylim)) c(miny[ipara], maxy[ipara]) else ylim
			do.call('plot', c(list(xindex[thinpoints$index], 
						pars[[1]][thinpoints$index, para], 
						main=paste(para, main.postfix),
						col=col[1], type='l'), plot.args))
			if (nr.chains > 1) {
				for (i in 2:nr.chains) {
					lines(xindex[thinpoints$index], 
						pars[[i]][thinpoints$index,para], col=col[i], type='l')
				}
			}
		}
		ipara <- ipara+1
	}
	#stop('')
}

tfr.partraces.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr.parameter.names(trans=TRUE), 
									nr.points=NULL, dev.ncol=5, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	do.plot.tfr.partraces(mcmc.list, 'load.tfr.parameter.traces', chain.ids=chain.ids, 
							nr.points=nr.points, par.names=par.names, dev.ncol=dev.ncol, ...)
}

tfr.partraces.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'),
									chain.ids=NULL, par.names=tfr.parameter.names.cs(trans=TRUE),
									nr.points=NULL, dev.ncol=3, low.memory=TRUE, ...) {

	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	mcmc.list <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.list[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	stop.if.country.not.DL(country.obj, mcmc.list[[1]]$meta)
	do.plot.tfr.partraces(mcmc.list, 'load.tfr.parameter.traces.cs', 
		main.postfix=paste('(',country.obj$name,')', sep=''), chain.ids=chain.ids, nr.points=nr.points, 
		country=country.obj$code, par.names=par.names, dev.ncol=dev.ncol, ...)
}

do.plot.tfr.pardensity <- function(mcmc.list, func, par.names, par.names.ext, main.postfix='', 
								func.args=NULL, chain.ids=NULL, burnin=NULL, dev.ncol=5, ...) {
	if(class(mcmc.list) == 'bayesTFR.prediction') {
		if(!is.null(burnin) && burnin != mcmc.list$burnin)
			warning('Prediction was generated with different burnin. Burnin set to ', mcmc.list$burnin,
					'.\n Use a bayesTFR.mcmc.set object as the first argument, if the original traces should be used.')
		burnin <- 0 # because burnin was already cut of the traces
		if (!is.null(chain.ids) && max(chain.ids) > 1) {
			warning('Thinned traces from all chains used for plotting density.\n For selecting individual chains, use a bayesTFR.mcmc.set object as the first argument.')
			chain.ids <- NULL
		}
	}
	if(is.null(burnin)) burnin <- 0

	mcmc.list <- get.mcmc.list(mcmc.list)
	if (!is.null(chain.ids)) mcmc.list <- mcmc.list[chain.ids]
	par.names.l <- length(par.names.ext)
	if (par.names.l < dev.ncol) {
		ncols <- par.names.l
		nrows <- 1
	} else {
		ncols <- dev.ncol
		nrows <- ceiling(par.names.l/dev.ncol)
	}
	args <- extract.plot.args(...)
	par(mfrow=c(nrows,ncols))
	for (para in par.names) {
		values <- eval(do.call(func, c(list(mcmc.list, par.names=para, burnin=burnin), func.args)))
		values <-  filter.traces(values, par.names)
		for (par.name in colnames(values)) {
			dens <- do.call('density', c(list(values[,par.name]), args$other.args))
			do.call('plot', c(list(dens, main=paste(par.name, main.postfix)), args$plot.args))
		}
	}
}

tfr.pardensity.plot <- function(mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr.parameter.names(trans=TRUE), 
									burnin=NULL, dev.ncol=5, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	do.plot.tfr.pardensity(mcmc.list, 'get.tfr.parameter.traces', chain.ids=chain.ids, par.names=par.names,
							par.names.ext=get.full.par.names(par.names, tfr.parameter.names.extended()),
							burnin=burnin, dev.ncol=dev.ncol, ...)
}

tfr.pardensity.cs.plot <- function(country, mcmc.list=NULL, sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									chain.ids=NULL, par.names=tfr.parameter.names.cs(trans=TRUE), 
									burnin=NULL, dev.ncol=3, low.memory=TRUE, ...) {
	if (is.null(mcmc.list))
		mcmc.list <- get.tfr.mcmc(sim.dir, low.memory=low.memory)
	mcmc.l <- get.mcmc.list(mcmc.list)
	country.obj <- get.country.object(country, mcmc.l[[1]]$meta)
	if (is.null(country.obj$name))
		stop('Country ', country, ' not found.')
	stop.if.country.not.DL(country.obj, mcmc.l[[1]]$meta)
	do.plot.tfr.pardensity(mcmc.list, 'get.tfr.parameter.traces.cs', chain.ids=chain.ids, par.names=par.names,
							par.names.ext=get.full.par.names.cs(par.names, 
											tfr.parameter.names.cs.extended(country.obj$code)),
							main.postfix=paste('(',country.obj$name,')', sep=''),
							func.args=list(country.obj=country.obj),
							burnin=burnin, dev.ncol=dev.ncol, ...)
}

".get.gamma.pars" <- function(pred, ...) UseMethod (".get.gamma.pars")
.get.gamma.pars.bayesTFR.prediction <- function(pred, ...) {
	# estimated by
	# library(MASS)
	# data <- pred$tfr_matrix_reconstructed[12,]
	# gd <- fitdistr(data-min(data)+0.05, densfun='gamma')
	# min(data) is 0.95
	return(list(gamma.pars=list(shape=1.87, rate=0.94), gamma.shift=0.95-0.05, min.value=0.7,
					max.value=NULL))
}

get.tfr.map.parameters <- function(pred, tfr.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, ...) {
	map.pars <- list(pred=pred, quantile=quantile, ...)
	if (same.scale) {
		gp <- .get.gamma.pars(pred)
		data <- pred$quantiles[,as.character(quantile),1]
		q <- if(is.null(tfr.range)) c(min(pmax(data,gp$gamma.shift)), max(data))-gp$gamma.shift
			 else c(max(gp$gamma.shift, tfr.range[1]-gp$gamma.shift), max(gp$gamma.shift, tfr.range[2]-gp$gamma.shift))
		min.max.q <- pgamma(q, shape=gp$gamma.pars[['shape']], rate=gp$gamma.pars[['rate']])
		quantiles <- seq(min.max.q[1], min.max.q[2], length=nr.cats)
		quant.values <- c(gp$min.value, 
				qgamma(quantiles, shape=gp$gamma.pars[['shape']], rate=gp$gamma.pars[['rate']])+gp$gamma.shift)
		if(!is.null(gp$max.value) && gp$max.value > max(quant.values)) quant.values <- c(quant.values, gp$max.value)

		if(!is.null(tfr.range)) {
			if(tfr.range[1] < quant.values[1])
				quant.values <- c(tfr.range[1], quant.values)
			if(tfr.range[1] > quant.values[1])
				quant.values <- quant.values[quant.values >= tfr.range[1]]
			last <- length(quant.values)
			if(tfr.range[2] > quant.values[last])
				quant.values <- c(quant.values, tfr.range[2])
			if(tfr.range[2] < quant.values[last])
				quant.values <- quant.values[quant.values <= tfr.range[2]]
		}
		col <- rainbow(500, start=0, end=0.67)[seq(500, 1, length=length(quant.values)-1)]
		map.pars$catMethod <- quant.values
	} else {
		col <- rainbow(500, start=0, end=0.67)[seq(500, 1, length=nr.cats)]
		map.pars$numCats <- nr.cats
	}
	map.pars$colourPalette <- col
	return(map.pars)		
}

bdem.map.all <- function(pred, output.dir, type='tfr', output.type='png', range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, file.prefix='TFRwrldmap_', ...) {
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.years <- dimnames(pred$quantiles)[[3]]
	output.args <- list()
	postfix <- output.type
	if(output.type=='postscript') postfix <- 'ps'
	filename.arg <- 'filename'
	if(output.type=='postscript' || output.type=='pdf') {filename.arg<-'file'}
	else{output.args[['width']] <- 1000}
	
	map.pars <- get.tfr.map.parameters(pred, tfr.range=range, nr.cats=nr.cats, same.scale=same.scale,
					quantile=quantile, ...)
	
	for (year in all.years) {
		output.args[[filename.arg]] <- file.path(output.dir, 
										paste(file.prefix, year, '.', postfix, sep=''))
		do.call(paste(type, '.map', sep=''), c(list(projection.year=year, device=output.type, 
							device.args=output.args), map.pars))
		dev.off()
	}
	cat('Maps written into', output.dir, '\n')
}

tfr.map.all <- function(pred, output.dir, output.type='png', tfr.range=NULL, nr.cats=50, same.scale=TRUE, 
						quantile=0.5, file.prefix='TFRwrldmap_', ...) {
	bdem.map.all(pred=pred, output.dir=output.dir, type='tfr', output.type=output.type, range=tfr.range,
						nr.cats=nr.cats, same.scale=same.scale, quantile=quantile, file.prefix=file.prefix, ...)
}

".map.main.default" <- function(pred, ...) UseMethod (".map.main.default")

.map.main.default.bayesTFR.prediction <- function(pred, ...) return('TFR: quantile')
	
tfr.map <- function(pred, quantile=0.5, projection.year=NULL, projection.index=1,  device='dev.new', 
						main=NULL, device.args=NULL, ...
				) {
	require(rworldmap)
	years <- dimnames(pred$quantiles)[[3]]
	if(!is.null(projection.year)) {
		projection.index <- (1:length(years))[is.element(years, as.character(projection.year))]
		if(length(projection.index) == 0) stop('Projection year ', projection.year, ' not found.\nAvailable: ', 
												paste(years, collapse=', '))
	}
	if(!is.element(as.character(quantile), dimnames(pred$quantiles)[[2]]))
		stop('Quantile ', quantile, ' not found.\nAvailable: ', paste(dimnames(pred$quantiles)[[2]], collapse=', '))
	tfr <- data.frame(cbind(un=pred$mcmc.set$meta$regions$country_code, 
							tfr=pred$quantiles[, as.character(quantile),projection.index]))
	mtfr <- joinCountryData2Map(tfr, joinCode='UN', nameJoinColumn='un')
	if(is.null(main)) main <- paste(years[projection.index], .map.main.default(pred), quantile)
	if (device != 'dev.cur')
		do.call('mapDevice', c(list(device=device), device.args))
	mapParams<-mapCountryData(mtfr, nameColumnToPlot='tfr', addLegend=FALSE, mapTitle=main, ...
	)
	do.call(addMapLegend, c(mapParams, legendWidth=0.5, legendMar=2, legendLabels='all'))
}