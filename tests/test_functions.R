
test.load.UNtfr <- function(wpp.year=2008) {
	# read the UN TFR input file
	tfr <- bayesTFR:::read.UNtfr(wpp.year)
	stopifnot(length(dim(tfr$data))==2)
	stopifnot(dim(tfr$data)[1] > 150)
	stopifnot(is.element('last.observed', colnames(tfr$data)))
	stopifnot(length(tfr$replaced) == 0)
	stopifnot(length(tfr$added) == 0)
	cat('\n===> Test of loading UN TFR file OK.\n')
}

test.load.UNlocations <- function(wpp.year=2008) {
	tfr <- bayesTFR:::read.UNtfr(wpp.year)
	locs <- bayesTFR:::read.UNlocations(tfr$data, wpp.year)
	stopifnot(length(dim(locs$loc_data)) == 2)
	stopifnot(all(is.element(c('country_code', 'include_code'), colnames(locs$loc_data))))
	stopifnot(dim(locs$loc_data)[1] > 200)
	stopifnot(all(is.element(intersect(c(0,1,2), locs$loc_data[,'include_code']), c(0,1,2))))
	cat('\n===> Test of loading WPP location file OK.\n')
}

test.create.tfr.matrix <- function(wpp.year=2008) {
	tfr <- bayesTFR:::read.UNtfr(wpp.year)
	locs <- bayesTFR:::read.UNlocations(tfr$data, wpp.year)
	tfr.and.regions <- bayesTFR:::get.TFRmatrix.and.regions(tfr$data, locs$loc_data, 
												present.year=2009)
	tfr.matrix <- tfr.and.regions$tfr_matrix
	stopifnot(dim(tfr.matrix)[1] == 12)
	stopifnot(rownames(tfr.matrix)[12] == '2008')
	cat('\n===> Test of creating TFR matrix OK.\n')
}

test.run.mcmc.simulation <- function() {
	sim.dir <- tempfile()
	# run MCMC
	m <- run.tfr.mcmc(iter=5, nr.chains=1, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 5)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 5)
	cat('\n===> Test of running MCMC OK.\n')
	# continue MCMC
	m <- continue.tfr.mcmc(iter=5, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 10)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 10)
	stopifnot(!is.element(900, m$meta$regions$country_code)) # 'World' should not be included
	cat('\n===> Test of continuing MCMC OK.\n')
	# run MCMC for an aggregation
	data.dir <- file.path(.find.package("bayesTFR"), 'data')
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, my.tfr.file=file.path(data.dir, 'my_tfr_template.txt'), burnin=0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	cat('\n===> Test of running MCMC for extra areas OK.\n')
	# run prediction
	pred <- tfr.predict(m, burnin=0, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 10)
	stopifnot(!is.element(903, pred$mcmc.set$regions$country_code))
	npred <- dim(pred$tfr_matrix_reconstructed)[2]
	cat('\n===> Test of running projections OK.\n')
	# run MCMC for another aggregation
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, countries=903, burnin=0)
	# run prediction only for the area 903
	pred <- tfr.predict.extra(sim.dir=sim.dir, verbose=FALSE)

	#stopifnot(is.element(903, pred$mcmc.set$regions$country_code))
	stopifnot(dim(pred$tfr_matrix_reconstructed)[2] == npred+1)
	cat('\n===> Test of running projections on extra areas OK.\n')
	
	projs <- summary(pred, country='Uganda')$projections
	tfr.median.shift(sim.dir, country='Uganda', shift=1.5, from=2051, to=2080)
	shifted.pred <- get.tfr.prediction(sim.dir)
	shifted.projs <- summary(shifted.pred, country='Uganda')$projections
	stopifnot(all(projs[10:15,c(1,3:dim(projs)[2])]+1.5 == shifted.projs[10:15,c(1,3:dim(projs)[2])]))
	stopifnot(all(projs[c(1:9, 16:19),c(1,3:dim(projs)[2])] == shifted.projs[c(1:9, 16:19),c(1,3:dim(projs)[2])]))
	cat('\n===> Test of shifting the median OK.\n')
	shifted.pred <- tfr.median.shift(sim.dir, country='Uganda', reset = TRUE)
	shifted.projs <- summary(shifted.pred, country='Uganda')$projections
	stopifnot(all(projs[,c(1,3:dim(projs)[2])] == shifted.projs[,c(1,3:dim(projs)[2])]))
	cat('\n===> Test of resetting the median OK.\n')
	expert.values <- c(2.3, 2.4, 2.4)
	shift <- expert.values - pred$quantiles[15, '0.5',4:6] # Uganda has index 15
	mod.pred <- tfr.median.set(sim.dir, country='Uganda', values=expert.values, years=2024)
	mod.projs <- summary(mod.pred, country='Uganda')$projections
	stopifnot(all(mod.projs[4:6, c(1,3:dim(projs)[2])]==projs[4:6, c(1,3:dim(projs)[2])]+shift))
	stopifnot(all(mod.projs[c(1:3,7:19), c(1,3:dim(projs)[2])]==projs[c(1:3,7:19), c(1,3:dim(projs)[2])]))
	cat('\n===> Test of setting the median OK.\n')
	unlink(sim.dir, recursive=TRUE)
}

test.thinned.simulation <- function() {
	sim.dir <- tempfile()
	# run MCMC
	m <- run.tfr.mcmc(iter=10, nr.chains=2, output.dir=sim.dir, thin=2)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 10)
	stopifnot(m$mcmc.list[[1]]$length == 6)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 20)
	stopifnot(get.stored.mcmc.length(m$mcmc.list, burnin=4) == 6)
	cat('\n===> Test of running thinned MCMC OK.\n')
	# continue MCMC
	m <- continue.tfr.mcmc(iter=10, output.dir=sim.dir)
	stopifnot(m$mcmc.list[[1]]$finished.iter == 20)
	stopifnot(m$mcmc.list[[1]]$length == 11)
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 40)
	stopifnot(get.stored.mcmc.length(m$mcmc.list, burnin=4) == 16)
	cat('\n===> Test of continuing thinned MCMC OK.\n')
	# run MCMC for an aggregation
	data.dir <- file.path(.find.package("bayesTFR"), 'data')
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, my.tfr.file=file.path(data.dir, 'my_tfr_template.txt'), burnin=0)
	stopifnot(is.element(900, m$meta$regions$country_code)) # 'World' should be included
	cat('\n===> Test of running thinned MCMC for extra areas OK.\n')
	# run prediction
	pred <- tfr.predict(m, burnin=5, verbose=FALSE)
	spred <- summary(pred)
	stopifnot(spred$nr.traj == 16) # 2x8
	stopifnot(pred$mcmc.set$mcmc.list[[1]]$finished.iter == 16)
	stopifnot(length(pred$mcmc.set$mcmc.list) == 1)
	npred <- dim(pred$tfr_matrix_reconstructed)[2]
	cat('\n===> Test of running thinned projections OK.\n')
	m <- run.tfr.mcmc.extra(sim.dir=sim.dir, countries=903, burnin=5)
	# run prediction only for the area 903
	pred <- tfr.predict.extra(sim.dir=sim.dir, verbose=FALSE)
	stopifnot(dim(pred$tfr_matrix_reconstructed)[2] == npred+1)
	cat('\n===> Test of running thinned projections on extra areas OK.\n')
	filename <- tempfile()
	png(filename=filename)
	DLcurve.plot(m, 903, burnin=5, nr.curves=5)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of plotting DL curves with thinned MCMC OK.\n')
	filename <- tempfile()
	png(filename=filename)
	tfr.trajectories.plot(pred, 900, nr.traj=5, pi=c(70,65))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of plotting TFR trajectories with thinned MCMC OK.\n')
	traces <- get.tfr.parameter.traces(m$mcmc.list, burnin=5, 
					thinning.index=c(1, 3, 10, 15))
	stopifnot(nrow(traces)==4)
	cat('\n===> Test of getting parameter traces with thinned MCMC OK.\n')
	filename <- tempfile()
	png(filename=filename)
	tfr.pardensity.plot(m, burnin=10)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of plotting parameter density OK.\n')
	unlink(sim.dir, recursive=TRUE)
}


test.run.mcmc.simulation.auto <- function() {
	sim.dir <- tempfile()
	# run MCMC
	m <- run.tfr.mcmc(iter='auto', output.dir=sim.dir,
					auto.conf=list(iter=10, nr.chains=2, thin=1, burnin=5))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 60)
	cat('\n===> Test of running auto MCMC OK.\n')
	m <- continue.tfr.mcmc(iter='auto', output.dir=sim.dir, auto.conf=list(max.loops=2))
	stopifnot(get.total.iterations(m$mcmc.list, 0) == 100)
	cat('\n===> Test of continuing auto MCMC OK.\n')
	unlink(sim.dir, recursive=TRUE)
}

test.existing.simulation <- function() {
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	m <- get.tfr.mcmc(sim.dir, low.memory=FALSE, burnin=25, chain.ids=c(1,2))
	stopifnot(length(m$mcmc.list)==2)
	stopifnot(dim(m$mcmc.list[[1]]$traces)[1]==25)
	cat('\n===> Test of retrieving MCMC results OK.\n')
}

test.DLcurve <- function() {
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	m <- get.tfr.mcmc(sim.dir)
	filename <- tempfile()
	png(filename=filename)
	DLcurve.plot(m, 'Burkina Faso')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of plotting DL curves OK.\n')
}

test.TFRtrajectories <- function() {
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	tfr.trajectories.plot(pred, 'Australia')
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of plotting TFR trajectories OK.\n')
	t <- tfr.trajectories.table(pred, 'Australia', pi=c(90, 80, 70))
	stopifnot(all(dim(t) == c(30, 9)))
	cat('\n===> Test of tabulating TFR trajectories OK.\n')
}

test.plot.density <- function() {
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	png(filename=filename)
	tfr.pardensity.cs.plot('Ireland', pred)
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of plotting parameter density OK.\n')
}

test.plot.map <- function() {
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	pred <- get.tfr.prediction(sim.dir=sim.dir)
	filename <- tempfile()
	tfr.map(pred, projection.year=2043, device='png', device.args=list(filename=filename))
	dev.off()
	size <- file.info(filename)['size']
	unlink(filename)
	stopifnot(size > 0)
	cat('\n===> Test of creating TFR maps OK.\n')
}
	
test.get.parameter.traces <- function() {
	sim.dir <- file.path(.find.package("bayesTFR"), "ex-data", 'bayesTFR.output')
	m <- get.tfr.mcmc(sim.dir, low.memory=TRUE)
	traces <- get.tfr.parameter.traces(m$mcmc.list, burnin=40, 
					thinning.index=c(4, 25, 29))
	stopifnot(nrow(traces)==3)
	m.check <- get.tfr.mcmc(sim.dir, low.memory=FALSE, burnin=40, chain.ids=c(1,3))
	stopifnot(traces[1,'chi']==m.check$mcmc.list[[1]]$traces[4,'chi'])
	stopifnot(all(traces[c(2,3),'chi']==m.check$mcmc.list[[2]]$traces[c(5,9),'chi']))
	
	traces <- get.tfr.parameter.traces(m$mcmc.list, burnin=40, thin=8)
	stopifnot(nrow(traces)==4)
	stopifnot(traces[2,'psi']==m.check$mcmc.list[[1]]$traces[9,'psi'])
	stopifnot(traces[4,'psi']==m.check$mcmc.list[[2]]$traces[5,'psi'])
	cat('\n===> Test of getting parameter traces OK.\n')
}


	
