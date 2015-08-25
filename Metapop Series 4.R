####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
#	-------------------------------------------
#		    SIMULATION SERIES
####################################################

# In this series I vary:
# - number of patches in a circular geometry
# - birth pulse phase variation (with a tight birth pulse in each patch)
# - patch size

# Constants:
# - All simulations initiated with 1 case in patch 1
# - span=1, s=100, tau.1=0.5, R0=4, IP=1/12, mu=0.01

# WARNING: tau.leap.par must be reduced for any simulations with high seasonality (0.001 seems a good ballpark)

source('Metapop SIR adaptivetau.R')



# ==================== Numerical Values ===========================

# Parameters taken from Peel et al (2014), Fig. 2b, with a CCS around 10,000.

par.4 <- list(span=1,s=0,tau=0.5,R0=4,IP=1/12,mu=0.01)

# Explore a range of values for s, tau, mu and N
meta.try <- c(2,4)
s.max.try <- c(0,1,10,100)
N.try <- c(1000,5000,10000)

par.try <- expand.grid(meta.try,s.max.try,N.try)
colnames(par.try) <- c('Patches','s.max','N')

save(par.try,file='par_try_series_4.RData')

# ===================== RUN SIMULATIONS ==================================

registerDoParallel(cores=15) # Register the parallel backend to be used by foreach

# Run series and return extinction times
sim.series.4 <- sapply(1:nrow(par.try),function(i) {
	par.i <- par.4
	np <- par.try[i,'Patches']
	par.i$s <- par.try[i,'s.max'] * (1-abs(1-(0:(np-1))*(2/np)))
	# varies s linearly from 0 in patch 1 to s.max in the opposite patch
	par.i$tau <- as.numeric(rep(par.4['tau'],np))
	print(par.try[i,])
	scale.N <- meta.SIR.init(par.i)
	init.i <- c(round(scale.N*par.try[i,'N']-1), 1, rep(0,np-1), rep(0,np))
	names(init.i) <- c(paste('S',1:np,sep=''),paste('I',1:np,sep=''),paste('R',1:np,sep=''))
	
	sims <- meta.SIR.ssa.thin(n.simul=1000, meta.matrix(np,'circle'), par.i, init.i, t.end=20, tau.leap.par=list(epsilon=0.01), dt=0.01)

	save(par.i,init.i,sims,file=paste('~/Documents/Repository/Metapopulation/Series_4/Sims_',i,'.RData',sep=''))
	sapply(sims, function(sim){
		# ADDED NOTE: The code below is incorrect: column numbers should be (np+2):(2*np+1)
		ext <- which(apply(sim[,(np+1):(2*np)],1,sum)==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
})

write.csv(sim.series.4,file='~/Documents/Repository/Metapopulation/Series_4/Extinctions.csv')

save.image('Metapop Series 4.RData')
