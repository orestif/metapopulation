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

source('Metapop SIR adaptivetau.R')



# ==================== Numerical Values ===========================

# Parameters taken from Peel et al (2014), Fig. 2b, with a CCS around 10,000.

par.3 <- list(span=1,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01)

# Explore a range of values for s, tau, mu and N
meta.try <- c(2,4,6,8)
tau.max.try <- c(0,0.25,0.5)
N.try <- c(1000,5000,10000)

par.try <- expand.grid(meta.try,tau.max.try,N.try)
colnames(par.try) <- c('Patches','tau.max','N')

# ===================== RUN SIMULATIONS ==================================

registerDoParallel(cores=15) # Register the parallel backend to be used by foreach

# Run series and return extinction times
sim.series.3 <- sapply(1:nrow(par.try),function(i) {
	par.i <- par.3
	np <- par.try[i,'Patches']
	par.i$s <- as.numeric(rep(par.3['s'],np))
	# varies tau linearly from 0.5 in patch 1 to 0.5-tau.max in the opposite patch
	par.i$tau <- 0.5 - (par.try[i,'tau.max'] * (1-abs(1-(0:(np-1))*(2/np))))
	print(par.i)
	scale.N <- meta.SIR.init(par.i)
	init.i <- c(round(scale.N*par.try[i,'N']-1), 1, rep(0,np-1), rep(0,np))
	names(init.i) <- c(paste('S',1:np,sep=''),paste('I',1:np,sep=''),paste('R',1:np,sep=''))
	
	sims <- meta.SIR.ssa.thin(n.simul=5000, meta.matrix(np,'circle'), par.i, init.i, t.end=20, tau.leap.par=list(epsilon=0.02), dt=0.01)

	save(par.i,init.i,sims,file=paste('~/Documents/Repository/Metapopulation/Series_3/Sims_',i,'.RData',sep=''))
	sapply(sims, function(sim){
		# ADDED NOTE: The code below is incorrect: column numbers should be (np+2):(2*np+1)
		ext <- which(apply(sim[,(np+1):(2*np)],1,sum)==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
})

write.csv(sim.series.3,file='~/Documents/Repository/Metapopulation/Series_3/Extinctions.csv')

save.image('Metapop Series 3.RData')
