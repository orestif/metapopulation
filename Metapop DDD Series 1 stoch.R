####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		WITH ANNUAL BIRTH PULSES
# 	  AND DENSITY-DEPENDENT DEATH RATE
#	-------------------------------------------
#		    SIMULATION SERIES 1
####################################################

# In this series I vary:
# - number of patches in a circular geometry
# - birth pulse phase variation (with a tight birth pulse in each patch)
# - patch size

# Constants:
# - All simulations initiated with 1 case in patch 1
# - span=1, s=100, tau.1=0.5, R0=4, IP=1/12, mu=0.01

# NOTE: Written for run on serversofia

source('Metapop DDD SIR adaptivetau.R')

library(foreach)

registerDoParallel(cores=32) # Register the parallel backend to be used by foreach


# ==================== Numerical Values ===========================

# Parameters taken from Peel et al (2014), Fig. 2b, with a CCS around 10,000.

def.par <- list(life.span=1,d=0,K=1000,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01)

# Explore a range of values for s, tau, mu and K
np.try <- c(2,4,6,8)
tau.max.try <- c(0,0.25,0.5)
K.try <- c(5000,10000)

par.try <- expand.grid(np.try,tau.max.try,K.try)
colnames(par.try) <- c('Patches','tau.max','K')

# ===================== RUN SIMULATIONS ==================================


# Run series and return extinction times
DDD.sim.series.1 <- sapply(1:nrow(par.try),function(i) {
	par.i <- def.par
	np <- par.try$Patches[i]
	par.i$K <- rep(par.try$K[i],np)
	par.i$s <- as.numeric(rep(def.par$s,np))
	# varies tau linearly from 0.5 in patch 1 to 0.5-tau.max in the opposite patch
	par.i$tau <- 0.5 - (par.try$tau.max[i] * (1-abs(1-(0:(np-1))*(2/np))))
	print(par.try[i,])
	init.i <- c(round(meta.DDD.SIR.init(par.i)-1), 1, rep(0,np-1), rep(0,np))
	names(init.i) <- c(paste('S',1:np,sep=''),paste('I',1:np,sep=''),paste('R',1:np,sep=''))
	
	sims <- meta.DDD.SIR.ssa.tau(n.simul=1000, meta.matrix(np,'circle'), par.i, init.i, t.end=20, thin=0.01, tau.leap.par=list(epsilon=0.01))

	save(par.i,init.i,sims,file=paste('DDD_Series_1/Sims_',i,'.RData',sep=''))
	sapply(sims, function(sim){
		ext <- which(apply(sim[,(np+2):(2*np+1)],1,sum)==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
})

write.csv(DDD.sim.series.1,file='DDD_Series_1/Extinctions.csv')

save.image('Metapop DDD Series 1.RData')
