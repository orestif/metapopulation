####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
#		   WITH ANNUAL BIRTH PULSES
# 		AND DENSITY-DEPENDENT DEATH RATE
#	------------------------------------------
#	TWO COMPETING STRAINS WITH CROSS-IMMUNITY
#	-------------------------------------------
#		    SIMULATION SERIES 1
####################################################

# In this series I vary:
# - number of patches in a circular geometry
# - birth pulse phase variation (with a tight birth pulse in each patch)
# - patch size

# Constants:
# - All simulations initiated with 1 case in patch 1
# - life.span=1, s=100, tau.1=0.5, R0=4, IP=1/12, mu=0.01

source('Metapop DDD SIR 2 adaptivetau.R')

library(foreach)

registerDoParallel(cores=15) # Register the parallel backend to be used by foreach


# ==================== Numerical Values ===========================

# Parameters taken from Peel et al (2014), Fig. 2b, with a CCS around 10,000.

def.par <- list(life.span=1,d=0,K=1000,s=100,tau=0.5,R0=c(4,4),IP=c(1/24,1/12),mu=0.01)

# Explore a range of values for s, tau, mu and K
np.try <- c(2,4,8)
s.try <- c(10,100)
K.try <- c(5000,10000)

par.try <- expand.grid(np.try,s.try,K.try)
colnames(par.try) <- c('Patches','s','K')

# ===================== RUN SIMULATIONS ==================================


# Run series and return extinction times
DDD.sim.series.1 <- sapply(1:nrow(par.try),function(i) {
	par.i <- def.par
	np <- par.try$Patches[i]
	par.i$K <- rep(par.try$K[i],np)
	par.i$s <- rep(par.try$s[i],np)
	par.i$tau <- rep(def.par$tau,np)
	print(par.try[i,])
	init.i <- c(round(meta.DDD.SIR.2.init(par.i)), rep(0,3*np))
	init.i[np+1] <- 1 # Ix in patch 1
	init.i[2*np+1+np/2] <- 1 # Iy[] in patch 1+np/2
	names(init.i) <- c(paste('S',1:np,sep=''),paste('Ix',1:np,sep=''),paste('Iy',1:np,sep=''),paste('R',1:np,sep=''))
	
	sims <- meta.DDD.SIR.2.ssa.tau(n.simul=1000, meta.matrix(np,'circle'), par.i, init.i, t.end=20, thin=0.01, tau.leap.par=list(epsilon=0.01))

	save(par.i,init.i,sims,file=paste('DDD_SIR_2_Series_1/Sims_',i,'.RData',sep=''))
	sapply(sims, function(sim){
		ext <- which(apply(sim[,(np+2):(2*np+1)],1,sum)==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
})

write.csv(DDD.sim.series.1,file='DDD_SIR_2_Series_1/Extinctions.csv')

save.image('Metapop DDD SIR 2 Series 1.RData')
