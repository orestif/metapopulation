####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
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

# Run on or226 25/11 - 2h40 for 84 series of 200 simulations with 16 cores.

source('Models/Metapop DDD MSIR adaptivetau.R')

library(foreach)
registerDoParallel(cores=16) # Register the parallel backend to be used by foreach


# ==================== Numerical Values ===========================

# Parameters taken from Peel et al (2014), Fig. 2b, with a CCS around 10,000.

def.par <- list(life.span=1,d=0,K=1000,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01,rho=1,MP=0.1)

# Explore a range of values for s, IP, rho and MP
s.try <- c(0,1,10,100)
IP.try <- c(7/365, 1/12, 1/2)
rho.try <- c(0,0.5,1)
MP.try <- c(0.1,0.5,1)

par.try <- expand.grid(s.try, IP.try, rho.try, MP.try)
colnames(par.try) <- c('s', 'IP', 'rho', 'MP')

# Remove duplicate rows with no MTA (rho=0) and variable MP
par.try <- par.try[-which(par.try$rho==0 & par.try$MP>0.1),]


# ===================== RUN SIMULATIONS ==================================


# Run series and return extinction times
DDD.MSIR.stoch.series.1 <- sapply(1:nrow(par.try),function(i) {
	par.i <- replace.par(def.par,par.try[i,])
	np <- 1
	print(par.try[i,])
	init.i <- c(S1=round(meta.DDD.MSIR.init(par.i)-1), I1=1L, R1=0, M1=0)
	sims <- meta.DDD.MSIR.ssa.tau(n.simul=200, meta.matrix(np,'circle'), par.i, init.i, t.end=20, thin=0.01, tau.leap.par=list(epsilon=0.01))

	save(par.i,init.i,sims,file=paste('Outputs/DDD_MSIR_Stoch_Series_1/Sims_',i,'.RData',sep=''))
	ext.sim <- sapply(sims, function(sim){
		ext <- which(sim[,3]==0)
		if(length(ext)==0) NA else sim[ext[1],1]
	})
	print(paste("Extinctions:",length(which(is.finite(ext.sim)))))
	return(ext.sim)
})

write.csv(DDD.MSIR.stoch.series.1,file='Outputs/DDD_MSIR_Stoch_Series_1/Extinctions.csv')

save.image('Outputs/DDD_MSIR_Stoch_Series_1/DDD_MSIR_Stoch_Series_1.RData')
