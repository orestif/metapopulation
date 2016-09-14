####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#	    SERIES 2 - DETERMINISTIC SOLUTIONS
# 		        SINGLE PATCH
####################################################

# RAN ON OR226 04/12/2015

library(foreach)
registerDoParallel(cores=15)

source('Models/Metapop DDD MSIR adaptivetau.R')

ode.def.par <- list(life.span=2,d=0,K=1000,s=100,tau=0.5,R0=4,IP=1/12,rho=1,MP=0.1,mu=0)

# Explore a range of values for s, IP, rho, MP and K
ode.par.list <- list(s=c(0,10,100), IP=c(1/12, 1/6, 1/3), rho=c(0,1), MP=c(1/12,1/6,1/3))
ode.par.try <- do.call(expand.grid,ode.par.list)

# Remove duplicate rows with no MTA (rho=0) and variable MP
ode.par.try <- ode.par.try[-which(ode.par.try$rho==0 & ode.par.try$MP>min(ode.par.list$MP)),]

# ===================== SOLVE ODES ==================================

ode.DDD.MSIR.series.2 <- foreach(i = 1:nrow(ode.par.try)) %dopar% {
	par.i <- replace.par(ode.def.par,ode.par.try[i,])
	np <- 1
	init.i <- c(S1=meta.DDD.MSIR.init(par.i), I1=1, R1=0, M1=0)
	meta.DDD.MSIR.ode(meta.matrix(1,'circle'), par.i, init.i, t.end=20, dt=0.02)
}

save(ode.DDD.MSIR.series.2, ode.par.try,file='Outputs/DDD_MSIR_Series_2_ODE.RData')
