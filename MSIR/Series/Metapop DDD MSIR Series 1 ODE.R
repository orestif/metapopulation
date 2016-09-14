####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#	    SERIES 1 - DETERMINISTIC SOLUTIONS
# 		        SINGLE PATCH
####################################################

# RAN ON OR226 24/11

library(foreach)
registerDoParallel(cores=15)

source('Models/Metapop DDD MSIR adaptivetau.R')

def.par <- list(life.span=1,d=0,K=1000,s=100,tau=0.5,R0=4,IP=1/12,mu=0.01,rho=1,MP=0.1)

# Explore a range of values for s, IP, rho and MP
s.try <- c(0,1,10,100)
IP.try <- c(7/365, 1/12, 1/2)
rho.try <- c(0,0.5,1)
MP.try <- c(0.1,0.5,1)

par.try <- expand.grid(s.try, IP.try, rho.try, MP.try)

# Remove duplicate rows with no MTA (rho=0) and variable MP
par.try <- par.try[-which(par.try$rho==0 & MP.try>0.1),]

colnames(par.try) <- c('s', 'IP', 'rho', 'MP')

# ===================== SOLVE ODES ==================================

ode.DDD.MSIR.series.1 <- foreach(i = 1:nrow(par.try)) %dopar% {
	par.i <- replace.par(def.par,par.try[i,])
	np <- 1
	print(par.try[i,])
	init.i <- c(S1=meta.DDD.MSIR.init(par.i), I1=1E-4, R1=0, M1=0)
	meta.DDD.MSIR.ode(meta.matrix(1,'circle'), par.i, init.i, t.end=20, dt=0.02)
}

save(ode.DDD.MSIR.series.1, par.try,file='Outputs/DDD_MSIR_Series_1_ODE.RData')
