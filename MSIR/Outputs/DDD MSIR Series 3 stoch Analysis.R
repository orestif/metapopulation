####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#	ANALYSIS OF RESULTS FROM SIMULATION SERIES 3
####################################################

# In this series, simulations were started with 5 infected case and R0=4, hence a low probability of initial fade-out

source('Models/Metapop DDD MSIR adaptivetau.R')
load('Outputs/DDD_MSIR_Stoch_Series_3.RData')

library(dplyr)
library(tidyr)
library(ggplot2)

# --------------- Analyse extinctions ----------------------------

# Record summary statistics
n.sim <- nrow(DDD.MSIR.stoch.series.3)
p.ext <- apply(DDD.MSIR.stoch.series.3,2, function(x) length(which(is.finite(x)))/n.sim)
t.ext.mean <- apply(DDD.MSIR.stoch.series.3,2, mean, na.rm=T)
t.ext.median <- apply(DDD.MSIR.stoch.series.3,2, median, na.rm=T)

range(p.ext)
range(DDD.MSIR.stoch.series.3,na.rm = T)

# Density plot of time to extinction from first series
plot(density(DDD.MSIR.stoch.series.3[,1],na.rm=T,kernel = "biweight"))

# Complete Dataframe
MSIR.series.ext <- cbind(par.try,p.ext=p.ext,t.ext.mean=t.ext.mean,t.ext.median=t.ext.median)
MSIR.series.ext <- MSIR.series.ext %>% mutate(MP.m=round(12*rho*MP))

# Heatmap for P(ext) - MP vs K
MSIR.series.ext %>% ggplot(aes(factor(MP.m),factor(K))) + geom_tile(aes(fill=p.ext)) + facet_grid(~s, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_gradient2(low=rgb(1,1,1),mid=rgb(1,1,0),high=rgb(1,0,0),midpoint=0.5) + xlab("Duration of MAb protection (months)") + ylab("Population size") + theme(text=element_text(size=18))

# s vs K
MSIR.series.ext %>% ggplot(aes(factor(s),factor(K))) + geom_tile(aes(fill=p.ext)) + facet_grid(~MP.m, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_gradient2(low=rgb(1,1,1),mid=rgb(1,1,0),high=rgb(1,0,0),midpoint=0.5) + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18))


# Contour plot
MSIR.series.ext %>% ggplot(aes(log10(1+s),log10(K))) + geom_contour(aes(z=p.ext)) + facet_grid(~MP.m, labeller=label_both) + labs(fill='P(extinction)') + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18))



# Heatmap for <T(ext)> - MP vs K
MSIR.series.ext %>% ggplot(aes(factor(MP.m),factor(K))) + geom_tile(aes(fill=t.ext.mean)) + facet_grid(~s, labeller=label_both) + labs(fill='T(extinction)') + scale_fill_continuous(low=rgb(0.9,0,0.5),high=rgb(0,0.9,0.5)) + xlab("Duration of MAb protection (months)") + ylab("Population size") + theme(text=element_text(size=18))

# Bubble plot combining P(ext) and <T(ext)> - s vs K
MSIR.series.ext %>% ggplot(aes(factor(s),factor(K))) + geom_point(aes(col=t.ext.mean,size=p.ext)) + facet_grid(~MP.m, labeller=label_both) + labs(fill='T(extinction)') + scale_color_continuous(low=rgb(1,1,0.2),high=rgb(0,0,0.8)) + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18)) + scale_size(range = c(2,10))


# ----------- Estimate CCS --------------
library(MASS)
# Fit a binomial glm to extinctions ~ log10(K) to each series, and use MASS::dose.p() to estimate the corresponding CCS_50.

MSIR.series.CCS <- MSIR.series.ext %>% group_by(s,MP.m) %>% mutate(N.ext = round(p.ext*n.sim), N.per = round((1-p.ext)*n.sim)) %>% summarise(CCS = as.numeric(10^dose.p(glm(cbind(N.per,N.ext)~log10(K),binomial))))

MSIR.series.CCS %>% ggplot(aes(factor(s),factor(MP.m))) + geom_tile(aes(fill=CCS)) + scale_fill_gradient(high=rgb(1,1,0.5),low=rgb(0.5,0,0))


# ================================================ DETERMINISTIC DYNAMICS ===============================

# load('Outputs/DDD_MSIR_Series_3_ODE.RData')

par(mfrow=c(2,2))

sel.1 <- which(ode.par.try$s==0 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[1] & ode.par.try$rho==0)
matplot(ode.DDD.MSIR.series.3[[sel.1]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.1,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))x

sel.2 <- which(ode.par.try$s==100 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[1] & ode.par.try$rho==0)
matplot(ode.DDD.MSIR.series.3[[sel.2]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.2,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

sel.3 <- which(ode.par.try$s==0 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[3] & ode.par.try$rho==1)
matplot(ode.DDD.MSIR.series.3[[sel.3]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.3,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

sel.4 <- which(ode.par.try$s==100 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[3] & ode.par.try$rho==1)
matplot(ode.DDD.MSIR.series.3[[sel.4]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.4,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

