####################################################
# 	METAPOPULATION MODEL FOR INFECTION DYNAMICS
# 		   WITH ANNUAL BIRTH PULSES
# 	       DENSITY-DEPENDENT DEATH RATE
#         AND MATERNAL TRANSFER OF ANTIBODIES
#	-------------------------------------------
#	ANALYSIS OF RESULTS FROM SIMULATION SERIES 2
####################################################

source('Models/Metapop DDD MSIR adaptivetau.R')
load('Outputs/DDD_MSIR_Stoch_Series_2.RData')

library(dplyr)
library(tidyr)
library(ggplot2)

# --------------- Analyse extinctions ----------------------------

# Record summary statistics
n.sim <- nrow(DDD.MSIR.stoch.series.2)
p.ext <- apply(DDD.MSIR.stoch.series.2,2, function(x) length(which(is.finite(x)))/n.sim)
t.ext.mean <- apply(DDD.MSIR.stoch.series.2,2, mean, na.rm=T)
t.ext.median <- apply(DDD.MSIR.stoch.series.2,2, median, na.rm=T)

range(p.ext)
range(DDD.MSIR.stoch.series.2,na.rm = T)


# Based on an R0 of 4, on average 1/4 of simulations should have an early extinction.
# Exclude the first 250 - Rather crude
p.ext.ex <- apply(DDD.MSIR.stoch.series.2,2, function(x) {
	y <- sort(x,na.last = T)[251:1000]
	length(which(is.finite(y)))/length(y)
	})

# This doesn't seem to be useful
t.ext.mean.ex <- apply(DDD.MSIR.stoch.series.2,2, function(x) {
	y <- sort(x,na.last = T)[251:1000]
	mean(y,na.rm=T)
})


plot(density(DDD.MSIR.stoch.series.2[,1],na.rm=T,kernel = "biweight"))

MSIR.series.ext <- cbind(par.try,p.ext=p.ext,t.ext.mean=t.ext.mean,t.ext.median=t.ext.median,p.ext.ex=p.ext.ex)

MSIR.series.ext <- MSIR.series.ext %>% mutate(MTA=rho*MP)

# Heatmap for P(ext) - MP vs K
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(MP.m),factor(K))) + geom_tile(aes(fill=p.ext)) + facet_grid(IP.m~s, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_gradient2(low=rgb(1,1,1),mid=rgb(1,1,0),high=rgb(1,0,0),midpoint=0.5) + xlab("Duration of MAb protection (months)") + ylab("Population size") + theme(text=element_text(size=18))

# s vs K
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(s),factor(K))) + geom_tile(aes(fill=p.ext)) + facet_grid(IP.m~MP.m, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_continuous(low=rgb(1,1,0.7),high=rgb(1,0,0)) + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18))

MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(s),factor(K))) + geom_tile(aes(fill=p.ext)) + facet_grid(IP.m~MP.m, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_gradient2(low=rgb(1,1,1),mid=rgb(1,1,0),high=rgb(1,0,0),midpoint=0.5) + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18))



# Contour plot
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(log10(1+s),log10(K))) + geom_contour(aes(z=p.ext)) + facet_grid(IP.m~MP.m, labeller=label_both) + labs(fill='P(extinction)') + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18))


# MP vs K - Exclude the first 250 extinctions from each set
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(MP.m),factor(K))) + geom_tile(aes(fill=p.ext.ex)) + facet_grid(IP.m~s, labeller=label_both) + labs(fill='P(extinction)') + scale_fill_continuous(low=rgb(1,1,0.7),high=rgb(0.5,0,0)) + xlab("Duration of MAb protection (months)") + ylab("Population size") + theme(text=element_text(size=18))


# Heatmap for <T(ext)> - MP vs K
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(MP.m),factor(K))) + geom_tile(aes(fill=t.ext.mean)) + facet_grid(IP.m~s, labeller=label_both) + labs(fill='T(extinction)') + scale_fill_continuous(low=rgb(0.9,0,0.5),high=rgb(0,0.9,0.5)) + xlab("Duration of MAb protection (months)") + ylab("Population size") + theme(text=element_text(size=18))

# Bubble plot combining P(ext) and <T(ext)> - s vs K
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(s),factor(K))) + geom_point(aes(col=t.ext.mean,size=p.ext)) + facet_grid(IP.m~MP.m, labeller=label_both) + labs(fill='T(extinction)') + scale_color_continuous(low=rgb(1,1,0.2),high=rgb(0,0,0.8)) + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18)) + scale_size(range = c(2,10))

# MP vs K
MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(MP.m),factor(K))) + geom_point(aes(col=t.ext.mean,size=p.ext)) + facet_grid(IP.m~s, labeller=label_both) + labs(fill='T(extinction)') + scale_color_continuous(low=rgb(1,1,0.2),high=rgb(0,0,0.8)) + xlab("Duration of MAb protection (months)") + ylab("Population size") + theme(text=element_text(size=18)) + scale_size(range = c(3,12))


MSIR.series.ext %>% mutate(IP.m=round(IP*12),MP.m=round(MTA*12)) %>% ggplot(aes(factor(s),factor(K))) + geom_point(aes(col=t.ext.mean,size=p.ext,shape=(p.ext>0.25))) + facet_grid(IP.m~MP.m, labeller=label_both) + labs(fill='T(extinction)') + scale_color_continuous(low=rgb(1,0.6,0.1),high=rgb(0.1,0.1,1)) + xlab("Tightness of birth pulse") + ylab("Population size") + theme(text=element_text(size=18)) + scale_shape_manual(values=c(15,18)) + scale_size_area(max_size = 8)



# -------------- Estimate CCS ---------------


# ================================================ DETERMINISTIC DYNAMICS ===============================

# load('Outputs/DDD_MSIR_Series_2_ODE.RData')

par(mfrow=c(2,2))

sel.1 <- which(ode.par.try$s==0 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[1] & ode.par.try$rho==0)
matplot(ode.DDD.MSIR.series.2[[sel.1]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.1,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

sel.2 <- which(ode.par.try$s==100 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[1] & ode.par.try$rho==0)
matplot(ode.DDD.MSIR.series.2[[sel.2]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.2,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

sel.3 <- which(ode.par.try$s==0 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[3] & ode.par.try$rho==1)
matplot(ode.DDD.MSIR.series.2[[sel.3]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.3,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

sel.4 <- which(ode.par.try$s==100 & ode.par.try$IP==ode.par.list$IP[1] & ode.par.try$MP==ode.par.list$MP[3] & ode.par.try$rho==1)
matplot(ode.DDD.MSIR.series.2[[sel.4]],lwd=2,main=paste(names(ode.par.try),round(as.numeric(ode.par.try[sel.4,]),3),sep="=",collapse=", "),log="y",ylim=c(1,1000),xlim=c(0,10))

