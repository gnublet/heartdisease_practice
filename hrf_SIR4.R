#######AUTHOR: KEVIN WONG  3/29/2016
#TODO: jittering, TR = 2, various presentation times

#######FUNCTIONS#####################
#auxiliary functions#########
#pure hrf
hrf=function(t, parss, onset=0){
  a1=parss$a1
  a2=parss$a2
  b1=parss$b1
  b2=parss$b2
  g=parss$g
  max.bold=parss$max.bold
  
  d1=a1*b1
  d2=a2*b2
  
  grid=t
  new.t=t-onset
  pos.t.idx=new.t>=0#nonnegative indices of tau (scaled time)
  pos.t=new.t[pos.t.idx==TRUE]#non-negative values of tau (scaled time)
  y.hrf=1/max.bold  * (((pos.t/d1)^a1)*exp(-(pos.t-d1)/b1) - g*(pos.t/d2)^a2*exp(-(pos.t-d2)/b2))
  grid[pos.t.idx==TRUE]=y.hrf
  grid[pos.t.idx==FALSE]=0
  return(grid)
}

#noisy hrf with noise parss$sigma?
rhrf_conv=function(t,parss,times,plot=FALSE){
  n.times=length(times)
  hrfs=array(NA,c(n.times,length(t)))#hrfs ~ all.hrfs?
  if(n.times==1)all.hrfs=hrf(t=t,parss=parss,times)
  if(n.times>1)all.hrfs=sapply(times,hrf,t=t,parss=parss)
  ms=all.hrfs%*%matrix(1,n.times,1)
  if(plot==TRUE){
    plot(t,ms,lwd=5,type="l")
    matlines(t,all.hrfs,type="l",lty=1)
    abline(v=times,lty=3)
  }
  rnorm(length(t), ms, .005)
}

#convolution
hrf_conv=function(t,parss,times,plot=FALSE){
  n.times=length(times)
  hrfs=array(NA,c(n.times,length(t)))
  if(n.times==1)all.hrfs=hrf(t=t,parss=parss,times)
  if(n.times>1)all.hrfs=sapply(times,hrf,t=t,parss=parss)
  ms=all.hrfs%*%matrix(1,n.times,1) #matrix multiplication
  if(plot==TRUE){
    plot(t,ms,lwd=5,type="l")
    matlines(t,all.hrfs,type="l",lty=1)
    abline(v=times,lty=3)
  }
  return(ms)
}

#simulates the true bold response
hrf_sim_truth = function(truth, tmax, tr, ts, tstep, pres_ts, acq_ts, n_acqs,parss){
  T=matrix(acq_ts,n_trials,n_acqs,byrow=TRUE)#trials
  t = T[1,]
  Y=array(NA,c(n_trials,n_acqs))
  ys=sapply(ts,hrf_conv,parss=parss,times=pres_ts)
  return(ys)
}

#simulates the noisy measurements
hrf_sim_measurements = function(truth, tmax, tr, ts, tstep, pres_ts, acq_ts, n_acqs,parss){
  n.trials=1 # number of trials
  times=matrix(acq_ts,n.trials,n_acqs,byrow=TRUE)#trials
  t = times[1,]
  Y=array(NA,c(n.trials,n_acqs))
  y = rep(NA, c(n_acqs))
  
  # this will generate random BOLD measures for the experiment
  #for(n in 1:n.trials){
  #  Y[n]=rhrf_conv(T[n,],parss,n_acqs)
  #}
  y = rhrf_conv(times[1,],parss,pres_ts)
  return(y)
}

########Sample Importance Resampling for estimation of parameters of digamma hrf##########
hrf_sir = function(N, y, truth = NULL, parss, prior, ts, acq_ts, pres_ts, sys_noise, PLOT = F){
  ######initialize####### 
  n=length(y) ##number of observations
  d=6         ##dimension
  
  ##(1)simulate particles
  #particles=matrix(0,nrow=N,ncol=d) ## particles: X_t^{(i)}
  particles = array(0, dim = c(N, d))
  #particles_history = matrix(0, nrow = N, ncol = n)
  particles_history = array(0, dim = c(N, d, n))
  #weight_history = matrix(NA, nrow = N, ncol = n)
  weight_history = array(0, dim = c(N, n))
  for(j in 1:d){
    particles[,j] = rnorm(N, prior[2*j-1], prior[2*j])
  } 
  particles_history[ , ,1] = particles[,]
  ##(2) calculate weights
  w = rep(NA, N)
  for(i in 1: N){ 
    temp_par = parss
    #temp_par$max.bold = particles[j,1]
    temp_par$max.bold = particles[i,1]
    temp_par$a1 = particles[i,2]
    temp_par$a2 = particles[i,3]
    temp_par$b1 = particles[i,4]
    temp_par$b2 = particles[i,5]
    temp_par$g = particles[i,6]
    w[i] = dnorm(y[1], hrf_conv(acq_ts, temp_par, pres_ts)[1,1], parss$sigma_meas)
    #if(w[j]<.00001) w[j]=.00001
  }
  
  #w = dnorm(y[1], hrf_conv(acq_ts[1], temp_par, pres_ts), par$sigma_meas )#hrf.conv=function(t,parss,times,plot=FALSE)
  
  ##(3) log of estimate of likelihood
  logl=log(mean(w))
  
  ##OPTIONAL -- STORAGE OF FEATURES OF FILTER -- MEAN AND VAR
  #M.st=matrix(0,nrow=n,ncol=d)
  #V.st=matrix(0,nrow=n,ncol=d)
  w=w/sum(w)
  weight_history[,1] = w
  #for(j in 1:d){
  #  M.st[1,]=sum(w*particles[,j])
  #  V.st[1,]=sum(w*particles[,j]^2)-M.st[1,]^2
  #}#deleted plot options
  
  ######rejuvenation######
  #Augmentation (resample)
  for(k in 2:n){##LOOP through acquisition times
    ##(1) Resample 
    index=sample(1:N,prob=w,size=N,rep=T)#sample with replacement (this leads to degeneracy)
    #hist(index, 100)
    ##particles=particles[index,]
    #particles = matrix(particles[index], ncol = 1)
    for(j in 1:d){particles[,j] = particles[index,j] } 
    
    ##(2) Propagate particles
    #fixed particle (probability of heads) since we have a static parameter
    for(j in 1:d){particles[,j]=particles[,j]+rnorm(N,0,sys_noise[j])}#rnorm(N,0,parss$sigma_sys)} 
    
    #particles[,1] = rbeta(N, prior[1], prior[2])
    particles_history[,,k] = particles[,] 
    #try to jitter particles
    #temp_particles = matrix(0,nrow=N,ncol=d)
    #temp_particles[,1] = abs(particles[,1] + rnorm(N, 0, .01) )#rbeta(N, par[1], par[2])
    ##for(k in 1:N)
    #{
    #  if (temp_particles[k,1]>1)
    #  {temp_particles[k,1] = 1 - ((temp_particles[k,1])-1)}
    #}
    #particles[,1] = temp_particles[,1]
    
    ##(3) Weight

    for(i in 1: N){ 
      temp_par = parss
      #temp_par$max.bold = particles[j,1]
      temp_par$max.bold = particles[i,1]
      temp_par$a1 = particles[i,2]
      temp_par$a2 = particles[i,3]
      temp_par$b1 = particles[i,4]
      temp_par$b2 = particles[i,5]
      temp_par$g = particles[i,6]
      w[i] = dnorm(y[k], (hrf_conv(acq_ts, temp_par, pres_ts))[k,1], parss$sigma_meas)
      #if(w[j]<.00001) w[j]=.00001#make sure weights are nonzero
    }
    #w = dnorm(y[i], hrf_conv(acq_ts[i], temp_par, pres_ts), par$sigma_meas )
    ##(4) update log of estimate of likelihood
    logl=logl+log(mean(w))
    
    ##OPTIONAL -- STORE FEATURES OF THE PARTICLE APPROXIMATION
    #normalise weights
    w=w/sum(w)
    weight_history[,k] = w
    #for(j in 1:d){
    #  M.st[i,]=sum(w*particles[,j])
    #  V.st[i,]=sum(w*particles[,j]^2)-M.st[i,]^2
    #}
    #plot option
    #if(PLOT){#kw: edited for debug
    #  par(mfrow=c(1,3))
    #  plot(range(particles[,1]),c(0,max(w)),xlab="X_t",ylab="Weight",type="n")
    #  for(j in 1:N) lines(c(particles[j,1],particles[j,1]),c(0,w[j]))
    #  if(length(truth>=1)) abline(v=truth[i],col=2)
    #  plot(density(particles[,1],weights=w),xlab="X_t",ylab="Density",main="")
    #  hist(particles)
    #  if(length(truth>=1)) abline(v=truth[1],col=2)
    #}
  }
  #Evolution (move) (SIR doesn't include this step.)
  
  
  return(list(phistory = particles_history, whistory = weight_history)) ##output summaries (should include more info)
  #mean=M.st,var=V.st,l=logl,
}

#######TIME###############
tstep = .1
tmax =200
ts = seq(0, tmax, tstep)

n_trials = 1#don't really need
tr = 1.5
acq_ts = seq(1, tmax, by = tr)#acquisition times
n_acqs = length(acq_ts)
trs = which(is.element(ts, acq_ts))#acquisition times in units of ts

n_pres = 25#number of presentations
#pres_ts = seq(0, tmax, by = n_pres)#presentation times
pres_ts=seq(0,tmax,tmax/n_pres)[-(n_pres+1)] # evenly spaced presentation

########parss############
#from lindquist's paper on digamma hrfs
#parss = NULL
#parss$a1 = 6
#parss$a2 = 16
#parss$b1 = 1
#parss$b2 = 1
#parss$c = 1/6
#parss$A = 1


parss=NULL
parss$max.bold=6
parss$a1=2
parss$a2=2
parss$b1=1
parss$b2=2
parss$g=.3

parss$sigma=0.01

parss$sigma_sys = .05
parss$sigma_meas =.005

prior = NULL
prior[1] =5.85#max.bold mean
prior[2] = .5#max.bold sd
prior[3] = 2.25#a1 mean...
prior[4] = .5
prior[5] = 2.15#a2
prior[6] = .5
prior[7] = 1.15 #b1
prior[8] = .3
prior[9] = 2.15 #b2
prior[10] = .5
prior[11] = .325 #g
prior[12] = .1
prior[13] = 0 #min
prior[14] = 4  #max

sys_noise = c(.01, .01, .01, .01, .01, .01)

#prior$A_sd = .05 #expected sd in in parameter A

#########GENERATE TRUTH##########
truth = hrf_sim_truth(truth, tmax, tr, ts, tstep, pres_ts, acq_ts, n_acqs, parss)
y = hrf_sim_measurements(truth, tmax, tr, ts, tstep, pres_ts, acq_ts, n_acqs, parss)

#########PLOT SIMULATION#####################
par(mfrow=c(1,1))
plot(ts,truth,type="l",lwd=4,xlab="Time",ylab="BOLD Response", xlim = c(0, 100))#black = truth
abline(h=0,lty=3)
abline(v=acq_ts,lty=3,col="red")#red lines = Acquisition (measurement) times
points(acq_ts,y,col="red",pch=16)#red pts = Measurements
abline(v=pres_ts,lwd=3,col="green")#green = presentation times

#########MAIN###########

N = 300 #number of particles

ptm = proc.time()
out = hrf_sir(N, y, truth, parss, prior, ts, acq_ts, pres_ts, sys_noise,PLOT = F)

proc.time() - ptm
#####RECONSTRUCT HRF FROM PARTICLES#######
bs = out$mean


########PLOT PARTICLE FILTER ESIMATES###########
parvals = as.numeric(parss)
parnames = c('max.bold', 'a1', 'a2', 'b1', 'b2', 'g')
par(mfcol = c(3,3), oma = c(0,0,2,0))
tt = c(1,2,5,length(trs)/6, 2*length(trs)/6, 3*length(trs)/6,4*length(trs)/6,5*length(trs)/6,6*length(trs)/6)
#trs
for(k in 1:6){
  for(t in tt){
    hist(out$phistory[,k,t],100,col = 1, main = paste('t=',t), xlim=c(min(out$phistory[,k,]), max(out$phistory[,k,])))
    abline(v = parvals[k], col = 2)
  }
  temp_string = paste('Particle histogram', parnames[k])
  title(main = temp_string, outer = TRUE)
}

###############Plot weights vs particles####################################

par(mfcol = c(3,3), oma = c(0,0,2,0))

for(k in 1:6){
  par(mfcol = c(3,3), oma = c(0,0,2,0))
  for(t in tt){
    plot(out$phistory[,k,t], out$whistory[,t],main = paste('t=',t), xlim=c(min(out$phistory[,k,]), max(out$phistory[,k,])) )
    abline(v = parvals[k], col = 2)
    }

  temp_string = paste('weights vs particles', parnames[k])
  title(main = temp_string, outer = TRUE)
}
###########Plot means##########


par(mfcol = c(2,3), oma = c(0,0,2,0))
ms = matrix(NA, nrow = 6, ncol = length(trs))
ql = matrix(NA, nrow = 6, ncol = length(trs))
qu = matrix(NA, nrow = 6, ncol = length(trs))

xs = 1:length(trs)
for(i in 1:length(trs)){
  for(k in 1:6){
    ms[k,i] = mean(out$phistory[,k,i])
    ql[k,i] = quantile(out$phistory[,k,i], .25)
    qu[k,i] = quantile(out$phistory[,k,i], .75)
  }
}

for(k in 1:6){
  temp_string = paste('particle means', parnames[k])

  plot(xs, ms[k,], ylim=c(min(out$phistory[,k,]), max(out$phistory[,k,])),  ylab = 'particle estimates', main = temp_string )#means
  par(new = T)
  plot(xs, ql[k,], ylim=c(min(out$phistory[,k,]), max(out$phistory[,k,])), pch = '-', ylab = 'particle estimates')
  par(new = T)
  plot(xs, qu[k,], ylim=c(min(out$phistory[,k,]), max(out$phistory[,k,])), pch = '-', ylab = 'particle estimates')
  #temp_string = paste('particle means', parnames[k])
  #title(main = temp_string, outer = TRUE)
  abline(parvals[k],0, col = 2)
  title(main = 'Means, 1st and 3rd quartiles', outer = TRUE)
}


######reconstruct hrf with new parameters
newpar = parss
endt = length(out$phistory[1,1,])
newpar$max.bold = mean(out$phistory[,1, endt])
newpar$a1 = mean(out$phistory[,2,endt])
#newpar$a1 = 3
newpar$a2 = mean(out$phistory[,3,endt])
newpar$b1 = mean(out$phistory[,4,endt])
newpar$b2 = mean(out$phistory[,5,endt])
newpar$g = mean(out$phistory[,6,endt])

hrf_est = hrf_sim_truth(truth, tmax, tr, ts, tstep, pres_ts, acq_ts, n_acqs, newpar)

par(mfrow=c(1,1))
plot(ts,truth,type="l",lwd=1,xlab="Time",ylab="BOLD Response", xlim = c(100, 200))#black = truth
abline(h=0,lty=3)
abline(v=acq_ts,lty=5,col="red")#red lines = Acquisition (measurement) times
points(acq_ts,y,col="red",pch=17)#red pts = Measurements
abline(v=pres_ts,lwd=3,col="green")#green = presentation times
points(ts, hrf_est, col = "blue", pch = 4)#blue = estimation
legend(180,.135,c('truth', 'acquisitions' ,'measurements', 'presentations', 'estimation'),lty=c(1,5,1,1,1),col=c('black', 'red', 'red', 'green', 'blue') )


