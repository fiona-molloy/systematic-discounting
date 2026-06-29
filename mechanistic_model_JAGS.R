# Load in required packages 
library("rjags")
library("runjags")
# #######################
# # Objects in the order they appear ##
# subjects is the number of subjects
# trials is a vector of length subjects containing the number of trials each subject completed
# Pd is a matrix of (maximum) trials by subjects containing the observed choice
# muPd is a matrix of (maximum) trials by subjects containing the probability of choosing a larger-later option
# m is a vector of length subjects containing each subject's m estimate
# Vd is a matrix of (maximum) trials by subjects containing the subjective (discounted) value of the larger-later option
# V1 is a matrix of (maximum) trials by subjects containing the subjective (discounted) value of the smaller-sooner option, if delay = 0, this is equivalent to the objective value
# r is a matrix of (maximum) trials by subjects containing the objective value of larger-later option
# k is a vector of length subjects containing the subject's estimated discounting rate (k)
# time is a matrix of (maximum) trials by subjects containing the delay in days of the larger-later option
# mean_pd is a vector of length subjects calculating the overall probability of choosing larger-later
# k_hyper and m_hyper are the estimated hyperparameters governing the means of k and m respectively
# k_sigma and m_sigma are the estimated hyperparameters governing the precisions of k and m respectively (Sm is freely estimated and used to compute sigma_m)
# ####################

# JAGS Model 
model.hyperbolic = "
model{

# Likelihood
for (s in 1:subjects){
for(t in 1:trials[s]){
Pd[t,s] ~ dbern(muPd[t,s]*.99999999+0.00000001)
muPd[t,s] = 1/(1.00000001+exp(-m[s]*(Vd[t,s]-V1[t,s]))) 
}
}

# Hyperbolic discounting function
for(s in 1:subjects){
for(t in 1:trials[s]){
Vd[t,s] = r[t,s]/(1+k[s]*time[t,s])
}
# Overall probability of choosing larger-later 
mean_pd[s]=sum(muPd[1:trials[s],s])/trials[s]

}


# Priors
for(s in 1:subjects){
k[s] <- exp(logK[s])
logK[s] ~ dnorm(k_hyper,sigma_k) 
m[s] ~ dnorm(m_hyper,sigma_m) T(0,)
}

k_hyper ~ dnorm(0,.01)
m_hyper ~ dunif(0,10)

sigma_k ~ dgamma(.01,.01)
sigma_m <- 1/pow(Sm,2)
Sm ~ dunif(0,20)

}
"

# Save observed data as a list
dat.hyperbolic = list(Pd=Pd,time=time,V1=V1,r=r,trials=trials,subjects=subjects)

# fit model (summarize = F selected for speed, can calculate later)
results <- run.jags(model=model.hyperbolic,data=dat.hyperbolic,n.chains=3,monitor= c("k_hyper","m_hyper","logK","m"),summarise=FALSE)



