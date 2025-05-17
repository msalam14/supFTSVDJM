# Number of subjects
n<-500
# Parameters for supervised component
Eta2<-c(1,1,1)

#diag(Eta2)<-c(1.85,2.50,1.3,1.60,1.25,3.5)
gam_par<-matrix(c(2.2,1.6,2.5,3.6,-1.50,2.6),ncol=3)
gam_par<-apply(gam_par,2,function(u){u/sqrt(sum(u^2))})

# Survival parameters
alp_par<-matrix(c(3.50,2.60,-0.15,-0.20,0.15),ncol=1)
lmd<-0.10
## Subject-specific covariates # ADAS13
set.seed(n)
Vmat<<-cbind(round(runif(n),2),round(rbeta(n,2.5,1.5),2))
EAval<-Vmat%*%gam_par
miv_subE<-sapply(Eta2,function(u){rnorm(n,mean=0,sd=sqrt(u))})
ZetaSL<-EAval+miv_subE

# Training data set
# Survival data generation
seed_n<-4
set.seed(seed_n*24) # for both survival and longitudinal data
surv_time<-(-log(runif(n,0,1)))/(lmd*exp(as.numeric(cbind(Vmat,
                                                          ZetaSL)%*%(alp_par))))
summary(surv_time)
censor_time<-runif(n,0,5)
survT<-apply(cbind(surv_time,censor_time),1,min)
#summary(survT)
cenI<-apply(cbind(surv_time,censor_time),1,which.min)-1
cen_indx<-which(survT>1)
survT[cen_indx]<-1
cenI[cen_indx]<-1

sl_val<-sl_lik(par=alp_par,Xmat = cbind(Vmat,ZetaSL),surv_time = surv_time,
       censor_status = cenI,bz_hat = rep(lmd,n),cum_bzhat = lmd*surv_time)
