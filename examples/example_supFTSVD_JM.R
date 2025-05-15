# Number of subjects
n<-500
# feature dimension
pdim<-c(200)
# No of rank-1 components
r<-3
model_rank<-r

# Weight for components
lmd_val<-c(5.20,4.80,3.35)

# noise variance in tensor
Tau2<-c(0.1)

# Parameters for supervised component
Eta2<-c(1,1,1)

#diag(Eta2)<-c(1.85,2.50,1.3,1.60,1.25,3.5)
gam_par<-matrix(c(2.2,1.6,2.5,3.6,-1.50,2.6),ncol=3)
gam_par<-apply(gam_par,2,function(u){u/sqrt(sum(u^2))})

# Grid of Time points
nres<-101
Time<-seq(0,1,length.out=nres)

# Singular Function
PhiFunc<-list(function(x){(8+(6*x^8)-(3*x^2)-(4*x^3))/sqrt(45.61)},
              function(x){(10*x^2/exp(x^5))/(sqrt(10-10*exp(-2)))},
              function(x){sqrt(2)*sin(2.5*pi*x)})

PhiF<-sapply(1:r,function(k){PhiFunc[[k]](Time)})

# Feature loading
set.seed(pdim)
Bval<<-sapply(1:r, function(b){runif(pdim)})
bval<<-Bval*outer(rep(1,pdim),1/apply(Bval,2,norm,type="2"))

# controlling the signal ot noise ratio
cmp_var<-((lmd_val^2)*Eta2)
CmpV<-Reduce(`+`,lapply(1:r,function(k){cmp_var[k]*(outer(bval[,k],PhiF[,k]))^2}))


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
surv_time<-(-log(runif(n,0,1)))/(lmd*exp(as.numeric(cbind(Vmat,ZetaSL)%*%(alp_par))))
summary(surv_time)
censor_time<-runif(n,0,5)
survT<-apply(cbind(surv_time,censor_time),1,min)
#summary(survT)
cenI<-apply(cbind(surv_time,censor_time),1,which.min)-1
cen_indx<-which(survT>1)
survT[cen_indx]<-1
cenI[cen_indx]<-1
cenP<-mean(cenI)
cenP
sum(as.numeric(survT==1))


# training m-omics data
m_i<-sample(5:10,n,replace = TRUE)
bl_time<-sapply(survT, function(u){runif(1,0,u)})
#tr_obsTIME<-lapply(1:n, function(i){c(bl_time[i],sort(runif(m_i[i]-2,bl_time[i],survT[i])),survT[i])})
tr_obsTIME<-lapply(1:n, function(i){c(0,sort(runif(m_i[i]-2,0,survT[i])),survT[i])})
gen_dataCFS<-omics_data_gen_surv(m_i = m_i,Zeta= ZetaSL,obsTIME = tr_obsTIME,Xi = bval,
                                 PsiF = PhiFunc,sing_val = lmd_val,Data_Var = Tau2,surv_time=survT)


fit_model<-supFTSVD_JM(datlist = gen_dataCFS$data,
                       response=Vmat, interval = c(0,1), r = model_rank,resolution=50, CVPhi=TRUE, K=5, cvT=5, smooth=round(exp(seq(-7,5,length.out=20)),3),
                       surv_time=survT,
                       censor_status=cenI,
                       maxiter=100, epsilon=1e-5,KInd=NULL,rsvd_seed=100,conv_criteria = "cond_lik",
                       survX = Vmat,scale = TRUE,constant_hazard = TRUE)

