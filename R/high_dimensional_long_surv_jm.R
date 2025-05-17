#' Supervised functional singular value decomposition
#'
#' This function performs low-rank decomposition of high-dimensional
#' multivariate function following the methodology proposed in Alam (2024+)
#'
#' @param datlist a list with n elements, each is a matrix of dimension p+1
#' times m_i. The first row represents time points where measurements were
#' obtained. Next p rows represent the function values observed at the time
#' points
#' @param response a matrix of dimension n by q contains the design matrix for
#' supervision of the subject loadings
#' @param interval represents the domain of the function. The default is NULL.
#' @param r number of low-rank components involved in the assumed model
#' @param resolution grid size on which the singular functions will be estimated
#' @param CVPhi a logical scalar representing whether cross-validation (CV) will be
#' performed for determining the smoothness parameters involved with singular
#' functions
#' @param K number of folds for the CV step
#' @param cvT number of initial iterations at which CV will be performed
#' @param smooth Smoothing parameter for RKHS norm when CVPhi is FALSE. With
#' CVPhi=TRUE, a vector of numeric values for smoothness parameters for grid
#' search
#' @param surv_time a vector of N elements  represening observed survival times
#' @param censor_status a vector of censoring indicators (1 indicates censoring)
#' @param survX time-independent covariates specific to the survival sub-model
#' @param mc_ss Monte Carlo sample size for the E-step
#' @param maxiter Maximum number of iteration. Default: 20.
#' @param epsilon Convergence threshold for objective function value
#' @param KInd is a vector indices to construct the folds for CV. The default is
#' NULL that indicates that folds will be constructed randomly inside
#' @param rsvd_seed seed for random SVD (default is 100)
#' @param conv_criteria on which quantity the algorithm would be terminated (the
#' default is "par", indicating maximum difference in parameters will be used)
#' @param scale a logical indicating whether unit norm and scale for subject
#' loading components will be used (default is TRUE)
#' @param constant_hazard a logical indicting whether a constant hazard will be
#' used (default is TRUE). If FALSE, a piecewise constant hazard will be fitted
#' @param change_points the change_points to be used when a piecewise constant
#' hazard function is fitted. The default is NULL, indicating deciles from the
#' observed survival times will be used as change points.
#' @return a list with following objects:
#' \itemize{
#'  \item A.hat : a matrix of dimension N by K containing the values of
#'  conditional means for subjects
#'  \item B.hat : a matrix of dimension P by K representing the subject loadings
#'  \item Phi.hat : a matrix of dimension resolution by r representing the
#'  estimated singular functions
#'  \item time : a vector of time points where Phi.hat values are estimated
#'  \item ConU : a matrix of N by K, containing conditional mean obtained at the
#'  final iteration
#'  \item Lambda : a vector of r lambda values for weights (CP decomposition) of
#'   low-rank components
#'  \item Gamma : a matrix of q by r representing the estimated beta parameters
#'  \item Xb=Xb : deterministic parts of the conditional means
#'  \item r2_stat : goodness of fit statistics, R^2
#'  \item Sigma2R : r estimated variance parameters involving with subject
#'   loading models
#'  \item Sigma2 : estimated tensor noise variance
#'  \item cnvr_thresh : thresh of the convergence criterion from the fintal step
#'  \item CoxPH_obj : estimated baselize hazard, cumulative hazard, estimated
#'  parameters, and other quantities related to the survival sub-model
#'  \item cnvrgnt : change of objective function value at different iterations
#'  \item Scale : value of the argument Scale (necessary for other functions)
#'  \item Constant_Hazard : value fo the Constant_Hazard argument (necessary fot
#'  prediction)
#' }
#' @example examples/example_supFTSVD_JM.R
#' @export
supFTSVD_JM <- function(datlist, response=NULL, interval = NULL, r = 3,
                        resolution=100, CVPhi=FALSE, K=5, cvT=5, smooth=1e-8,
                        surv_time,censor_status,survX=NULL,mc_ss=1000,maxiter=20,
                        epsilon=1e-4,KInd=NULL,rsvd_seed=100,conv_criteria="par",scale=TRUE,
                        constant_hazard=TRUE,change_points=NULL){
  n = length(datlist)
  p = nrow(datlist[[1]])-1
  event_status<-1-censor_status

  # response in matrix form
  if(is.null(response)){
    response<-matrix(rep(1,n),ncol=1)
  } else if(!is.matrix(response)){
    response<-matrix(response,ncol=1)
  }

  # storing original response
  ires<-response


  # Calculate range.
  timestamps.all = do.call(c,lapply(datlist, FUN=function(u){u[1,]}))

  timestamps.all = sort(unique(timestamps.all))
  if (is.null(interval)){
    interval = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  }

  # rescale the time to 0-1.
  input.time.range = c(timestamps.all[1], timestamps.all[length(timestamps.all)])
  for (i in 1:n){
    datlist[[i]][1,] = (datlist[[i]][1,] - input.time.range[1]) / (input.time.range[2] - input.time.range[1])
  }
  interval = (interval - input.time.range[1]) / (input.time.range[2] - input.time.range[1])

  ti = lapply(1:n,FUN=function(u){
    temp = 1 + round((resolution-1) * (datlist[[u]][1,] - interval[1]) / (interval[2] - interval[1]))
    temp[which(temp<=0 | temp>resolution)] = 0
    temp
  })

  y0=do.call(c,lapply(1:n,FUN=function(u){
    as.vector(t(datlist[[u]][-1,which(ti[[u]]>0)]))
  }))

  y<-y0

  tipos<-lapply(ti,function(u){u>0})

  Lt<-lapply(datlist,FUN = function(u){u[1,]})
  ind_vec<-rep(1:n,sapply(Lt,length))

  tm <- unlist(Lt)
  Kmat <- bernoulli_kernel(tm, tm)
  Kmat_output <- bernoulli_kernel(seq(interval[1],interval[2],length.out = resolution), tm)

  # Initial values
  ## initialization of b
  data.unfold = do.call(cbind,lapply(datlist,FUN=function(u){u[-1,]}))

  set.seed(rsvd_seed)
  b.initials <- rsvd(data.unfold, k=r)$u
  b.hat = as.matrix(b.initials[,1:r])

  ## initialization of a
  if(r==1){
    a.hat<-as.matrix(sapply(1:n, function(i){
      colMeans(t(datlist[[i]][-1,])%*%b.hat)
    }))
  } else{
    a.hat<-t(sapply(1:n, function(i){
      colMeans(t(datlist[[i]][-1,])%*%b.hat)
    }))
  }


  ## initialize U
  Ahat<-a.hat

  ## Initialize for gamma
  Resp<-response #cbind(1,response)

  if(scale){
    gammaP<-apply(sapply(1:r,function(k){
      as.numeric((solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(Ahat[,k],ncol = 1)))
    }),2,function(u){u/sqrt(sum(u^2))})
  } else{
    gammaP<-sapply(1:r,function(k){
      as.numeric((solve(t(Resp)%*%Resp))%*%(t(Resp)%*%matrix(Ahat[,k],ncol = 1)))
    })
  }


  # for all subject, rows correspond to subjects
  Xb<-Resp%*%gammaP

  # Prediction of conditional mean using initial values
  if(scale){
    LAhat<-sqrt(apply(Ahat-Xb,2,var))
    conU<-sapply(1:r,function(k){(Ahat[,k]-Xb[,k])/LAhat[k]})
  } else{
    LAhat<-rep(1,r)
    conU<-(Ahat-Xb)
  }

  Ahat<-Xb+conU


  ## initialize for sigma_r
  if(scale){
    sigR<-rep(1,r)
    sigma_R<-diag(r)
    diag(sigma_R)<-1/sigR
  } else{
    sigR<-colMeans((Ahat-Xb)^2)
    sigma_R<-diag(r)
    diag(sigma_R)<-1/sigR
  }


  # Initialize xi function
  # compress the data and apply function PCA.
  phi.hat<-sapply(1:r, function(k){
    Ly<-lapply(1:n,FUN = function(i){
      (LAhat[k]*Ahat[i,k])*as.numeric(b.hat[,k]%*%datlist[[i]][2:(p+1),])
    })

    if(CVPhi){
      phiH = cv_freg_rkhs(Ly, Ahat[,k], ind_vec, Kmat, Kmat_output, smooth=smooth, kfold = K,KInd=KInd)[[2]]
    } else{
      phiH = freg_rkhs(Ly, Ahat[,k], ind_vec, Kmat, Kmat_output, smooth=smooth)
    }
    #phiH = freg_rkhs(Ly, a.hat, ind_vec, Kmat, Kmat_output, smooth=smooth)
    (phiH / sqrt(sum(phiH^2)))*sqrt(resolution)
  })

  t_arg<-seq(interval[1],interval[2],length.out = resolution)

  allT<-do.call(c,Lt)
  pXI<-apply(phi.hat,2,function(u){
    approx(t_arg,u,xout = allT,rule=2)$y
  })


  # initial values for error variance
  mi<-as.numeric(sapply(Lt,length))
  Mi<-cumsum(mi)

  prdXI<-lapply(1:n,function(w){
    if(w ==1){
      as.matrix(pXI[1:Mi[1],])
    } else{
      as.matrix(pXI[(Mi[w-1]+1):(Mi[w]),])
    }
  })

  # Initialization of Sigma2
  cwY2<-sapply(1:r, function(k){
    as.numeric(do.call(cbind,lapply(1:n, function(i){
      outer(b.hat[,k],prdXI[[i]][,k])*(Ahat[i,k]*LAhat[k])
    })))
  })

  #LAhat<-matrix(as.numeric(coef(lm(yr~.-1,data=data.frame("yr"=as.numeric(data.unfold),cwY2)))))

  #prdY2<-rowSums(hadamard.prod(cwY2,outer(rep(1,nrow(cwY2)),as.numeric(LAhat))))
  prdY2<-rowSums(cwY2)

  sig<-mean((as.numeric(data.unfold)-prdY2)^2)


  # Necessary matrices
  Hmat<-lapply(1:n, function(i){
    sapply(1:r, function(k){
      LAhat[k]*as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))
    })
  })

  # for all subject, rows correspond to subjects
  HiB<-Resp%*%gammaP


  ## Subject wise super X matrix
  Hbeta<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(HiB[i,],ncol=1))
  })

  # Subject wise vectorized data
  Yvec<-lapply(1:n, function(i){
    as.numeric(datlist[[i]][-1,])
  })

  # Subject-specific conditional covariances
  sw_vcovY<-lapply(Hmat,function(u){
    solve(crossprod(u)*(1/sig)+sigma_R)
  })

  if(r==1){
    varU<-as.matrix(sapply(1:n, function(i){
      diag(sw_vcovY[[i]])
    }))
  } else{
    varU<-t(sapply(1:n, function(i){
      diag(sw_vcovY[[i]])
    }))
  }

  # another necessary matrix
  HUmat<-lapply(1:n, function(i){
    as.numeric(Hmat[[i]]%*%matrix(conU[i,],ncol=1))
  })


  # Initialization of survival parameters
  dp_dim<-ncol(survX)
  if(!is.null(survX)){
    surv_xmat<-cbind(survX,Ahat)
  } else{
    surv_xmat<-Ahat
  }
  srv_dat<-data.frame(surv_time,1-censor_status,surv_xmat)
  if(!is.null(survX)){
    colnames(srv_dat)<-c("Stime","Sevent",paste('X',1:dp_dim,sep="_"),paste('Component', 1:r,sep = "_"))
  } else{
    colnames(srv_dat)<-c("Stime","Sevent",paste('Component', 1:r,sep = "_"))
  }
  fit_s<-coxph(Surv(Stime,Sevent)~.,data = srv_dat)
  # unique event time
  evn_time<-surv_time[event_status==1]
  #ev_time<-sort(unique(evn_time))
  #ev_n<-sapply(ev_time, function(u){length(which(evn_time==u))})
  if(constant_hazard){
    ev_time<-sort(unique(evn_time))
  } else{
    if(is.null(change_points)){
      ev_time<-unique(as.numeric(c(0,quantile(surv_time,probs = seq(0.1,0.9,0.1)),max(surv_time)))) #sort(unique(evn_time))
    } else{
      ev_time<-change_points
    }

    bhz_time<-t(sapply(surv_time, function(u){
      a1<-which.max(ev_time[ev_time<u])
      an<-length(ev_time)
      if(a1==1){
        c(u,rep(0,an-2))
      } else{
        c(ev_time[2:a1]-ev_time[1:(a1-1)],u-ev_time[a1],rep(0,an-a1-1))
      }
    }))
  }

  an<-length(ev_time)
  ev_n<-sapply(2:length(ev_time), function(u){length(which(evn_time>ev_time[u-1]&evn_time<=ev_time[u]))})
  rsk_ev<-lapply(ev_time[-an], function(u){which(surv_time>=u)})
  if(constant_hazard){
    bz_hat<-coxph.detail(fit_s)$hazard
    hzr<-as.numeric(coef(lm(cumsum(bz_hat)~coxph.detail(fit_s)$time-1)))
    fitbz<-data.frame("hazard"=hzr*ev_time,"time"=ev_time)#basehaz(fit_s)
    fitbz_dat<-data.frame("hazard"=hzr*surv_time,"time"=surv_time)
    bzH<-rep(hzr,n)
  } else{
    bz_hat<-sapply(2:an, function(u){
      ev_n[u-1]/sum(exp((as.matrix(srv_dat[rsk_ev[[u-1]],-c(1,2)])%*%matrix(as.numeric(coef(fit_s)))))*(pmin(surv_time[rsk_ev[[u-1]]],ev_time[u])-ev_time[u-1]))
    })
    fitbz<-data.frame("hazard"=cumsum(bz_hat*(ev_time[2:an]-ev_time[1:(an-1)])),"time"=ev_time[-1])#basehaz(fit_s)
    survPOS<-sapply(surv_time, function(u){
      #a<-which.min((u-fitbz$time)^2)
      which.max(ev_time[ev_time<u])
    })
    fitbz_dat<-data.frame("hazard"=bhz_time%*%bz_hat,"time"=surv_time)
    bzH<-bz_hat[survPOS]
  }
  sg_par<-as.numeric(fit_s$coefficients)
  if(!is.null(survX)){
    sg_parD<-sg_par[1:dp_dim]
    sg_parL<-sg_par[(dp_dim+1):(dp_dim+r)]
  }



  sl_lik_val<-sl_lik(par=sg_par,Xmat=surv_xmat,surv_time=surv_time,censor_status=censor_status,
                     bz_hat=bzH,cum_bzhat=fitbz_dat[,1])


  ## objective function value using initial values
  obj_val<-Reduce(`+`,lapply(1:n, function(i){
    t1<-(-((p*mi[i])/2)*log(sig))-sum((1/2)*log(sigR))
    t2<-sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2/(sig))
    t3<-sum(diag((crossprod(Hmat[[i]])*(1/sig))%*%sw_vcovY[[i]]))
    t4<-sum(diag(sigma_R%*%sw_vcovY[[i]]))
    t5<-sum((conU[i,]^2/sigR))
    t1-(t2+t3+t4+t5)/2
  }))+sl_lik_val

  # initial return object
  ret_obj<-list("A.hat" = Ahat,
                "B.hat" = b.hat,
                "Phi.hat" = phi.hat,
                "Gamma"=gammaP,
                "Xb"=Xb,
                "conU"= conU,
                "Lambda"=LAhat,
                "Sigma2R"=sigR,
                "Sigma2"=sig,
                "cnvrgnt"=NULL,
                "CoxPH_obj"=list("hazard"=bz_hat,"cum_hazard"=fitbz,"coefficients"=sg_par,"infor_mat"=NULL,"surv_design"=surv_xmat))


  # Starting of the EM based estimation
  t<-1
  iter_dif<-NULL
  while(t<=(maxiter+cvT)){
    dif<-NULL
    # E step for EM algorithm
    # conditional mean for all subjects
    if(r==1){
      tconU<-as.matrix(sapply(1:n,FUN = function(i){
        sw_vcovY[[i]]%*%(t(Hmat[[i]]*(1/sig))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))

      mc_obj<-lapply(1:n,function(i){
        if(!is.null(survX)){
          xb_i<-exp(sum(survX[i,]*sg_parD)+sum(Xb[i,]*sg_parL))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_parL)))
        } else{
          xb_i<-exp(sum(Xb[i,]*sg_par))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_par)))
        }
        srv_ll<-sapply(1:mc_ss,function(u){
          aa<-exp(-rw_sam[u]*fitbz_dat[i,1]*xb_i)
          bb<-ifelse(aa<1e-6,1e-6,aa)
          ((bzH[i]*xb_i*rw_sam[u])^event_status[i])*(bb)
        })
        wt<-srv_ll/mean(srv_ll)
        sL<-colMeans(hadamard.prod(sam,outer(wt,rep(1,r))))
        svcov<-Reduce(`+`,lapply(1:mc_ss, function(j){
          outer(sam[j,]-sL,sam[j,]-sL)*wt[j]
        }))/mc_ss
        exp_rw<-mean(rw_sam*wt)
        rw_exp_rw<-colMeans(hadamard.prod(sam,outer((wt*rw_sam),rep(1,r))))
        xx_exp_rw<-Reduce(`+`,lapply(1:mc_ss, function(j){
          outer(sam[j,],sam[j,])*fitbz_dat[i,1]*rw_sam[j]*wt[j]
        }))/mc_ss
        list(sL,exp_rw,rw_exp_rw,xx_exp_rw,svcov)
      })
      conU<-as.matrix(sapply(1:n,function(i){mc_obj[[i]][[1]]}))
      expRW<-sapply(1:n,function(i){mc_obj[[i]][[2]]})
      uexpRW<-as.matrix(sapply(1:n,function(i){mc_obj[[i]][[3]]}))
      ixE<-lapply(1:n,function(i){mc_obj[[i]][[4]]})
      sw_vcov<-lapply(1:n,function(i){mc_obj[[i]][[5]]})

    } else{
      tconU<-t(sapply(1:n,FUN = function(i){
        sw_vcovY[[i]]%*%(t(Hmat[[i]]*(1/sig))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))

      mc_obj<-lapply(1:n,function(i){
        if(!is.null(survX)){
          xb_i<-exp(sum(survX[i,]*sg_parD)+sum(Xb[i,]*sg_parL))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_parL)))
        } else{
          xb_i<-exp(sum(Xb[i,]*sg_par))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_par)))
        }
        srv_ll<-sapply(1:mc_ss,function(u){
          aa<-exp(-rw_sam[u]*fitbz_dat[i,1]*xb_i)
          bb<-ifelse(aa<1e-6,1e-6,aa)
          ((bzH[i]*xb_i*rw_sam[u])^event_status[i])*(bb)
        })
        wt<-srv_ll/mean(srv_ll)
        sL<-colMeans(hadamard.prod(sam,outer(wt,rep(1,r))))
        svcov<-Reduce(`+`,lapply(1:mc_ss, function(j){
          outer(sam[j,]-sL,sam[j,]-sL)*wt[j]
        }))/mc_ss
        exp_rw<-mean(rw_sam*wt)
        rw_exp_rw<-colMeans(hadamard.prod(sam,outer((wt*rw_sam),rep(1,r))))
        xx_exp_rw<-Reduce(`+`,lapply(1:mc_ss, function(j){
          outer(sam[j,],sam[j,])*fitbz_dat[i,1]*rw_sam[j]*wt[j]
        }))/mc_ss
        list(sL,exp_rw,rw_exp_rw,xx_exp_rw,svcov)
      })
      conU<-t(sapply(1:n,function(i){mc_obj[[i]][[1]]}))
      expRW<-sapply(1:n,function(i){mc_obj[[i]][[2]]})
      uexpRW<-t(sapply(1:n,function(i){mc_obj[[i]][[3]]}))
      ixE<-lapply(1:n,function(i){mc_obj[[i]][[4]]})
      sw_vcov<-lapply(1:n,function(i){mc_obj[[i]][[5]]})
    }

    if(r==1){
      varU<-as.matrix(sapply(1:n, function(i){
        diag(sw_vcov[[i]])
      }))
    } else{
      varU<-t(sapply(1:n, function(i){
        diag(sw_vcov[[i]])
      }))
    }

    # M step
    # update of beta_k for all k
    eta_val<-NULL
    Anew<-conU+Xb
    for(k in 1:r){
      if(r>1){
        indx<-c(1:r)[-k]
      }

      # residual data
      if(r==1){
        rY<-lapply(1:n, function(i){
          datlist[[i]][-1,]
        })
      } else{
        rY<-lapply(1:n, function(i){
          datlist[[i]][-1,]-Reduce(`+`,lapply(indx, function(l){
            (LAhat[l]*Anew[i,l])*outer(b.hat[,l],prdXI[[i]][,l])
          }))
        })
      }

      # Beta update
      ## joint model beta_k estimation
      bt_scr<-t(sapply(1:n, function(i){

        xi<-do.call(rbind,lapply((prdXI[[i]][,k]), function(j){
          t(outer(Resp[i,],b.hat[,k])*j*LAhat[k])
        })) #%*%gammaP[,k]

        xib<-xi%*%gammaP[,k]


        zi<-as.numeric(rY[[i]]-((LAhat[k]*conU[i,k])*outer(b.hat[,k],prdXI[[i]][,k])))

        aa<-Reduce('+',lapply(seq_len(length(zi)), function(b){
          xi[b,]*(zi[b]-xib[b])*(1/sig)
        }))

        if(!is.null(survX)){
          bb<-Resp[i,]*event_status[i]*sg_parL[k]
          dd<-as.numeric(expRW[i]*fitbz_dat[i,1]*sg_parL[k]*exp(sum(survX[i,]*sg_parD)+sum(Xb[i,]*sg_parL)))*Resp[i,]
        } else{
          bb<-Resp[i,]*event_status[i]*sg_par[k]
          dd<-as.numeric(expRW[i]*fitbz_dat[i,1]*sg_par[k]*exp(sum(Xb[i,]*sg_par)))*Resp[i,]
        }

        aa + bb-dd
      }))

      btSD<-Reduce('+',lapply(1:n, function(i){

        xi<-do.call(rbind,lapply((prdXI[[i]][,k]), function(j){
          t(outer(Resp[i,],b.hat[,k])*j*LAhat[k])
        })) #%*%gammaP[,k]

        aa<-Reduce('+',lapply(seq_len(nrow(xi)), function(b){
          outer(xi[b,],xi[b,])*(1/sig)
        }))

        if(!is.null(survX)){
          dd<-as.numeric(expRW[i]*fitbz_dat[i,1]*(sg_parL[k]^2)*exp(sum(survX[i,]*sg_parD)+sum(Xb[i,]*sg_parL)))*outer(Resp[i,],Resp[i,])
        } else{
          dd<-as.numeric(expRW[i]*fitbz_dat[i,1]*(sg_par[k]^2)*exp(sum(Xb[i,]*sg_par)))*outer(Resp[i,],Resp[i,])
        }

        -aa-dd
      }))


      ngamma<-(gammaP[,k]-solve(btSD,colSums(bt_scr)))
      if(scale){
        new_gamma<-ngamma/sqrt(sum(ngamma^2))
        dif<-c(dif,sum((new_gamma-gammaP[,k])^2)/sum(gammaP[,k]^2))
        gammaP[,k]<-new_gamma
      } else{
        dif<-c(dif,sum((ngamma-gammaP[,k])^2)/sum(gammaP[,k]^2))
        gammaP[,k]<-ngamma
      }

      # updated Xb
      Xb<-Resp%*%gammaP
      Anew<-conU+Xb


      # Scaling factor
      scD<-sqrt(((LAhat[k]*Anew[,k])^2)+sapply(sw_vcov,function(u){u[k,k]}))

      if(r==1){
        srY<-lapply(1:n,function(i){
          (matrixcalc::hadamard.prod(outer(rep(LAhat[k]*Anew[i,k],p),rep(1,nrow(prdXI[[i]]))),rY[[i]]))*(1/scD[i])
        })
      } else{
        srY<-lapply(1:n,function(i){
          (matrixcalc::hadamard.prod(outer(rep(LAhat[k]*Anew[i,k],p),rep(1,nrow(prdXI[[i]]))),rY[[i]])-Reduce(`+`,lapply(indx,function(l){
            outer(b.hat[,l],prdXI[[i]][,l])*(sw_vcov[[i]][k,l]*LAhat[l])
          })))*(1/scD[i])
        })
      }

      # update of xi_b k

      # scaled Residual data
      sY<-sapply(1:p,function(b){
        do.call(c,lapply(srY, function(u){u[b,]}))
      })


      nbhat<-NULL
      for(b in 1:p){
        bX<-do.call(c,lapply(1:n,function(i){scD[i]*prdXI[[i]][,k]}))

        nbhat<-c(nbhat,sum(sY[,b]*bX)/sum(bX^2))
      }
      n_bhat<-nbhat/sqrt(sum(nbhat^2))
      dif<-c(dif,sum((n_bhat-b.hat[,k])^2))
      b.hat[,k]<-n_bhat

      # Update of psi function
      Lty = lapply(1:n,function(i){
        list(scD[i]*as.numeric(b.hat[,k]%*%srY[[i]]))
        #(matrixcalc::hadamard.prod(scF[[i]],rY[[i]]))))
      })

      if(CVPhi & t<=cvT){
        cvfit<-cv_freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=smooth,kfold = K,KInd=KInd)
        Smv<-cvfit[[1]]
        phiH = cvfit[[2]]
        eta_val<-c(eta_val,Smv)
      } else{
        phiH = freg_rkhs(Lty,scD, ind_vec, Kmat, Kmat_output, smooth=ifelse(CVPhi,Smv,smooth))
        eta_val<-c(eta_val,ifelse(CVPhi,Smv,smooth))
      }
      phi.new = (phiH / sqrt(sum(phiH^2)))*sqrt(resolution)
      dif<-c(dif,sum((phi.hat[,k] - phi.new)^2))
      phi.hat[,k] <- phi.new

      kpXI<-approx(t_arg,phi.new,xout = allT,rule=2)$y

      for(i in 1:n){
        if(i ==1){
          prdXI[[i]][,k]<-kpXI[1:Mi[1]]
        } else{
          prdXI[[i]][,k]<-kpXI[(Mi[i-1]+1):(Mi[i])]
        }
      }
    }

    ## update of survival parameters
    if(!is.null(survX)){
      Xbs<-exp(as.numeric(survX%*%matrix(sg_parD))+as.numeric(Xb%*%matrix(sg_parL)))
    } else{
      Xbs<-exp(as.numeric(Xb%*%matrix(sg_par)))
    }
    nexpRW<-Xbs*expRW
    #nexpRW<-exp(Anew%*%matrix(sg_par))
    if(constant_hazard){
      hzr<-sum(ev_n)/sum(surv_time*nexpRW)
      nbz_hat<-rep(hzr,length(ev_time))
      bz_dif<-sqrt(sum((nbz_hat-bz_hat)^2))/sqrt(sum((bz_hat)^2))
      bz_hat<-nbz_hat
      bzH<-rep(hzr,n)
      fitbz<-data.frame("hazard"=hzr*ev_time,"time"=ev_time)#basehaz(fit_s)
      fitbz_dat<-data.frame("hazard"=hzr*surv_time,"time"=surv_time)
    } else{
      nbz_hat<-sapply(2:length(ev_time), function(u){
        ev_n[u-1]/sum(nexpRW[rsk_ev[[u-1]]]*(pmin(surv_time[rsk_ev[[u-1]]],ev_time[u])-ev_time[u-1]))
      })

      bz_dif<-sqrt(sum((nbz_hat-bz_hat)^2))/sqrt(sum((bz_hat)^2))
      bz_hat<-nbz_hat
      bzH<-bz_hat[survPOS]
      fitbz<-data.frame("hazard"=cumsum(bz_hat*(ev_time[2:an]-ev_time[1:(an-1)])),"time"=ev_time[-1])#basehaz(fit_s)
      fitbz_dat<-data.frame("hazard"=bhz_time%*%bz_hat,"time"=surv_time)
      # bzH<-bz_hat[survPOS]
      # fitbz<-data.frame("hazard"=cumsum(bz_hat),"time"=ev_time)#basehaz(fit_s)
      # fitbz_dat<-fitbz[survPOS,]
    }


    if(r==1){
      nuexpRW<-matrix(sapply(1:n, function(i){
        (Xb[i,]*Xbs[i]*expRW[i])+(Xbs[i]*uexpRW[i,])
      }))
    } else{
      nuexpRW<-t(sapply(1:n, function(i){
        (Xb[i,]*Xbs[i]*expRW[i])+(Xbs[i]*uexpRW[i,])
      }))
    }


    if(!is.null(survX)){
      surv_scr<-t(sapply(1:n, function(i){
        (c(survX[i,],Anew[i,])*event_status[i])-(fitbz_dat[i,1]*c(((Xbs[i]*expRW[i])*survX[i,]),nuexpRW[i,]))
      }))
    } else{
      if(r==1){
        surv_scr<-matrix(sapply(1:n, function(i){
          (Anew[i,]*event_status[i])-(fitbz_dat[i,1]*nuexpRW[i,])
        }))
      } else{
        surv_scr<-t(sapply(1:n, function(i){
          (Anew[i,]*event_status[i])-(fitbz_dat[i,1]*nuexpRW[i,])
        }))
      }
    }

    #hadamard.prod(Anew,outer(event_status,rep(-1,r)))-hadamard.prod(nuexpRW,outer(fitbz_dat[,1],rep(-1,r)))


    Ix<-Reduce(`+`,lapply(1:n,function(i){outer(surv_scr[i,],surv_scr[i,])}))

    Ixy<-outer(colSums(surv_scr),colSums(surv_scr))

    Imat<-Ix-Ixy


    nsg_par<-sg_par+solve(Ix,colSums(surv_scr))

    sg_dif<-sqrt(sum((nsg_par-sg_par)^2))/sqrt(sum(sg_par^2))

    sg_par<-nsg_par

    if(!is.null(survX)){
      sg_parD<-sg_par[1:dp_dim]
      sg_parL<-sg_par[(dp_dim+1):(dp_dim+r)]
    }

    if(!is.null(survX)){
      surv_xmat<-cbind(survX,Anew)
    } else{
      surv_xmat<-Anew
    }

    sl_lik_val<-sl_lik(par=sg_par,Xmat=surv_xmat,surv_time=surv_time,censor_status=censor_status,
                       bz_hat=bzH,cum_bzhat=fitbz_dat[,1])


    bz_ll<-sum(log(bzH)*event_status)


    # update of sigma2_k
    if(scale){
      sk2_new<-rep(1,r)
    } else{
      sk2_new<-colMeans(conU^2+varU)
    }


    # to avoid sigma2k=0
    if(any(sk2_new<1e-10)){
      sk2_new[sk2_new<1e-10]<-1e-10
    }

    if(!scale){
      dif<-c(dif,((sk2_new-sigR)^2)/sigR^2)
    }


    sigR<-sk2_new

    diag(sigma_R)<-1/sigR

    ## update of lambda
    cwY2<-sapply(1:r, function(k){
      as.numeric(do.call(cbind,lapply(1:n, function(i){
        outer(b.hat[,k],prdXI[[i]][,k])*(Anew[i,k])
      })))
    })

    # adjustment for trace involving lambda in the conditional likelihood
    HmatL<-lapply(1:n, function(i){
      sapply(1:r, function(k){
        as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))
      })
    })


    xx<-(t(cwY2)%*%cwY2)+Reduce(`+`,lapply(1:n,function(i){
      hadamard.prod((t(HmatL[[i]])%*%HmatL[[i]]),sw_vcov[[i]])
    }))
    xy<-t(cwY2)%*%matrix(as.numeric(data.unfold))

    LAhat<-as.numeric(solve(xx)%*%xy)

    #LAhat<-as.numeric(coef(lm(yr~.-1,data=data.frame("yr"=as.numeric(data.unfold),cwY2))))


    # update of sigma2
    # Necessary matrices (updated)
    Hmat<-lapply(1:n, function(i){
      sapply(1:r, function(k){
        as.numeric(outer(b.hat[,k],prdXI[[i]][,k]))*LAhat[k]
      })
    })

    # for all subject, rows correspond to subjects
    HiB<-Resp%*%gammaP


    ## Subject wise super X matrix
    Hbeta<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(HiB[i,],ncol=1))
    })

    HUmat<-lapply(1:n, function(i){
      as.numeric(Hmat[[i]]%*%matrix(conU[i,],ncol=1))
    })

    nsig2<-sum(sapply(1:n, function(i){
      sum(diag((t(Hmat[[i]])%*%Hmat[[i]])%*%sw_vcov[[i]]))+
        sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2)
    }))/sum(sapply(Yvec,length))

    dif<-c(dif,((nsig2-sig)^2)/sig^2)
    sig<-nsig2

    # objective function
    nobj_val<-Reduce(`+`,lapply(1:n, function(i){
      t1<-(-((p*mi[i])/2)*log(sig))-sum((1/2)*log(sigR))
      t2<-sum((Yvec[[i]]-Hbeta[[i]]-HUmat[[i]])^2/(sig))
      t3<-sum(diag((crossprod(Hmat[[i]])*(1/sig))%*%sw_vcov[[i]]))
      t4<-sum(diag(sigma_R%*%sw_vcov[[i]]))
      t5<-sum((conU[i,]^2/sigR))
      t1-(t2+t3+t4+t5)/2
    }))+sum(eta_val)+sl_lik_val

    # stopping criterion
    par_dif<-max(dif)

    if(CVPhi & t<=cvT+1){
      relC<-epsilon+0.05
    } else{
      relC<-(nobj_val-obj_val)/abs(obj_val)
    }

    if(conv_criteria=="par"){
      c_crit<-par_dif
    } else{
      c_crit<-relC
    }


    # update of all quantities by the results of current iteration
    dif<-c(dif,sg_dif,bz_dif,bz_ll,sl_lik_val,nobj_val)
    iter_dif<-rbind(iter_dif,dif)
    obj_val<-nobj_val
    ret_obj<-list("A.hat" = Anew,
                  "B.hat" = b.hat,
                  "Phi.hat" = phi.hat,
                  "Gamma"=gammaP,
                  "Xb"=Xb,
                  "conU"= conU,
                  "Lambda"=LAhat,
                  "Sigma2R"=sigR,
                  "Sigma2"=sig,
                  "cnvrgnt"=iter_dif,
                  "CoxPH_obj"=list("hazard"=bz_hat,"cum_hazard"=fitbz,"coefficients"=sg_par,"infor_mat"=Ix,"surv_design"=surv_xmat))


    # Testing the convergences and termination
    if(t>cvT & c_crit<=epsilon)
      break

    # Required quantities for next iteration if not terminated
    t <- t+1

    # Subject-specific conditional covariances
    sw_vcovY<-lapply(Hmat,function(u){
      solve(crossprod(u)*(1/sig)+sigma_R)
    })

  }

  # elements of returning object
  PCname <- paste('Component', 1:r)
  colnames(ret_obj$A.hat) = PCname
  colnames(ret_obj$Xb) = PCname
  colnames(ret_obj$B.hat) = PCname
  colnames(ret_obj$Phi.hat) = PCname
  rownames(ret_obj$A.hat) = names(datlist)
  rownames(ret_obj$Xb) = names(datlist)
  rownames(ret_obj$B.hat) = rownames(datlist[[1]])[-1]

  ## Time point where phi functions are estiamted

  time.return = seq(interval[1],interval[2],length.out = resolution)
  time.return = time.return * (input.time.range[2] - input.time.range[1]) + input.time.range[1]

  # Goodness of fit summary using the fitted object
  # vectorized observed data
  obs_y<-do.call(c,Yvec)
  tss<-sum((obs_y-mean(obs_y))^2)
  # observed value of estimated singular function
  subSF<-lapply(Lt, function(t){
    do.call(cbind,lapply(1:r, function(k){
      approx(x=time.return,y=ret_obj$Phi.hat[,k],xout = t,rule = 2)$y
    }))
  })

  # using Y_hat
  fit_y<-sapply(1:r, function(k){
    do.call(c,lapply(1:n,function(i){
      (ret_obj$Lambda[k]*ret_obj$A.hat[i,k])*as.numeric((outer(ret_obj$B.hat[,k],subSF[[i]][,k])))
    }))
  })

  err_y<-sapply(1:r,function(k){
    if(k==1){
      sum((obs_y-fit_y[,1])^2)
    } else{
      sum((obs_y-rowSums(fit_y[,1:k]))^2)
    }
  })

  r2_y<-1-(err_y/tss)

  crr_y<-sapply(1:r,function(k){
    if(k==1){
      (cor(obs_y,fit_y[,1]))^2
    } else{
      (cor(obs_y,rowSums(fit_y[,1:k])))^2
    }
  })


  r2lm_y<-sapply(1:r,function(k){
    if(k==1){
      summary(lm(obs_y~fit_y[,1]-1))$r.square
    } else{
      summary(lm(obs_y~fit_y[,1:k]-1))$r.square
    }
  })

  # using Xb
  fit_xb<-sapply(1:r, function(k){
    do.call(c,lapply(1:n,function(i){
      (ret_obj$Lambda[k]*ret_obj$Xb[i,k])*as.numeric((outer(ret_obj$B.hat[,k],subSF[[i]][,k])))
    }))
  })

  err_xb<-sapply(1:r,function(k){
    if(k==1){
      sum((obs_y-fit_xb[,1])^2)
    } else{
      sum((obs_y-rowSums(fit_xb[,1:k]))^2)
    }
  })

  r2_xb<-1-(err_xb/tss)

  crr_xb<-sapply(1:r,function(k){
    if(k==1){
      (cor(obs_y,fit_xb[,1]))^2
    } else{
      (cor(obs_y,rowSums(fit_xb[,1:k])))^2
    }
  })

  r2lm_xb<-sapply(1:r,function(k){
    if(k==1){
      summary(lm(obs_y~fit_xb[,1]-1))$r.square
    } else{
      summary(lm(obs_y~fit_xb[,1:k]-1))$r.square
    }
  })

  r2_stat<-data.frame("R2_Yhat"=r2_y,
                      "r_Yhat"=crr_y,
                      "R2lm_Yhat"=r2lm_y,
                      "R2_Xb"=r2_xb,
                      "r_Xb"=crr_xb,
                      "R2lm_Xb"=r2lm_xb)

  # return object
  results = list("A.hat" = ret_obj$A.hat, "B.hat" = ret_obj$B.hat,
                 "Phi.hat" = ret_obj$Phi.hat, "time" = time.return,
                 "ConU" = ret_obj$conU,
                 "Lambda" = ret_obj$Lambda,
                 "Gamma"=ret_obj$Gamma,
                 "Xb"=ret_obj$Xb,
                 "r2_stat" = r2_stat,
                 "Sigma2R"=ret_obj$Sigma2R,
                 "Sigma2"=ret_obj$Sigma2,
                 "cnvrgnt"=ret_obj$cnvrgnt,
                 "cnvr_thresh" = c_crit,
                 "CoxPH_obj"=ret_obj$CoxPH_obj,
                 "Scale"=scale,
                 "Constant_Hazard"=constant_hazard)

  return(results)
}



#' Full likelihood function induced by the cox-ph model
#'
#' This function calculates the full likelihood for observed survival data with
#' right censoring.
#'
#' @param par a vector of parameters
#' @param Xmat a matrix of covariates
#' @param surv_time observed survival time
#' @param censor_status a indicator vector (1/0 means censoring or not)
#' @param bz_hat a vector with baselize hazard values at the observed times
#' @param cum_bzhat a vector with cumulative baselize hazard values
#' at the observed times
#' @returns value of the log-likelihood for the observed data
#' @example examples/example_sl_lik.R
#' @export
sl_lik<-function(par,Xmat,surv_time,censor_status,bz_hat,cum_bzhat){
  event_status<-1-censor_status
  n<-nrow(Xmat)
  sum(sapply(1:n, function(i){
    rw<-sum(Xmat[i,]*par)
    if(event_status[i]==1){
      log(bz_hat[i])+rw-(exp(rw)*cum_bzhat[i])
    } else{
      -(exp(rw)*cum_bzhat[i])
    }
  }))
}


#' Prediction using the fitted objects from supFTSVD_JM
#'
#' @param obj fitted object by supFTSVD
#' @param designM matrix of covariates for all subjects (rows represent subject)
#' @param mc_ss Monte Carlo sample size for the E-step approximations of
#' conditional expectations
#' @param surv_time a vector of N elements  represening observed survival times
#' @param censor_status a vector of censoring indicators (1 indicates censoring)
#' @param survX time-independent covariates specific to the survival sub-model
#' @param new_dat is a list objects; each element is observed data for a subjects.
#' The first row contains the time-points where measurements were taken. Thus,
#' each element has (p+1) rows. The default is NULL; it is only required when
#' subject specific prediction is the interest.
#' @param TimeG Grid in the time domain where prediction will be made.
#' @return a list with following objects:
#' \itemize{
#'  \item meanTRJ : predicted mean trajectories (using Xb)
#'  \item subTRJ : predicted subject-feature-specific trajectories (using zeta)
#'  \item subEFF : Ahat
#'  \item Chazard : data.frame("hazard"=grid_bz,"time"=TimeG)
#'  \item Surv_Prob : pred_survP
#'  \item Surv_Prov_Var : var_sP
#' }
#' @example examples/example_predict_supFTSVD_JM.R
#' @export
predict.supFTSVD_JM<-function(obj,designM,mc_ss=1000,surv_time,
                              censor_status,survX=NULL,new_dat=NULL,TimeG=NULL){

  # number of components in the fitted model
  r<-ncol(obj$A.hat)
  # number of subjects for which prediction will be made
  ns<-nrow(designM)
  # Prediction of subject loading
  Xb<-designM%*%obj$Gamma
  # Indicator of event
  event_status<-1-censor_status


  # Observed time points for new subjects
  if(!is.null(new_dat)){
    newT<-lapply(1:ns, function(i){
      new_dat[[i]][1,]
    })

    newSF<-lapply(newT, function(nT){
      do.call(cbind,lapply(1:r, function(k){
        approx(x=obj$time,y=obj$Phi.hat[,k],xout = nT,rule = 2)$y
      }))
    })

    Hmat<-lapply(1:ns, function(i){
      as.matrix(sapply(1:r,function(k){
        obj$Lambda[k]*as.numeric(outer(obj$B.hat[,k],newSF[[i]][,k]))
      }))
    })

    Hbeta<-lapply(1:ns, function(i){
      as.numeric(Hmat[[i]]%*%matrix(Xb[i,],ncol=1))
    })

    sigma_R<-diag(r)
    diag(sigma_R)<-1/obj$Sigma2R

    sw_vcovY<-lapply(Hmat,function(u){
      solve(crossprod(u)*(1/obj$Sigma2)+sigma_R)
    })


    # Subject wise vectorized data
    Yvec<-lapply(1:ns, function(i){
      as.numeric(new_dat[[i]][-1,])
    })

    # survival parameters
    scale<-obj$Scale
    constant_hazard<-obj$Constant_Hazard
    cfv_mat<-solve(obj$CoxPH_obj$infor_mat) # estimated covariance matrix of survival model parameters
    sg_par<-as.numeric(obj$CoxPH_obj$coefficients)
    if(constant_hazard){
      fitbz<-obj$CoxPH_obj$cum_hazard
      hzr<-unique(obj$CoxPH_obj$hazard)
      fitbz_dat<-data.frame("hazard"=hzr*surv_time,"time"=surv_time)
      bzH<-rep(hzr,ns)
    } else{
      ev_time<-c(0,obj$CoxPH_obj$cum_hazard$time)
      an<-length(ev_time)
      survPOS<-sapply(surv_time, function(u){
        which.max(ev_time[ev_time<u])
      })
      bhz_time<-t(sapply(surv_time, function(u){
        a1<-which.max(ev_time[ev_time<u])
        if(a1==1){
          c(u,rep(0,an-2))
        } else{
          c(ev_time[2:a1]-ev_time[1:(a1-1)],u-ev_time[a1],rep(0,an-a1-1))
        }
      }))

      bz_hat<-obj$CoxPH_obj$hazard
      bzH<-bz_hat[survPOS]
      fitbz<-data.frame("hazard"=cumsum(bz_hat*(ev_time[2:an]-ev_time[1:(an-1)])),"time"=ev_time[-1])#basehaz(fit_s)
      fitbz_dat<-data.frame("hazard"=bhz_time%*%bz_hat,"time"=surv_time)
    }

    dp_dim<-ncol(survX)
    if(!is.null(survX)){
      sg_parD<-sg_par[1:dp_dim]
      sg_parL<-sg_par[(dp_dim+1):(dp_dim+r)]
    }


    # conditional mean for all subjects
    if(r==1){
      tconU<-as.matrix(sapply(1:ns,FUN = function(i){
        sw_vcovY[[i]]%*%(t(Hmat[[i]]*(1/obj$Sigma2))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))


      conU<-as.matrix(sapply(1:ns,function(i){
        if(!is.null(survX)){
          xb_i<-exp(sum(survX[i,]*sg_parD)+sum(Xb[i,]*sg_parL))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_parL)))
        } else{
          xb_i<-exp(sum(Xb[i,]*sg_par))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_par)))
        }
        srv_ll<-sapply(1:mc_ss,function(u){
          aa<-exp(-rw_sam[u]*fitbz_dat[i,1]*xb_i)
          bb<-ifelse(aa<1e-6,1e-6,aa)
          ((bzH[i]*xb_i*rw_sam[u])^event_status[i])*(bb)
        })

        wt<-srv_ll/mean(srv_ll)
        colMeans(hadamard.prod(sam,outer(wt,rep(1,r))))
      }))
    } else{
      tconU<-t(sapply(1:ns,FUN = function(i){
        sw_vcovY[[i]]%*%(t(Hmat[[i]]*(1/obj$Sigma2))%*%matrix(Yvec[[i]]-Hbeta[[i]],ncol=1))
      }))

      conU<-t(sapply(1:ns,function(i){
        if(!is.null(survX)){
          xb_i<-exp(sum(survX[i,]*sg_parD)+sum(Xb[i,]*sg_parL))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_parL)))
        } else{
          xb_i<-exp(sum(Xb[i,]*sg_par))
          sam<-mvrnorm(n=mc_ss,mu = tconU[i,],Sigma = sw_vcovY[[i]])
          rw_sam<- as.numeric(exp(sam%*%matrix(sg_par)))
        }
        srv_ll<-sapply(1:mc_ss,function(u){
          aa<-exp(-rw_sam[u]*fitbz_dat[i,1]*xb_i)
          bb<-ifelse(aa<1e-6,1e-6,aa)
          ((bzH[i]*xb_i*rw_sam[u])^event_status[i])*(bb)
        })
        wt<-srv_ll/mean(srv_ll)
        colMeans(hadamard.prod(sam,outer(wt,rep(1,r))))
      }))
    }
    Ahat<-Xb+conU
  } else{
    Ahat<-Xb
  }

  # Initialization of survival parameters
  if(!is.null(survX)){
    surv_xmat<-cbind(survX,Ahat)
  } else{
    surv_xmat<-Ahat
  }

  if(!is.null(surv_time)){
    varQ<-apply(surv_xmat,1,function(u){u%*%(cfv_mat%*%matrix(u))})
  } else{
    varQ<-apply(obj$CoxPH_obj$surv_design,1,function(u){u%*%(cfv_mat%*%matrix(u))})
  }


  if(!is.null(surv_time)&!is.null(TimeG)){
    if(constant_hazard){
      grid_bz<-hzr*TimeG
    } else{
      bhz_grid<-t(sapply(TimeG, function(u){
        a1<-which.max(ev_time[ev_time<u])
        if(a1==1){
          c(u,rep(0,an-2))
        } else {
          c(ev_time[2:a1]-ev_time[1:(a1-1)],u-ev_time[a1],rep(0,an-a1-1))
        }
      }))

      grid_bz<-bhz_grid%*%bz_hat
    }
    pred_survP<-sapply(exp(surv_xmat%*%matrix(sg_par)),function(u){exp(-u*grid_bz)})
    var_sP<-matrixcalc::hadamard.prod(sapply(exp(surv_xmat%*%matrix(sg_par)),function(u){((exp(-u*grid_bz)*(-u*grid_bz))^2)}),sapply(varQ,rep,times=length(grid_bz)))
  } else if(!is.null(surv_time)&is.null(TimeG)){
    pred_survP<-sapply(exp(surv_xmat%*%matrix(sg_par)),function(u){exp(-u*fitbz$hazard)}) #survival::survfit(obj$CoxPH_obj,newdata=srv_dat)$surv
    var_sP<-matrixcalc::hadamard.prod(sapply(exp(surv_xmat%*%matrix(sg_par)),function(u){((exp(-u*fitbz$hazard)*(-u*fitbz$hazard))^2)}),sapply(varQ,rep,times=length(fitbz$hazard)))
  } else{
    pred_survP<-sapply(exp(obj$CoxPH_obj$surv_design%*%matrix(sg_par)),function(u){exp(-u*fitbz$hazard)}) #survival::survfit(obj$CoxPH_obj,newdata = data.frame(obj$A.hat))$surv
    var_sP<-matrixcalc::hadamard.prod(sapply(exp(obj$CoxPH_obj$surv_design%*%matrix(sg_par)),function(u){((exp(-u*fitbz$hazard)*(-u*fitbz$hazard))^2)}),sapply(varQ,rep,times=length(fitbz$hazard)))
  }

  if(is.null(TimeG)){
    list("meanTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],obj$Phi.hat[,k])*(Xb[i,k]*obj$Lambda[k])
      }))
    }),
    "subTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],obj$Phi.hat[,k])*(Ahat[i,k]*obj$Lambda[k])
      }))
    }),
    "subEFF" = Ahat,
    "Chazard" = fitbz,
    "Surv_Prob" = pred_survP,
    "Surv_Prov_Var" = var_sP
    )
  } else{
    predSF<-sapply(1:r, function(k){
      approx(x=obj$time,y=obj$Phi.hat[,k],xout = TimeG,rule = 2)$y
    })
    list("meanTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],predSF[,k])*(Xb[i,k]*obj$Lambda[k])
      }))
    }),
    "subTRJ"=lapply(1:ns,function(i){
      Reduce(`+`,lapply(1:r, function(k){
        outer(obj$B.hat[,k],predSF[,k])*(Ahat[i,k]*obj$Lambda[k])
      }))
    }),
    "subEFF" = Ahat,
    "Chazard" = data.frame("hazard"=grid_bz,"time"=TimeG),
    "Surv_Prob" = pred_survP,
    "Surv_Prov_Var" = var_sP
    )
  }
}

#' Compute the conditional probability given that subjects survived upto time t
#'
#' @param Obj object returned by predict.supFTSVD_JM
#' @param lm_time landmark time point where the subjects were alive
#' @export
cond_surv_prob<-function(Obj,lm_time){
  nsF<-ncol(Obj$Surv_Prob)
  nsT<-nrow(Obj$Surv_Prob)
  scF<-sapply(1:nsF, function(u){(Obj$Surv_Prov_Var[,u]/(((log(Obj$Surv_Prob[,u]))*Obj$Surv_Prob[,u])^2))[1]})
  lapply(lm_time, function(u){
    lmP<-which.min((Obj$Chazard$time-u)^2)
    cond_survP<-sapply(1:nsF,function(v){Obj$Surv_Prob[(lmP):nsT,v]/Obj$Surv_Prob[lmP,v]})
    var_csP<-sapply(1:nsF, function(v){((cond_survP[,v]*log(cond_survP[,v]))^2)*scF[v]})
    list("TimeG"=Obj$Chazard$time[(lmP):nsT],
         "cond_survP"=cond_survP,
         "var_csF"=var_csP)
  })
}

#' Select a rank for supFTSVD by a modified version of cross-validation approach
#' in Lock and Li (2018)
#'
#' @param datlist a list with n elements, each is a matrix of dimension p+1
#' times m_i. The first row represents time points where measurements were
#' obtained. Next p rows represent the function values observed at the time
#' points
#' @param response a matrix of dimension n by q contains the design matrix for
#' supervision of the subject loadings
#' @param interval represents the domain of the function. The default is NULL.
#' @param ranks a vector of ranks on which grid search will be made
#' @param resolution grid size on which the singular functions will be estimated
#' @param CVPhi a logical scalar representing whether cross-validation (CV) will be
#' performed for determining the smoothness parameters involved with singular
#' functions
#' @param K number of folds for the CV step
#' @param cvT number of initial iterations at which CV will be performed
#' @param smooth Smoothing parameter for RKHS norm when CVPhi is FALSE. With
#' CVPhi=TRUE, a vector of numeric values for smoothness parameters for grid
#' search
#' @param maxiter Maximum number of iteration. Default: 20.
#' @param epsilon Convergence threshold for objective function value
#' @param KInd is a vector indices to construct the folds for CV. The default is
#' NULL that indicates that folds will be constructed randomly inside
#' @param rsvd_seed supFTSVD uses random svd for computational advantage while
#' setting inial values for feature loading vector. This seed number will be used
#' before setting the initial values. It confirms the reproducibility of your results.
#' @param conv_criteria the default is par which sets the threshold on the
#' difference between updated and current estimate at M step to stop the iteration.
#' It is computationally expensive; thus we recommend to use "cond_lik" instead of "par".
#' The result is not sensitive.
#' @param rank_fold number of folds to be used in the cross-validation
#' @param rfolds_seed seed number to be used in random splitting of the data
#' @param rc_thresh threshold on test-data likelihood for rank selection
#' @return a list with following objects:
#' \itemize{
#'  \item cv_res : detailed results for different folds
#'  \item opt_r : optimal rank
#' }
#' @example examples/example_cv_rank_supFTSVD.R
#' @export
cv_rank_supFTSVD_JM<-function(datlist, response, interval = NULL,
                           ranks = c(1:5),resolution=100, CVPhi=FALSE, K=5, cvT=5,
                           smooth=1e-8,maxiter=20, epsilon=1e-4,KInd=NULL,
                           surv_time,censor_status,survX=NULL,mc_ss=1000,
                           rsvd_seed=100,conv_criteria="par",
                           scale=TRUE,constant_hazard=TRUE,change_points=NULL,
                           rank_fold=5,rfolds_seed=50,rc_thresh=0.05,stratum=NULL,
                           par_computing=TRUE,n_core=5){
  n<-length(datlist)
  if(par_computing){
    library(doParallel)
    library(foreach)
    library(parallel)
    library(doSNOW)
    cls<-parallel::makeCluster(n_core)
    registerDoParallel(cls)
  }
  if(!is.null(stratum)){
    grp_frq<-as.numeric(table(stratum))
    lv_grp<-as.character(unique(stratum))
    set.seed(rfolds_seed)
    SJ_ind<-do.call(c,lapply(seq_len(length(lv_grp)),function(i){
      sample(1:grp_frq[i],grp_frq[i],replace=FALSE)
    }))
    rfolds_indx<-do.call(c,lapply(seq_len(length(lv_grp)),function(i){
      sj_ind<-sample(1:grp_frq[i],grp_frq[i],replace=FALSE)
      c(rep(1:rank_fold,each=grp_frq[i]%/%rank_fold),sample(1:rank_fold,grp_frq[i]%%rank_fold,replace = FALSE))
    }))
  } else{
    set.seed(rfolds_seed)
    SJ_ind<-sample(1:n,n,replace=FALSE)
    rfolds_indx<-c(rep(1:rank_fold,each=n%/%rank_fold),sample(1:rank_fold,n%%rank_fold,replace = FALSE))
  }
  rank_res<-NULL
  for(r in ranks){
    if(par_computing){
      crank<-foreach(g=1:rank_fold,
                     .combine = "c",
                     .packages = c("rsvd","matrixcalc","survival","MASS"),
                     .export = c("supFTSVD_JM",
                                 "predict.supFTSVD_JM",
                                 "norm_llik_supFTSVD",
                                 "sl_lik",
                                 "bernoulli_kernel",
                                 "cv_freg_rkhs",
                                 "freg_rkhs"),
                     .errorhandling = "remove") %dopar% {
                       train_indx<-sort(SJ_ind[which(rfolds_indx!=g)])
                       test_indx<-sort(SJ_ind[which(rfolds_indx==g)])
                       train_dat<-datlist[train_indx]
                       test_dat<-datlist[test_indx]
                       fit_model<-supFTSVD_JM(datlist=train_dat,
                                              response=response[train_indx,],
                                              interval = interval,
                                              r = r,
                                              resolution=resolution,
                                              CVPhi=CVPhi,
                                              K=K,
                                              cvT=cvT,
                                              smooth=smooth,
                                              maxiter=maxiter,
                                              epsilon=epsilon,
                                              KInd=KInd,
                                              surv_time=surv_time[train_indx],
                                              censor_status=censor_status[train_indx],
                                              survX=survX[train_indx,],
                                              mc_ss=mc_ss,
                                              rsvd_seed=rsvd_seed,
                                              conv_criteria=conv_criteria,
                                              scale=scale,
                                              constant_hazard=constant_hazard,
                                              change_points=change_points)
                       predFIT<-predict.supFTSVD_JM(obj = fit_model,
                                                    designM = response[test_indx,],
                                                    new_dat = test_dat,
                                                    mc_ss = mc_ss,
                                                    surv_time=surv_time[test_indx],
                                                    censor_status=censor_status[test_indx],
                                                    survX=survX[test_indx,],
                                                    TimeG=NULL)

                       test_Xb<-response[test_indx,]%*%fit_model$Gamma
                       test_conU<-predFIT$subEFF-test_Xb

                       lng_lik_val<-norm_llik_supFTSVD(datlist = test_dat,
                                                       x_matrix = response[test_indx,],
                                                       comp_beta = fit_model$Gamma,
                                                       subj_dev = test_conU,
                                                       feat_loading = fit_model$B.hat,
                                                       sing_func = fit_model$Phi.hat,
                                                       sing_func_arg = fit_model$time,
                                                       sing_val = fit_model$Lambda,
                                                       sigR = fit_model$Sigma2R,
                                                       sig = fit_model$Sigma2)

                       # Quautities associated with survival model
                       sg_par<-as.numeric(fit_model$CoxPH_obj$coefficients)
                       fitbz<-fit_model$CoxPH_obj$cum_hazard
                       appbz<-approx(x=fitbz$time,y=fitbz$hazard,xout=surv_time[test_indx],rule=2)
                       fitbz_dat<-data.frame("hazard"=appbz$y,"time"=appbz$x)
                       bzH<-approx(x=fit_model$CoxPH_obj$cum_hazard$time, y=fit_model$CoxPH_obj$hazard,xout=surv_time[test_indx],rule = 2)$y

                       if(!is.null(survX)){
                         test_survX<-cbind(survX[test_indx,],predFIT$subEFF)
                       } else{
                         test_survX<-predFIT$subEFF
                       }

                       sl_lik_val<-sl_lik(par=sg_par,
                                          Xmat=test_survX,
                                          surv_time=surv_time[test_indx],
                                          censor_status=censor_status[test_indx],
                                          bz_hat=bzH,
                                          cum_bzhat=fitbz_dat[,1])

                       lng_lik_val + sl_lik_val
                     }
    } else{
      crank<-sapply(1:rank_fold, function(g){
        train_indx<-sort(SJ_ind[which(rfolds_indx!=g)])
        test_indx<-sort(SJ_ind[which(rfolds_indx==g)])
        train_dat<-datlist[train_indx]
        test_dat<-datlist[test_indx]
        fit_model<-supFTSVD_JM(datlist=train_dat, response=response[train_indx,], interval = interval,
                               r = r,resolution=resolution, CVPhi=CVPhi, K=K, cvT=cvT,
                               smooth=smooth,maxiter=maxiter, epsilon=epsilon,KInd=KInd,
                               surv_time=surv_time[train_indx],
                               censor_status=censor_status[train_indx],
                               survX=survX[train_indx,],
                               mc_ss=mc_ss,
                               rsvd_seed=rsvd_seed,
                               conv_criteria=conv_criteria,scale=scale,
                               constant_hazard=constant_hazard,change_points=change_points)
        predFIT<-predict.supFTSVD_JM(obj = fit_model,
                                     designM = response[test_indx,],
                                     new_dat = test_dat,
                                     mc_ss = mc_ss,
                                     surv_time=surv_time[test_indx],
                                     censor_status=censor_status[test_indx],
                                     survX=survX[test_indx,],
                                     TimeG=NULL)

        test_Xb<-response[test_indx,]%*%fit_model$Gamma
        test_conU<-predFIT$subEFF-test_Xb

        lng_lik_val<-norm_llik_supFTSVD(datlist = test_dat,
                                        x_matrix = response[test_indx,],
                                        comp_beta = fit_model$Gamma,
                                        subj_dev = test_conU,
                                        feat_loading = fit_model$B.hat,
                                        sing_func = fit_model$Phi.hat,
                                        sing_func_arg = fit_model$time,
                                        sing_val = fit_model$Lambda,
                                        sigR = fit_model$Sigma2R,
                                        sig = fit_model$Sigma2)

        # Quautities associated with survival model
        sg_par<-as.numeric(fit_model$CoxPH_obj$coefficients)
        fitbz<-fit_model$CoxPH_obj$cum_hazard
        appbz<-approx(x=fitbz$time,y=fitbz$hazard,xout=surv_time[test_indx],rule=2)
        fitbz_dat<-data.frame("hazard"=appbz$y,"time"=appbz$x)
        bzH<-approx(x=fit_model$CoxPH_obj$cum_hazard$time, y=fit_model$CoxPH_obj$hazard,xout=surv_time[test_indx],rule = 2)$y

        if(!is.null(survX)){
          test_survX<-cbind(survX[test_indx,],predFIT$subEFF)
        } else{
          test_survX<-predFIT$subEFF
        }

        sl_lik_val<-sl_lik(par=sg_par,
                           Xmat=test_survX,
                           surv_time=surv_time[test_indx],
                           censor_status=censor_status[test_indx],
                           bz_hat=bzH,
                           cum_bzhat=fitbz_dat[,1])

        lng_lik_val + sl_lik_val
      })
    }
    if(r==1){
      rc<-1
    } else{
      rc<-(mean(crank)-mean(orank))/abs(mean(orank))
    }
    if(rc<rc_thresh)
      break
    orank<-crank
    rank_res<-rbind(rank_res,crank)
  }
  if(par_computing){
    rm(cls)
  }
  list("cv_res"=rank_res,"opt_r"=ranks[nrow(rank_res)])
}

#' Function for generating multi-omics (M-modal) data
#'
#' @param N number of subjects
#' @param n_i number of measurements (repeated data) from subject i
#' @param M number of modalities
#' @param XMatrix design matrix for subject specific mean
#' @param Beta regression parameters for subject specific mean
#' @param Xi feature singular vector for all modalities ( a list of M elements)
#' @param PsiF is a list of M elements; each element is also a list with r
#' elements which are functions
#' @param SubE_Var a list of M elements; each is a vector of r elements
#' represent the variance parameters for subject-loading model
#' @param Data_Var a vector of M elements represent measurement error variances
#' for different data modalities
#' @returns a list 3  elements
#' \itemize{
#'  \item data : a list of N elements, each is a matrix of P+1 by m_i. The first
#'  row is the observed time points.
#'  \item zeta : subject-specific loading used to data generation
#'  \item singF : singular functions
#' }
#' for subject-specific data
#' @export
omics_data_gen_surv<-function(m_i,Zeta,Xi,PsiF,sing_val,obsTIME=NULL,
                               Data_Var,Tgrid=seq(0,1,length.out=101),surv_time=NULL){

  n<-length(m_i)

  r<-ncol(Xi)

  if(is.null(obsTIME)){
    obs_time<-lapply(1:n,function(i){
      sort(runif(m_i[i],0,surv_time[i]))
    })
  } else{
    obs_time<-obsTIME
  }


  singF<-lapply(1:n, function(i){
    sapply(1:r,function(k){
      PsiF[[k]](obs_time[[i]]) #PsiF[[m]][[k]](obs_time[[i]])
    })
  })



  dat<-lapply(1:n,function(i){
    rbind(obs_time[[i]],
          Reduce(`+`,lapply(1:r,function(k){
            (sing_val[k]*Zeta[i,k])*outer(Xi[,k],singF[[i]][,k])
          })) + rnorm(nrow(Xi)*length(obs_time[[i]]),mean = 0,sd = sqrt(Data_Var)))
  })


  list("data"=dat,
       "zeta"=Zeta,
       "singF"=singF)
}


#' Full likelihood estimation of Cox-PH model
#'
#' This function maximizes full likelihood function induced by the Cox
#' Proportional hazard model. Essentially it takes the estimates from coxph of
#' survival package as initial estimates and obtain piecewise or constant hazard
#' estimates. It repeats this steps until the likelihood increases.
#' @param surv_time a vector of N elements  represening observed survival times
#' @param censor_status a vector of censoring indicators (1 indicates censoring)
#' @param survX time-independent covariates specific to the survival sub-model
#' @param maxiter Maximum number of iteration. Default: 25.
#' @param cnv_type on which quantity the algorithm would be terminated (the
#' default is "par", indicating maximum difference in parameters will be used)
#' @param cnv_crit threshold for iteration termination
#' @param constant_hazard a logical indicting whether a constant hazard will be
#' used (default is TRUE). If FALSE, a piecewise constant hazard will be fitted
#' @param change_points the change_points to be used when a piecewise constant
#' hazard function is fitted. The default is NULL, indicating deciles from the
#' observed survival times will be used as change points.
#' @returns a list with following items
#' \itemize{
#'  \item CoxPH_obj : estimated baselize hazard, cumulative hazard, estimated
#'  parameters, and other quantities related to the survival sub-model
#'  \item cnvrgn : change in parameter estimates and likelihood at different
#'  iterations
#'  \item rc_thresh : final threshold used to terminate the algorithm
#' }
#' @example examples/example_coxph_mle_fllik.R
#' @export
coxph_mle_fllik<-function(surv_time,censor_status,Xmat,max_iter=25,
                          cnv_type="Param",cnv_crit=1e-6,constant_hazard=TRUE,
                          change_points=NULL){
  n<-length(surv_time)
  p<-ncol(Xmat)
  # survival parameter initialization
  event_status<-1-censor_status
  evn_time<-surv_time[event_status==1]
  #ev_time<-sort(unique(evn_time))
  if(constant_hazard){
    ev_time<-sort(unique(evn_time))
  } else{
    if(is.null(change_points)){
      ev_time<-unique(as.numeric(c(0,quantile(surv_time,probs = seq(0.1,0.9,0.1)),max(surv_time)))) #sort(unique(evn_time))
    } else{
      ev_time<-change_points
    }
    bhz_time<-t(sapply(surv_time, function(u){
      a1<-which.max(ev_time[ev_time<u])
      an<-length(ev_time)
      if(a1==1){
        c(u,rep(0,an-2))
      } else{
        c(ev_time[2:a1]-ev_time[1:(a1-1)],u-ev_time[a1],rep(0,an-a1-1))
      }
    }))
  }

  an<-length(ev_time)
  ev_n<-sapply(2:length(ev_time), function(u){length(which(evn_time>ev_time[u-1]&evn_time<=ev_time[u]))})
  rsk_ev<-lapply(ev_time[-an], function(u){which(surv_time>=u)})
  surv_data<-data.frame("obsT"=surv_time,
                        "ev_stats"=1-censor_status,
                        Xmat)
  fit_s<-coxph(Surv(obsT,ev_stats)~.,data = surv_data)
  if(constant_hazard){
    bz_hat<-coxph.detail(fit_s)$hazard
    hzr<-as.numeric(coef(lm(cumsum(bz_hat)~coxph.detail(fit_s)$time-1)))
    fitbz<-data.frame("hazard"=hzr*ev_time,"time"=ev_time)#basehaz(fit_s)
    fitbz_dat<-data.frame("hazard"=hzr*surv_time,"time"=surv_time)
    bzH<-rep(hzr,n)
  } else{
    bz_hat<-sapply(2:an, function(u){
      ev_n[u-1]/sum(exp((as.matrix(surv_data[rsk_ev[[u-1]],-c(1,2)])%*%matrix(as.numeric(coef(fit_s)))))*(pmin(surv_time[rsk_ev[[u-1]]],ev_time[u])-ev_time[u-1]))
    })
    fitbz<-data.frame("hazard"=cumsum(bz_hat*(ev_time[2:an]-ev_time[1:(an-1)])),"time"=ev_time[-1])#basehaz(fit_s)
    survPOS<-sapply(surv_time, function(u){
      which.max(ev_time[ev_time<u])
    })
    # fitbz_dat<-fitbz[survPOS,]
    # bzH<-bz_hat[survPOS]
    fitbz_dat<-data.frame("hazard"=bhz_time%*%bz_hat,"time"=surv_time)
    bzH<-bz_hat[survPOS]
  }

  yy<-log(fitbz_dat$hazard)-log(bzH)
  sg_par<-as.numeric(coef(fit_s)) #as.numeric(coef(lm(yy~Xmat-1)))

  isl_lik<-sl_lik(par = sg_par,Xmat = Xmat,surv_time = surv_time,censor_status = censor_status,bz_hat = bzH,cum_bzhat = fitbz_dat[,1])

  iter<-1
  iter_dif<-NULL
  ret_obj<-list("CoxPH_obj"=list("hazard"=bz_hat,
                                 "cum_hazard"=fitbz,
                                 "coefficients"=sg_par),
                "cnvrgn" = iter_dif,
                "rc_thresh"=NULL)
  repeat{
    dif<-NULL
    # update of survival parameters
    expRW<-as.numeric(exp(Xmat%*%matrix(sg_par)))
    uexpRW<-hadamard.prod(Xmat,outer(expRW,rep(1,p)))

    if(constant_hazard){
      hzr<-sum(ev_n)/sum(surv_time*expRW)
      nbz_hat<-rep(hzr,length(ev_time))
      bz_dif<-sqrt(sum((nbz_hat-bz_hat)^2))/sqrt(sum((bz_hat)^2))
      bz_hat<-nbz_hat
      bzH<-rep(hzr,n)
      fitbz<-data.frame("hazard"=hzr*ev_time,"time"=ev_time)#basehaz(fit_s)
      fitbz_dat<-data.frame("hazard"=hzr*surv_time,"time"=surv_time)
    } else{
      nbz_hat<-sapply(2:length(ev_time), function(u){
        ev_n[u-1]/sum(expRW[rsk_ev[[u-1]]]*(pmin(surv_time[rsk_ev[[u-1]]],ev_time[u])-ev_time[u-1]))
      })

      bz_dif<-sqrt(sum((nbz_hat-bz_hat)^2))/sqrt(sum((bz_hat)^2))
      bz_hat<-nbz_hat
      bzH<-bz_hat[survPOS]
      fitbz<-data.frame("hazard"=cumsum(bz_hat*(ev_time[2:an]-ev_time[1:(an-1)])),"time"=ev_time[-1])#basehaz(fit_s)
      fitbz_dat<-data.frame("hazard"=bhz_time%*%bz_hat,"time"=surv_time)
    }

    # update of baseline hazard at failure time
    # nbz_hat<-sapply(seq_len(length(ev_time)),function(i){
    #   ev_n[i]/sum(expRW[rsk_ev[[i]]])
    # })
    # bz_dif<-sqrt(sum((nbz_hat-bz_hat)^2))/sqrt(sum((bz_hat)^2))
    # bz_hat<-nbz_hat
    # bzH<-bz_hat[survPOS]
    # fitbz<-data.frame("hazard"=cumsum(bz_hat),"time"=ev_time)#basehaz(fit_s)
    # fitbz_dat<-fitbz[survPOS,]

    # score function
    surv_scr<-hadamard.prod(Xmat,outer(event_status,rep(1,p)))-hadamard.prod(uexpRW,outer(fitbz_dat$hazard,rep(1,p))) #(uexpRW*fitbz_dat[,1])

    Ix<-Reduce('+',lapply(1:n,function(i){outer(surv_scr[i,],surv_scr[i,])}))

    nsg_par<-sg_par+as.numeric(solve(Ix,colSums(surv_scr)))

    sg_dif<-sqrt(sum((nsg_par-sg_par)^2))/sqrt(sum(sg_par^2))

    dif<-c(dif,bz_dif,sg_dif)

    sg_par<-nsg_par

    sllik<-sl_lik(par = sg_par,Xmat = Xmat,surv_time = surv_time,censor_status = censor_status,bz_hat = bzH,cum_bzhat = fitbz_dat[,1])

    ll_dif<-(sllik-isl_lik)/abs(isl_lik)

    iter_dif<-rbind(iter_dif,c(dif,sllik,ll_dif))

    if(cnv_type=="Param"){
      rc_thrsh<-max(dif)
    } else{
      rc_thrsh<-ll_dif
    }

    if(iter>max_iter | rc_thrsh<=cnv_crit)
      break
    iter<-iter+1
    isl_lik<-sllik
    ret_obj<-list("CoxPH_obj"=list("hazard"=bz_hat,
                                   "cum_hazard"=fitbz,
                                   "coefficients"=sg_par),
                  "cnvrgn" = iter_dif,
                  "rc_thresh"=rc_thrsh)
  }

  ret_obj
}

#' Provides similar summary as survival::summary.coxph
#'
#' @param object takes the fitted model object obtained by fitting supFTSVD_JM
#' @return it returns a data frame that contains estimated coefficient, hazard
#' ratio, and associated Wald-type test for the null hypothesis of zero
#' @example examples/example_summary_coxph_supFTSVD_JM.R
#' @export
summary_coxph_supFTSVD_JM<-function(object){
  cf_est<-as.numeric(object$CoxPH_obj$coefficients)
  cf_se<-as.numeric(diag(solve(object$CoxPH_obj$infor_mat)))
  res_dat<-data.frame("coef"=cf_est,
             "exp(coef)"=exp(cf_est),
             "se(coef)" = sqrt(cf_se),
             "z" = cf_est/sqrt(cf_se),
             "P-value(two_sided)"=2*pnorm(abs(cf_est/sqrt(cf_se)),mean = 0,sd = 1,lower.tail = FALSE))
  if(ncol(object$CoxPH_obj$surv_design)==ncol(object$A.hat)){
    rownames(res_dat)<-paste("Component",1:ncol(object$A.hat),sep = "_")
  } else{
    surv_cov<-ncol(object$CoxPH_obj$surv_design)-ncol(object$A.hat)
    rownames(res_dat)<-c(colnames(object$CoxPH_obj$surv_design)[1:surv_cov],
                         paste("Component",1:ncol(object$A.hat),sep = "_"))
  }
  res_dat
}

#' It splits data into g-fold and fits supFTSVD-JM, separately treating each
#' fold as test data.
#'
#' @param datlist a list with n elements, each is a matrix of dimension p+1
#' times m_i. The first row represents time points where measurements were
#' obtained. Next p rows represent the function values observed at the time
#' points
#' @param response a matrix of dimension n by q contains the design matrix for
#' supervision of the subject loadings
#' @param interval represents the domain of the function. The default is NULL.
#' @param r rank of the supFTSVD model
#' @param resolution grid size on which the singular functions will be estimated
#' @param CVPhi a logical scalar representing whether cross-validation (CV) will be
#' performed for determining the smoothness parameters involved with singular
#' functions
#' @param K number of folds for the CV step
#' @param cvT number of initial iterations at which CV will be performed
#' @param smooth Smoothing parameter for RKHS norm when CVPhi is FALSE. With
#' CVPhi=TRUE, a vector of numeric values for smoothness parameters for grid
#' search
#' @param maxiter Maximum number of iteration. Default: 20.
#' @param epsilon Convergence threshold for objective function value
#' @param KInd is a vector indices to construct the folds for CV. The default is
#' NULL that indicates that folds will be constructed randomly inside
#' @param rsvd_seed supFTSVD uses random svd for computational advantage while
#' setting inial values for feature loading vector. This seed number will be used
#' before setting the initial values. It confirms the reproducibility of your results.
#' @param conv_criteria the default is par which sets the threshold on the
#' difference between updated and current estimate at M step to stop the iteration.
#' It is computationally expensive; thus we recommend to use "cond_lik" instead of "par".
#' The result is not sensitive.
#' @param rank_fold number of folds to be used in the cross-validation
#' @param rfolds_seed seed number to be used in random splitting of the data
#' @return a list with following objects:
#' \itemize{
#'  \item fold_res : detailed results for different folds
#'  \item rank : rank of the fitted model
#' }
#' @example examples/example_g_fold_supFTSVD_JM.R
#' @export
gfold_supFTSVD_JM<-function(datlist, response, interval = NULL,
                             r = 1,resolution=100, CVPhi=FALSE, K=5, cvT=5,
                             smooth=1e-8,maxiter=20, epsilon=1e-4,KInd=NULL,
                             surv_time,censor_status,survX=NULL,mc_ss=1000,
                             rsvd_seed=100,conv_criteria="par",
                             rank_fold=5,rfolds_seed=50,stratum=NULL,
                             par_computing=TRUE,n_core=5,scale=TRUE,
                            constant_hazard=TRUE,change_points=NULL){
  n<-length(datlist)
  if(par_computing){
    library(doParallel)
    library(foreach)
    library(parallel)
    library(doSNOW)
    cls<-parallel::makeCluster(n_core)
    registerDoParallel(cls)
  }
  if(!is.null(stratum)){
    grp_frq<-as.numeric(table(stratum))
    lv_grp<-as.character(unique(stratum))
    set.seed(rfolds_seed)
    SJ_ind<-do.call(c,lapply(seq_len(length(lv_grp)),function(i){
      sample(1:grp_frq[i],grp_frq[i],replace=FALSE)
    }))
    rfolds_indx<-do.call(c,lapply(seq_len(length(lv_grp)),function(i){
      sj_ind<-sample(1:grp_frq[i],grp_frq[i],replace=FALSE)
      c(rep(1:rank_fold,each=grp_frq[i]%/%rank_fold),sample(1:rank_fold,grp_frq[i]%%rank_fold,replace = FALSE))
    }))
  } else{
    set.seed(rfolds_seed)
    SJ_ind<-sample(1:n,n,replace=FALSE)
    rfolds_indx<-c(rep(1:rank_fold,each=n%/%rank_fold),sample(1:rank_fold,n%%rank_fold,replace = FALSE))
  }
  if(par_computing){
    crank<-foreach(g=1:rank_fold,
                   .combine = "c",
                   .packages = c("rsvd","matrixcalc","survival","MASS"),
                   .export = c("supFTSVD_JM",
                               "predict.supFTSVD_JM",
                               "norm_llik_supFTSVD",
                               "sl_lik",
                               "bernoulli_kernel",
                               "cv_freg_rkhs",
                               "freg_rkhs"),
                   .errorhandling = "remove") %dopar% {
                     train_indx<-sort(SJ_ind[which(rfolds_indx!=g)])
                     test_indx<-sort(SJ_ind[which(rfolds_indx==g)])
                     train_dat<-datlist[train_indx]
                     test_dat<-datlist[test_indx]
                     fit_model<-supFTSVD_JM(datlist=train_dat, response=response[train_indx,], interval = interval,
                                            r = r,resolution=resolution, CVPhi=CVPhi, K=K, cvT=cvT,
                                            smooth=smooth,maxiter=maxiter, epsilon=epsilon,KInd=KInd,
                                            surv_time=surv_time[train_indx],
                                            censor_status=censor_status[train_indx],
                                            survX=survX[train_indx,],
                                            mc_ss=mc_ss,
                                            rsvd_seed=rsvd_seed,
                                            conv_criteria=conv_criteria,
                                            scale=scale,
                                            constant_hazard=constant_hazard,
                                            change_points=change_points)
                     predFIT<-predict.supFTSVD_JM(obj = fit_model,
                                                  designM = response[test_indx,],
                                                  new_dat = test_dat,
                                                  mc_ss = mc_ss,
                                                  surv_time=surv_time[test_indx],
                                                  censor_status=censor_status[test_indx],
                                                  survX=survX[test_indx,],
                                                  TimeG=NULL)

                     test_Xb<-response[test_indx,]%*%fit_model$Gamma
                     test_conU<-predFIT$subEFF-test_Xb

                     lng_lik_val<-norm_llik_supFTSVD(datlist = test_dat,
                                                     x_matrix = response[test_indx,],
                                                     comp_beta = fit_model$Gamma,
                                                     subj_dev = test_conU,
                                                     feat_loading = fit_model$B.hat,
                                                     sing_func = fit_model$Phi.hat,
                                                     sing_func_arg = fit_model$time,
                                                     sing_val = fit_model$Lambda,
                                                     sigR = fit_model$Sigma2R,
                                                     sig = fit_model$Sigma2)

                     # Quautities associated with survival model
                     sg_par<-as.numeric(fit_model$CoxPH_obj$coefficients)
                     fitbz<-fit_model$CoxPH_obj$cum_hazard
                     appbz<-approx(x=fitbz$time,y=fitbz$hazard,xout=surv_time[test_indx],rule=2)
                     fitbz_dat<-data.frame("hazard"=appbz$y,"time"=appbz$x)
                     bzH<-approx(x=fit_model$CoxPH_obj$cum_hazard$time, y=fit_model$CoxPH_obj$hazard,xout=surv_time[test_indx],rule = 2)$y

                     if(!is.null(survX)){
                       test_survX<-cbind(survX[test_indx,],predFIT$subEFF)
                     } else{
                       test_survX<-predFIT$subEFF
                     }

                     sl_lik_val<-sl_lik(par=sg_par,
                                        Xmat=test_survX,
                                        surv_time=surv_time[test_indx],
                                        censor_status=censor_status[test_indx],
                                        bz_hat=bzH,
                                        cum_bzhat=fitbz_dat[,1])

                     lng_lik_val + sl_lik_val
                   }
  } else{
    crank<-sapply(1:rank_fold, function(g){
      train_indx<-sort(SJ_ind[which(rfolds_indx!=g)])
      test_indx<-sort(SJ_ind[which(rfolds_indx==g)])
      train_dat<-datlist[train_indx]
      test_dat<-datlist[test_indx]
      fit_model<-supFTSVD_JM(datlist=train_dat, response=response[train_indx,], interval = interval,
                             r = r,resolution=resolution, CVPhi=CVPhi, K=K, cvT=cvT,
                             smooth=smooth,maxiter=maxiter, epsilon=epsilon,KInd=KInd,
                             surv_time=surv_time[train_indx],
                             censor_status=censor_status[train_indx],
                             survX=survX[train_indx,],
                             mc_ss=mc_ss,
                             rsvd_seed=rsvd_seed,
                             conv_criteria=conv_criteria,
                             scale=scale,
                             constant_hazard=constant_hazard,
                             change_points=change_points)
      predFIT<-predict.supFTSVD_JM(obj = fit_model,
                                   designM = response[test_indx,],
                                   new_dat = test_dat,
                                   mc_ss = mc_ss,
                                   surv_time=surv_time[test_indx],
                                   censor_status=censor_status[test_indx],
                                   survX=survX[test_indx,],
                                   TimeG=NULL)

      test_Xb<-response[test_indx,]%*%fit_model$Gamma
      test_conU<-predFIT$subEFF-test_Xb

      lng_lik_val<-norm_llik_supFTSVD(datlist = test_dat,
                                      x_matrix = response[test_indx,],
                                      comp_beta = fit_model$Gamma,
                                      subj_dev = test_conU,
                                      feat_loading = fit_model$B.hat,
                                      sing_func = fit_model$Phi.hat,
                                      sing_func_arg = fit_model$time,
                                      sing_val = fit_model$Lambda,
                                      sigR = fit_model$Sigma2R,
                                      sig = fit_model$Sigma2)

      # Quautities associated with survival model
      sg_par<-as.numeric(fit_model$CoxPH_obj$coefficients)
      fitbz<-fit_model$CoxPH_obj$cum_hazard
      appbz<-approx(x=fitbz$time,y=fitbz$hazard,xout=surv_time[test_indx],rule=2)
      fitbz_dat<-data.frame("hazard"=appbz$y,"time"=appbz$x)
      bzH<-approx(x=fit_model$CoxPH_obj$cum_hazard$time, y=fit_model$CoxPH_obj$hazard,xout=surv_time[test_indx],rule = 2)$y

      if(!is.null(survX)){
        test_survX<-cbind(survX[test_indx,],predFIT$subEFF)
      } else{
        test_survX<-predFIT$subEFF
      }

      sl_lik_val<-sl_lik(par=sg_par,
                         Xmat=test_survX,
                         surv_time=surv_time[test_indx],
                         censor_status=censor_status[test_indx],
                         bz_hat=bzH,
                         cum_bzhat=fitbz_dat[,1])

      lng_lik_val + sl_lik_val
    })
  }

  if(par_computing){
    rm(cls)
  }
  list("fold_res"=crank,"rank"=r)
}


#' It calculates R2, square of Pearson correlation coefficient, and R2 if a
#' linear model would have been fitted using the fitted components.
#'
#' @param obj fitted object obtained by supFTSVD_JM
#' @param obs_data observed data used in the fitted object
#' @export
r2_supFTSVD_JM<-function(obj,obs_data){
  n<-length(obs_data)
  r<-ncol(obj$A.hat)
  obs_time<-lapply(obs_data, function(u){u[1,]})
  newSF<-lapply(obs_time, function(nT){
    do.call(cbind,lapply(1:r, function(k){
      approx(x=obj$time,y=obj$Phi.hat[,k],xout = nT,rule = 2)$y
    }))
  })

  obs_y<-do.call(c,lapply(1:n,function(i){
    as.numeric(obs_data[[i]][-1,])
  }))

  tss<-sum((obs_y-mean(obs_y))^2)

  # R2 square of supFTSVD for longitudinal data
  # using Y_hat
  fit_y<-sapply(1:r, function(k){
    do.call(c,lapply(1:n,function(i){
      (obj$Lambda[k]*obj$A.hat[i,k])*as.numeric((outer(obj$B.hat[,k],newSF[[i]][,k])))
    }))
  })

  err_y<-sapply(1:r,function(k){
    if(k==1){
      sum((obs_y-fit_y[,1])^2)
    } else{
      sum((obs_y-rowSums(fit_y[,1:k]))^2)
    }
  })


  crr_y<-sapply(1:r,function(k){
    if(k==1){
      (cor(obs_y,fit_y[,1]))^2
    } else{
      (cor(obs_y,rowSums(fit_y[,1:k])))^2
    }
  })

  r2_y<-1-(err_y/tss)

  r2lm_y<-sapply(1:r,function(k){
    if(k==1){
      summary(lm(obs_y~fit_y[,1]-1))$r.square
    } else{
      summary(lm(obs_y~fit_y[,1:k]-1))$r.square
    }
  })


  # using Xb
  fit_xb<-sapply(1:r, function(k){
    do.call(c,lapply(1:n,function(i){
      (obj$Lambda[k]*obj$Xb[i,k])*as.numeric((outer(obj$B.hat[,k],newSF[[i]][,k])))
    }))
  })


  err_xb<-sapply(1:r,function(k){
    if(k==1){
      sum((obs_y-fit_xb[,1])^2)
    } else{
      sum((obs_y-rowSums(fit_xb[,1:k]))^2)
    }
  })

  r2_xb<-1-(err_xb/tss)

  crr_xb<-sapply(1:r,function(k){
    if(k==1){
      (cor(obs_y,fit_xb[,1]))^2
    } else{
      (cor(obs_y,rowSums(fit_xb[,1:k])))^2
    }
  })

  r2lm_xb<-sapply(1:r,function(k){
    if(k==1){
      summary(lm(obs_y~fit_xb[,1]-1))$r.square
    } else{
      summary(lm(obs_y~fit_xb[,1:k]-1))$r.square
    }
  })

  data.frame("R2_Yhat"=r2_y,
             "r_Yhat"=crr_y,
             "R2lm_Yhat"=r2lm_y,
             "R2_Xb"=r2_xb,
             "r_Xb"=crr_xb,
             "R2lm_Xb"=r2lm_xb)
}
