StabCont.surv<- function(data, timevar, censorvar, splitvar, 
                          time.dist="exponential", cens.dist="NA", event.ind=1, print=FALSE){

CDF.D <- function(x)
	{
	k <- seq(1:20)
	1 + 2 * sum((-1)^k * exp(- 2 * k^2 * x^2))
	}

  #------------Checks
  if(!exists(as.character(substitute(data)), envir=sys.frame(-1L))) stop("Dataset does not exist\n")
  if(!is.data.frame(data)) stop("Dataset does not exist\n")

  #-------- Checks for timevar variable
  if(!timevar %in% colnames(data))  stop("The column ", timevar, " containing follow-up time is missing in    dataset.\n")
  data$timevar<- data[[timevar]]

  #-------- Checks for censorvar variable
  if(!censorvar %in% colnames(data))  stop("The column ", censorvar, " containing censoring status is missing in    dataset.\n")
  data$censorvar<- I(data[[censorvar]]==event.ind)*1

  #--------- Check whether splitvar exist or not
  if(!splitvar %in% colnames(data)) stop("The column ", splitvar, " is missing in dataset.\n")
  data$splitvar<- data[[splitvar]]

  #--------- Check time to event distribution name
  if(I(time.dist!="exponential") & I(time.dist!="weibull") & 
     I(time.dist!="lognormal") & I(time.dist!="gaussian")) stop("Please check value for time.dist\n")

  #--------- Check censoring distribution name
  if(I(cens.dist!="exponential") & I(cens.dist!="weibull") & I(cens.dist!="NA") &
     I(cens.dist!="lognormal") & I(cens.dist!="gaussian")) stop("Please check value for cens.dist\n")

  #-------------Drop observations with missing timevar, censorvar & splitvar
  data.t<- data[!is.na(data[[timevar]])& !is.na(data[[censorvar]])& !is.na(data[[splitvar]]),]


message("Stability Test for Continuous grouping variable \n", appendLF=FALSE)

T<- data.t$timevar
delta<- data.t$censorvar

D<- sum(delta) #number of events
N=nrow(data.t) #number of subjects
G<- length(unique(data.t$splitvar)) #number of unique groups
Total.time<- sum(T)
ncensor=N-D #number of censored

if(print){
if(D>0){
	km.obj1=survfit(Surv(T, delta == 1) ~ 1, conf.type = "log-log")
	median1.obs=summary(km.obj1)$table["median"]
      }
if(ncensor>0){
km.obj2=survfit(Surv(T, delta == 0) ~ 1, conf.type = "log-log")
median2.obs=summary(km.obj2)$table["median"]
     }
}

if(D>0 & time.dist=="exponential"){
   lambda1=D/Total.time
   median1<- (log(2)/lambda1)
   if(print) message("\nEvent time: Estimated lambda=", round(lambda1, 5), 
                     ",  Estimated Median=", round(median1, 3), 
                     "  (KM estimate=", round(median1.obs, 3), ")\n", appendLF=FALSE)
}else
if(D>0 & time.dist=="weibull"){
   fit1<- survreg(Surv(timevar, censorvar) ~ 1, data = data.t, dist = time.dist)
   a1<- 1/fit1$scale  #shape in rweibull()
   b1<- exp(fit1$coeff) #scale in rweibull()
   alpha1<- a1
   lambda1<- (1/b1)^a1
   median1<- (log(2)/lambda1)^(1/alpha1)
   if(print) message("\nEvent time: Estimated lambda=", round(lambda1, 5),
                     ",  Estimated alpha=", round(alpha1, 5), 
                     ",  Estimated Median=", round(median1, 3), 
                     "  (KM estimate=", round(median1.obs, 3), ")\n", appendLF=FALSE)
}
if(D>0 & (time.dist=="lognormal"| time.dist=="gaussian")){
   fit1<- survreg(Surv(timevar, censorvar) ~ 1, data = data.t, dist = time.dist)
   mu1<- fit1$coeff  #mulog
   sd1<- fit1$scale #sdlog
   if(time.dist=="lognormal") y1<- (log(T)-mu1)/sd1
   if(time.dist=="gaussian") y1<- (T-mu1)/sd1
   h1<- dnorm(y1)/pnorm(-y1)
   if(time.dist=="lognormal") median1<- exp(mu1)
   if(time.dist=="gaussian") median1<- mu1 
   if(print) message("\nEvent time: Estimated mu=", round(mu1, 3), ", Estimated sd=", round(sd1, 3), 
                     ",  Estimated Median=", round(median1, 3), 
                     "  (KM estimate=", round(median1.obs, 3), ")\n", appendLF=FALSE)
}


if(ncensor>0 & cens.dist=="exponential"){
   lambda2=(N-D)/Total.time
   median2<- (log(2)/lambda2)
   if(print) message("Censoring time: Estimated lambda=", round(lambda2, 5), 
                     ",  Estimated Median=", round(median2, 3), 
                     "  (KM estimate=", round(median2.obs, 3), ")\n", appendLF=FALSE)
}else
if(ncensor>0 & cens.dist=="weibull"){
   fit2<- survreg(Surv(timevar, 1-censorvar) ~ 1, data = data.t, dist = cens.dist)
   a2<- 1/fit2$scale  #shape in rweibull()
   b2<- exp(fit2$coeff) #scale in rweibull()
   alpha2<- a2
   lambda2<- (1/b2)^a2
   median2<- (log(2)/lambda2)^(1/alpha2)
   if(print) message("Censoring time: Estimated lambda=", round(lambda2, 5),
                     ",  Estimated alpha=", round(alpha2, 5), 
                     ",  Estimated Median=", round(median2, 3), 
                     "  (KM estimate=", round(median2.obs, 3), ")\n", appendLF=FALSE)
}else
if(ncensor>0 & (cens.dist=="lognormal"|cens.dist=="gaussian")){
   fit2<- survreg(Surv(timevar, 1-censorvar) ~ 1, data = data.t, dist = cens.dist)
   mu2<- fit2$coeff  #mulog
   sd2<- fit2$scale #sdlog
   if(cens.dist=="lognormal") y2<- (log(T)-mu2)/sd2
   if(cens.dist=="gaussian") y2<- (T-mu2)/sd2
   h2<- dnorm(y2)/pnorm(-y2)
   if(cens.dist=="lognormal") median2<- exp(mu2)
   if(cens.dist=="gaussian") median2<- mu2 
   if(print) message("Censoring time: Estimated mu=", round(mu2, 3), ", Estimated sd=", round(sd2, 3), 
                     ",  Estimated Median=", round(median2, 3),
                     "  (KM estimate=", round(median2.obs, 3), ")\n", appendLF=FALSE)
}


#-score function
p<- 0; type<- NULL; u<- NULL; J.sqrt.inv<- NULL

#---- Event to time part
if(D>0){   #score function & J related to event time
   if(time.dist=="exponential"){
	u1<- matrix(delta/lambda1-T, ncol=1)         
      J11.sqrt.inv<- matrix(sqrt(N/D)*lambda1, 1,1)
      colnames(J11.sqrt.inv)<- rownames(J11.sqrt.inv)<- colnames(u1)<- c("lambda.T")
      p<- p+1; type<- c(type, 1); u<- cbind(u, u1)
   } else
   if(time.dist=="weibull"){
	u11<- delta/lambda1-T^alpha1  
	u12<- delta*(1/alpha1 + log(T))-lambda1*log(T)*(T^alpha1)
      u1<- cbind(u11, u12)
      J11.11<- (D/N)*lambda1^(-2)        
      J11.12<- mean((T^alpha1)*log(T))
      J11.22<- (D/N)*alpha1^(-2) + lambda1*mean((log(T))^2*(T^alpha1))     
      J11<- matrix(c(J11.11, J11.12, J11.12, J11.22), ncol=2, nrow=2)
      e.J11 <- eigen(J11)
      J11.sqrt<- e.J11$vectors %*% diag(sqrt(e.J11$values)) %*% t(e.J11$vectors)
      J11.sqrt.inv<- solve(J11.sqrt)
      colnames(J11.sqrt.inv)<- rownames(J11.sqrt.inv)<- colnames(u1)<- c("lambda.T", "alpha.T")
      p<- p+2; type<- c(type, 1, 1); u<- cbind(u, u1)
   } else
   if(time.dist=="lognormal"|time.dist=="gaussian"){
      u11<- (delta*y1 + (1-delta)*h1)/sd1 
	u12<- delta*(y1^2-1) + (1-delta)*y1*h1
      u1<- cbind(u11, u12)
      J11.inv<- N*fit1$var
      e1<- eigen(J11.inv)
      J11.sqrt.inv<- e1$vectors %*% diag(sqrt(e1$values)) %*% t(e1$vectors)
      colnames(J11.sqrt.inv)<- rownames(J11.sqrt.inv)<- colnames(u1)<- c("mu.T", "log(SD.T)")
      p<- p+2; type<- c(type, 1, 1); u<- cbind(u, u1)
   }
   if(!is.null(J.sqrt.inv)){
   J.sqrt.inv<- adiag(J.sqrt.inv,J11.sqrt.inv)
   }  else
   if(is.null(J.sqrt.inv)) J.sqrt.inv<- J11.sqrt.inv
}

#---- Censoring part
if(ncensor>0 & cens.dist!="NA"){   #score function & J related to event time
   if(cens.dist=="exponential"){
	u2<- matrix((1-delta)/lambda2-T, ncol=1)	    
      J22.sqrt.inv<- matrix(sqrt(N/(N-D))*lambda2,1,1)
      colnames(J22.sqrt.inv)<- rownames(J22.sqrt.inv)<- colnames(u2)<- c("lambda.C")
      p<- p+1; type<- c(type, 2); u<- cbind(u, u2)
   }
   if(cens.dist=="weibull"){
	u21<- (1-delta)/lambda2-T^alpha2  
	u22<- (1-delta)*(1/alpha2 + log(T))-lambda2*log(T)*(T^alpha2)
	u2<- cbind(u21, u22)
      J22.11<- (1-D/N)*lambda2^(-2)       
      J22.12<- mean((T^alpha2)*log(T))
      J22.22<- (1-D/N)*alpha2^(-2) + lambda2*mean((log(T))^2*(T^alpha2))  
      J22<- matrix(c(J22.11, J22.12, J22.12, J22.22), ncol=2, nrow=2)
      e.J22 <- eigen(J22)
      J22.sqrt<- e.J22$vectors %*% diag(sqrt(e.J22$values)) %*% t(e.J22$vectors)
      J22.sqrt.inv<- solve(J22.sqrt)
      colnames(J22.sqrt.inv)<- rownames(J22.sqrt.inv)<- colnames(u2)<-  c("lambda.C", "alpha.C")
      p<- p+2; type<- c(type, 2, 2); u<- cbind(u, u2)
   }
   if(cens.dist=="lognormal"|cens.dist=="gaussian"){
      u21<- ((1-delta)*y2 + delta*h2)/sd2 
	u22<- (1-delta)*(y2^2-1) + delta*y2*h2
	u2<- cbind(u21, u22)
      J22.inv<- N*fit2$var
      e2<- eigen(J22.inv)
      J22.sqrt.inv<- e2$vectors %*% diag(sqrt(e2$values)) %*% t(e2$vectors)
      colnames(J22.sqrt.inv)<- rownames(J22.sqrt.inv)<- colnames(u2)<- c("mu.C", "log(SD.C)")
      p<- p+2; type<- c(type, 2, 2); u<- cbind(u, u2)
   }
   if(!is.null(J.sqrt.inv)){
      J.sqrt.inv<- adiag(J.sqrt.inv,J22.sqrt.inv)
   }  else
   if(is.null(J.sqrt.inv)) J.sqrt.inv<- J22.sqrt.inv
}

if(print){
message("\nCheck: Sum of score function should be 0\n", appendLF=FALSE)
message(paste0(apply(u, 2, sum), " "), "\n", appendLF=FALSE)
message("\nMaximum of absolute value of score function\n", appendLF=FALSE)
message(paste0(apply(abs(u), 2, max), " "), "\n", appendLF=FALSE)
message("\nJ^(-1/2)\n", appendLF=FALSE)
message(paste0(capture.output(J.sqrt.inv), collapse = "\n"), "\n\n", appendLF=FALSE)
}

u.grp<- W<- W.thetahat<- matrix(NA, nrow=G, ncol=p)
for(pi in 1:p) {
	u.grp[,pi]<- tapply(u[,pi], data.t$splitvar, sum)
	W.thetahat[,pi]<- cumsum(u.grp[,pi])/sqrt(N)
	}
M.thetahat<- J.sqrt.inv%*%t(W.thetahat)

#Test statistic & p-value
Test.stat<- apply(abs(M.thetahat),1,max)
out.p<- 1- unlist(lapply(Test.stat, FUN=CDF.D))
out.p.adj<- p.adjust(out.p, method = "hochberg")

message('Test.statistic=', paste(" ", round(Test.stat, digits=3), sep=''), 
        ',    Adj. p-value=', paste(" ", round(out.p.adj, digits=3), sep=''), ' \n\n',  appendLF=FALSE)
if(which.max(Test.stat)==1) message('Greater evidence of heterogeneity in time-to-event distribution \n',   appendLF=FALSE)
if(which.max(Test.stat)==2) message('Greater evidence of heterogeneity in censoring distribution \n',  appendLF=TRUE)
ret<- list(pval=min(out.p.adj), type=type[which.max(Test.stat)]) 
ret
}


#StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
#               time.dist="exponential", cens.dist="exponential", event.ind=2, print=TRUE)

#StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
#               time.dist="weibull", cens.dist="weibull", event.ind=2, print=TRUE)

#StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
#               time.dist="lognormal", cens.dist="lognormal", event.ind=2, print=TRUE)

#StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
#               time.dist="gaussian", cens.dist="gaussian", event.ind=2, print=TRUE)

#StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
#               time.dist="gaussian", cens.dist="NA", event.ind=2, print=TRUE)

#StabCont.surv(data=lung, timevar="time", censorvar="status", splitvar="age", 
#               time.dist="gaussian", cens.dist="NA", event.ind=2, print=TRUE)




