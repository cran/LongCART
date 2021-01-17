
PoS<- function(type, nsamples, null.value=NULL, alternative="greater", 
        N = NULL, D = NULL, a = 1,
        succ.crit = "trial", Z.crit.final = 1.96, alpha.final = 0.025,  clin.succ.threshold = NULL, 
        se.exp=NULL,
        meandiff.prior = NULL, mean.prior = NULL, sd.prior = NULL, propdiff.prior = NULL, prop.prior = NULL, hr.prior = NULL, D.prior = NULL) 
{
par(mar=c(1,1,1,1))
plot.new()
dev.off()

if(nsamples>2 | nsamples<1) stop("Number of samples should be either 1 or 2")
if(succ.crit=="clinical" & is.null(clin.succ.threshold)) stop("For clinical success, clin.succ.threshold must be specified.\n")
if(succ.crit=="trial" & is.null(Z.crit.final) & is.null(alpha.final)) stop("For trial success, Z.crit.final or alpha.final must be specified.\n")
if(type!="cont" & type!="bin" & type!="surv") stop("Endpoint type can be either cont (for continuous endpoint) or 
                                                    bin (for binary endpoint) or surv (for survival endpoint)\n")
if(type=="cont"){
	if(is.null(N) ) stop("For continuous endpoint, N must be specified.\n")
	if(is.null(sd.prior)) stop("For continuous case, sd.prior must be specified.\n")
	if(is.null(se.exp) ) stop("For continuous endpoint, se.exp must be specified.\n")
	if(nsamples==1 & is.null(mean.prior)) stop("For one sample continuous case, mean.prior must be specified.\n")
	if(nsamples==2 & is.null(meandiff.prior)) stop("For two samples continuous case, meandiff.prior must be specified.\n")
	}
if(type=="bin"){
	if(is.null(N) ) stop("For binary endpoint, N must be specified.\n")
	if(is.null(sd.prior)) stop("For binary case, sd.prior must be specified.\n")
	if(nsamples==2 & is.null(se.exp)) stop("For two sample binary case, se.exp must be specified.\n")
	if(nsamples==1 & is.null(prop.prior)) stop("For one sample binary case, prop.prior must be specified.\n")
	if(nsamples==2 & is.null(propdiff.prior)) stop("For two samples binary case, propdiff.prior must be specified.\n")
	}
if(type=="surv"){
	if(is.null(D) ) stop("For survival endpoint, D must be specified.\n")
	if(nsamples==1 ) stop("For survival endpoint, one sample case is not supported.\n")
	if(is.null(sd.prior) & is.null(D.prior)) stop("For survival case, either sd.prior or D.prior must be specified.\n")
	if(nsamples==2 & is.null(hr.prior)) stop("For survival endpoint, hr.prior must be specified.\n")
	}

if(is.null(alternative)) alternative=ifelse(type=="surv", "less", "greater")

Z1.crit=ifelse(!is.na(Z.crit.final), Z.crit.final, abs(qnorm(alpha.final)))

#--- Determine "r"
r=ifelse(nsamples==1, 1, (a+1)/sqrt(a))

#--- Determine "k.tilde"
if(type=="bin" & nsamples==1) se.exp=sqrt(prop.prior*(1-prop.prior)/N)
if(type=="surv" & nsamples==2) se.exp=r/sqrt(D)
k.tilde=se.exp

#--- Determine "sigma.0.sq"
if(type=="bin" & nsamples==1) se.exp=sqrt(prop.prior*(1-prop.prior)/N)
if(type=="surv" & nsamples==2) {
	sigma.0.sq=ifelse(!is.null(sd.prior), sd.prior^2, 4/D.prior)
      } else {
	sigma.0.sq=sd.prior^2
      }


#--- null.value
if(is.null(null.value)) null.value=ifelse(type=="surv",1,0)

#---------------------------------------------------------#
#    Determine "theta.min"                                #
#---------------------------------------------------------#
if(succ.crit=="clinical"){
	if(type=="surv"){
		theta.min=log(clin.succ.threshold) - log(null.value)
      } else theta.min=clin.succ.threshold - null.value
}

#---------------------------------------------------------#
#    Determine "theta.prior"                              #
#---------------------------------------------------------#

theta.prior<- NULL

#--- Continuous
if(type=="cont" & nsamples==1){ 
	parm.nam="mean"    
      theta.prior.raw=mean.prior 
	theta.prior=mean.prior - null.value
      } 
if(type=="cont" & nsamples==2){
            parm.nam="mean difference"
            theta.prior.raw=meandiff.prior     
		theta.prior=meandiff.prior - null.value
      }

#--- Binary endpoint
if(type=="bin" & nsamples==1){
      parm.nam="proportion"  
      theta.prior.raw=prop.prior  
	theta.prior=prop.prior - null.value
      } 
if(type=="bin" & nsamples==2){
      parm.nam="difference between two proportions"    
      theta.prior.raw=propdiff.prior
	theta.prior=propdiff.prior - null.value
      }

#--- Survival endpoint
if(type=="surv"){
      parm.nam="HR"    
      theta.prior.raw=log(hr.prior)
	theta.prior=log(hr.prior) - log(null.value)
    	} 

#---------------------------------------------------------#
# Adjust depending on direction on alternate hypothesis   #
#---------------------------------------------------------#

if(alternative=="less"){
   theta.prior=-theta.prior
   if(succ.crit=="clinical") theta.min=-theta.min
}


#---------------------------------------------------------#
#          Determine "gamma"                              #
#---------------------------------------------------------#

gamma=ifelse(succ.crit=="clinical", theta.min/k.tilde, Z1.crit)

#---------------------------------------------------------#
#              Printing text                              #
#---------------------------------------------------------#
alt.sign<- ifelse(alternative=="greater", ">", "<")
parm<- bquote(mu)
parm<- ifelse(type=="bin", bquote(pi), parm)
parm<- ifelse(type=="surv", "HR", parm)

if(nsamples==1| type=="surv"){
      parm.ext<- parm
      line1<- bquote("Test of hypothesis." ~ H[0] ~ ":" ~ .(parm) == .(null.value) ~ 
                               "vs." ~ H[1] ~ ":" ~ .(parm) ~ .(alt.sign) ~ .(null.value) )
	} else if(nsamples==2){
	parm.ext<- bquote(.(parm)[1] ~ "-" ~.(parm)[2])
	line1<- bquote("Test of hypothesis." ~ H[0] ~ ":" ~ .(parm.ext) == .(null.value) ~ 
                                     "vs." ~ H[1] ~ ":" ~ .(parm.ext) ~ .(alt.sign) ~ .(null.value) )
	}
plot.title<- bquote("Predictive distribution of " ~ .(parm.ext) ~ " given the prior distribution")

test<- "Approximate Z"
line1.1<- paste0("   (Statistical test: ", test, " test)")

prior.dist.val<- ifelse(type=='surv',
                   paste0('Log-normal (mean=log(', hr.prior, '), sd(log-scale)=', round(sqrt(sigma.0.sq), digits=4), ')'),
                   paste0('Normal (mean=', theta.prior+null.value, ', sd=', round(sqrt(sigma.0.sq), digits=4), ')'))
line3<- bquote("Prior distribution: " ~ .(parm.ext)  %~%  ~ .(prior.dist.val) ~ '' ) 


succ.nam<- ifelse(succ.crit == "trial", 
                  paste0("trial success (i.e., achieving statistical significance at 1-sided ", round(1-pnorm(Z1.crit), digits=4)," level)"),
                  paste0("clinical success (i.e., achieving threshold of ", clin.succ.threshold, ")"))
line4<- bquote("Evaluating for " ~ .(succ.nam))




#---------------------------------------------------------#
#       Calculate PoS                                     #
#---------------------------------------------------------#
sd.pred<- sqrt(sigma.0.sq + k.tilde^2)
quantile.pos= (theta.prior - k.tilde*gamma)/sd.pred

cat("Based on the prior distribution the probability of success is calculated as ", 
     round(pnorm(quantile.pos), digits=3), "\n") 

#---------------------------------------------------------#
#      Plotting                                           #
#---------------------------------------------------------#
mean.pred<- theta.prior.raw 
x <- seq(mean.pred-4*sd.pred, mean.pred+4*sd.pred, length=100)
hx <- dnorm(x, mean=mean.pred, sd=sd.pred)
if(type=="surv") x=exp(x)  #updated on 17-Jan-2021

par(mar=par()$mar + c(3,0,3,0), xpd=TRUE)
ymax=max(hx); ymin=min(hx); xmax=max(x); xmin=min(x)
plot(x, hx, type="l", lty=1, col=1, lwd=2, xlab="value", ylab="Density")
  
text(xmin , ymax +0.6*(ymax-ymin), bquote(.(line1) ~ .(line1.1)), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.5*(ymax-ymin), line3, cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.4*(ymax-ymin), line4, cex = 0.8, adj=c(0,0))

text(xmin, ymax +0.1*(ymax-ymin),plot.title, cex =0.9, adj=c(0,0))
}

