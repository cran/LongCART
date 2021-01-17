
succ_ia<- function(type, nsamples, null.value=NULL, alternative=NULL, 
        N = NULL, n = NULL, D = NULL, d = NULL, a = 1,
        meandiff.ia = NULL, mean.ia = NULL, propdiff.ia = NULL, prop.ia = NULL, hr.ia = NULL, stderr.ia = NULL, sd.ia = NULL, 
        succ.crit = "trial", Z.crit.final = 1.96, alpha.final = 0.025,  clin.succ.threshold = NULL, 
        meandiff.exp = NULL, mean.exp = NULL, propdiff.exp = NULL, prop.exp = NULL, hr.exp = NULL, 
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
	if(is.null(N) | is.null(n)) stop("For continuous endpoint, both N and n must be specified.\n")
	if(is.null(stderr.ia) & is.null(sd.ia)) stop("For continuous endpoint, either stderr.ia or sd.ia must be specified.\n")
	if(is.null(meandiff.ia) & is.null(mean.ia)) stop("For continuous endpoint, either meandiff.ia or mean.ia must be specified.\n")
	if(nsamples==1 & is.null(mean.ia)) stop("For one sample continuous case, mean.ia must be specified.\n")
	if(nsamples==2 & is.null(meandiff.ia)) stop("For two samples continuous case, meandiff.ia must be specified.\n")
	}
if(type=="bin"){
	if(is.null(N) | is.null(n)) stop("For binary endpoint, both N and n must be specified.\n")
	if(nsamples==2 & is.null(stderr.ia)) stop("For two sample binary case, stderr.ia must be specified.\n")
	if(is.null(propdiff.ia) & is.null(prop.ia)) stop("For binary endpoint, either propdiff.ia or prop.ia must be specified.\n")
	if(nsamples==1 & is.null(prop.ia)) stop("For one sample binary case, prop.ia must be specified.\n")
	if(nsamples==2 & is.null(propdiff.ia)) stop("For two samples binary case, propdiff.ia must be specified.\n")
	}
if(type=="surv"){
	if(is.null(D) | is.null(d)) stop("For survival endpoint, both D and d must be specified.\n")
	if(is.null(hr.ia)) stop("For survival endpoint, hr.ia must be specified.\n")
	if(nsamples==1 ) stop("For survival endpoint, one sample case is not supported.\n")
	if(nsamples==2 & is.null(hr.ia)) stop("For survival endpoint, hr.ia must be specified.\n")
	}

#--- updated on Jan 17, 2021
if(is.null(alternative)) alternative=ifelse(type=="surv", "less", "greater")

#--- updated on Jan 17, 2021
Z1.crit=ifelse(!is.na(Z.crit.final), Z.crit.final, abs(qnorm(alpha.final)))

#--- Determine "r"
r=ifelse(nsamples==1, 1, (a+1)/sqrt(a))

#--- Determine "t"
t=ifelse(type=="surv", d/D, n/N)

#--- Determine "k.ia" #--- Updated on Jan 17, 2021
if(type=="cont") stderr.ia=ifelse(!is.null(stderr.ia), stderr.ia, r*sd.ia/sqrt(n))
if(type=="bin" & nsamples==1) stderr.ia=sqrt(prop.ia*(1-prop.ia)/n)
if(type=="surv" & nsamples==2) stderr.ia=r/sqrt(d)
k.ia=sqrt(t)*stderr.ia

#--- null.value
if(is.null(null.value)) null.value=ifelse(type=="surv",1,0)

#--- Determine "theta.min"
if(succ.crit=="clinical"){
	if(type=="surv"){
		theta.min=log(clin.succ.threshold) - log(null.value)
      } else theta.min=clin.succ.threshold - null.value
}

#---------------------------------------------------------#
#    Determine "theta.est", "theta.exp", "theta.prior"    #
#---------------------------------------------------------#

psi=NULL
theta.exp<- NULL

#--- Continuous
if(type=="cont"){
	if(nsamples==1){ 
		parm.nam="mean"               #Added on 07-Jan-2021
            theta.est.raw=mean.ia         #Updated on 17-Jan-2021
		theta.est=mean.ia - null.value
      	if(!is.null(mean.exp)) theta.exp=mean.exp-null.value
      } else if(nsamples==2){
            parm.nam="mean difference"    #Added on 07-Jan-2021
            theta.est.raw=meandiff.ia     #Updated on 17-Jan-2021
		theta.est=meandiff.ia - null.value
      	if(!is.null(meandiff.exp)) theta.exp=meandiff.exp-null.value
      }
	if(nsamples==1 & !is.null(mean.prior) & !is.null(sd.prior) ){
            theta.prior.raw=mean.prior  #Updated on 17-Jan-2021
   		theta.prior=mean.prior-null.value; psi=sd.prior^2/(sd.prior^2+(k.ia^2/t)) 
    	} else if(nsamples==2 & !is.null(meandiff.prior) & !is.null(sd.prior)){
            theta.prior.raw=meandiff.prior #Updated on 17-Jan-2021
   		theta.prior=meandiff.prior-null.value; psi=sd.prior^2/(sd.prior^2+(k.ia^2/t)) 
      }
}

#--- Binary endpoint
if(type=="bin"){
	if(nsamples==1){ 
            parm.nam="proportion"    #Added on 07-Jan-2021
            theta.est.raw=prop.ia    #Updated on 17-Jan-2021
		theta.est=prop.ia - null.value
      	if(!is.null(prop.exp)) theta.exp=prop.exp-null.value
      } else if(nsamples==2){
            parm.nam="difference between two proportions"    #Added on 07-Jan-2021
            theta.est.raw=propdiff.ia    #Updated on 17-Jan-2021
		theta.est=propdiff.ia - null.value
      	if(!is.null(propdiff.exp)) theta.exp=propdiff.exp-null.value
      }
	if(nsamples==1 & !is.null(prop.prior) & !is.null(sd.prior) ){
            theta.prior.raw=prop.prior  #Updated on 17-Jan-2021
   		theta.prior=prop.prior-null.value; psi=sd.prior^2/(sd.prior^2+(k.ia^2/t)) 
    	} else if(nsamples==2 & !is.null(propdiff.prior) & !is.null(sd.prior)){
            theta.prior.raw=propdiff.prior  #Updated on 17-Jan-2021
   		theta.prior=propdiff.prior-null.value; psi=sd.prior^2/(sd.prior^2+(k.ia^2/t)) 
      }
}


#--- Survival endpoint
if(type=="surv"){
      parm.nam="HR"    #Added on 07-Jan-2021
      theta.est.raw=log(hr.ia)    #Updated on 17-Jan-2021
	theta.est=log(hr.ia) - log(null.value)
      if(!is.null(hr.exp)) theta.exp=log(hr.exp)-log(null.value)
	if(!is.null(hr.prior) & (!is.null(D.prior)|!is.null(sd.prior)) ){
            theta.prior.raw=log(hr.prior)  #Updated on 17-Jan-2021
   	 	theta.prior=log(hr.prior)-log(null.value)
            sd.prior=ifelse(!is.null(sd.prior), sd.prior, 2/sqrt(D.prior))
            psi=sd.prior^2/(sd.prior^2+(k.ia^2/t)) 
    	} 
}

#--- Adjust depending on direction on alternate hypothesis
if(alternative=="less"){
   theta.est=-theta.est
   if(succ.crit=="clinical") theta.min=-theta.min
   if(!is.null(theta.exp)) theta.exp=-theta.exp
   if(!is.null(psi)) theta.prior=-theta.prior
}



#----- Added on 07-Jan-2021: Start
alt.sign<- ifelse(alternative=="greater", ">", "<")
parm<- bquote(mu)
parm<- ifelse(type=="bin", bquote(pi), parm)
parm<- ifelse(type=="surv", "HR", parm)

#--- updated on 17-Jan-2021
if(type=="surv"|nsamples==1){
      parm.ext<- parm
	} else if (nsamples==2){
	parm.ext<- bquote(.(parm)[1] ~ "-" ~.(parm)[2])
	}
line1<- bquote("Test of hypothesis." ~ H[0] ~ ":" ~ .(parm.ext) == .(null.value) ~ 
               "vs." ~ H[1] ~ ":" ~ .(parm.ext) ~ .(alt.sign) ~ .(null.value) )
plot.title<- bquote("Predictive distribution of " ~ .(parm.ext) ~ " given the interim results")

test<- "Approximate Z"
line1.1<- paste0("   (Statistical test: ", test, " test)")

ia.val<- ifelse(type=='surv', hr.ia, theta.est+null.value)
se.nam<- ifelse(type=='surv', "standard error (log HR)", "standard error")
line2<- bquote("Interim results: Estimated " ~ .(parm.nam)  ==  .(ia.val) ~ "; " ~ .(se.nam) ==  .(round(stderr.ia, digits=4))) 
line2.1<- paste0('[post-interim trend: similar to interim]')
if(!is.null(theta.exp)){
	exp.val<- ifelse(type=='surv', hr.exp, theta.exp+null.value)
	line2.1<- bquote("[post-interim trend: expected " ~ .(parm.nam)  ==  .(exp.val) ~ "]") 
	}

line3<- bquote('Prior distribution: None' ~ '')
if(!is.null(psi)){
prior.dist.val<- ifelse(type=='surv',
                   paste0('Log-normal (mean=log(', hr.prior, '), sd(log-scale)=', round(sd.prior, digits=4), ')'),
                   paste0('Normal (mean=', theta.prior+null.value, ', sd=', round(sd.prior, digits=4), ')'))
line3<- bquote("Prior distribution: " ~ .(parm.ext)  %~%  ~ .(prior.dist.val) ~ '(' ~ psi == .(round(psi, digits=4)) ~ ')' ) 
}

succ.nam<- ifelse(succ.crit == "trial", 
                  paste0("trial success (i.e., achieving statistical significance at 1-sided ", round(1-pnorm(Z1.crit), digits=4)," level)"),
                  paste0("clinical success (i.e., achieving threshold of ", clin.succ.threshold, ")"))
line4<- bquote("Evaluating for " ~ .(succ.nam))

#----- Added on 07-Jan-2021: End


#--- Determine "gamma"
gamma=ifelse(succ.crit=="clinical", theta.min/k.ia, Z1.crit)

quantile.cp.est= ((theta.est/k.ia) - gamma)/sqrt(1-t)

cat("With the interim trend for post interim data, the condintional power is", 
     round(pnorm(quantile.cp.est), digits=3), "\n") 

#--- modified on Jan 16 - start
if(!is.null(theta.exp)){
    quantile.cp= ((theta.est*t+(1-t)*theta.exp)/k.ia - gamma)/sqrt(1-t)
    cat("With the specified trend for post interim data, the condintional power is", 
     round(pnorm(quantile.cp), digits=3), "\n") }
#--- modified on Jan 16 - end

cat("With the interim data, the predictive probability of success is", 
     round(pnorm(quantile.cp.est*sqrt(t)), digits=3), "\n") 

mean.pred1<- theta.est.raw       #updated on 17-Jan-2021
sd.pred1<- k.ia^2*(1/(1-t) + 1/t)
x1 <- seq(mean.pred1-4*sd.pred1, mean.pred1+4*sd.pred1, length=100)
hx1 <- dnorm(x1, mean=mean.pred1, sd=sd.pred1)
if(type=="surv") x1=exp(x1)  #updated on 17-Jan-2021

if(!is.null(psi)){ 
    quantile.ppos.prior= sqrt(t/(1-t))*
                   (1/sqrt((1-psi)*t + psi))*
                   ((theta.est/k.ia)*(t*(1-psi)+psi) + (1-t)*(1-psi)*(theta.prior/k.ia) - gamma)
    cat("psi=", psi, "\n")
    cat("Incoroprating prior information to the interim data, the predictive probability of success is", 
     round(pnorm(quantile.ppos.prior), digits=3), "\n") 

	mean.pred2<- psi*theta.est.raw + (1-psi)*theta.prior.raw    #updated on 17-Jan-2021
	sd.pred2<- k.ia^2*(1/(1-t) + psi/t)
	x2 <- seq(mean.pred2-4*sd.pred2, mean.pred2+4*sd.pred2, length=100)
	hx2 <- dnorm(x2, mean=mean.pred2, sd=sd.pred2)
      if(type=="surv") x2=exp(x2)      #updated on 17-Jan-2021
      }


par(mar=par()$mar + c(3,0,3,0), xpd=TRUE)
if(!is.null(psi)){
      ymax=max(hx1, hx2)
      ymin=min(hx1, hx2)
      xmax=max(x1, x2)
      xmin=min(x1, x2)
	plot(x1, hx1, type="l", lty=1, col=1, lwd=2,
           xlim=range(x1, x2), ylim=range(hx1, hx2),
           xlab="value", ylab="Density")
	points(x2, hx2, type="l", lty=2, col=2, lwd=2)
      legend(x="bottom", horiz=TRUE, inset=-0.5, bty = "n", cex=0.9,
             legend=c("Without prior distribution", "With prior distribution"),  lty=1:2, col=1:2, lwd=2)
     } else {
      ymax=max(hx1)
      ymin=min(hx1)
      xmax=max(x1)
      xmin=min(x1)
	plot(x1, hx1, type="l", lty=1, col=1, lwd=2, 
           xlab="value", ylab="Density")
      legend(x="bottom", horiz=TRUE, inset=-0.5, bty = "n", cex=0.9,
             legend=c("Without prior distribution"),  lty=1, col=1, lwd=2)
    }
text(xmin , ymax +0.6*(ymax-ymin), bquote(.(line1) ~ .(line1.1)), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.5*(ymax-ymin), bquote(.(line2) ~ .(line2.1)), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.4*(ymax-ymin), line3, cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.3*(ymax-ymin), line4, cex = 0.8, adj=c(0,0))

text(xmin, ymax +0.1*(ymax-ymin),plot.title, cex =0.9, adj=c(0,0))
}


