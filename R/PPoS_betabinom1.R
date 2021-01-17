
succ_ia_betabinom_one<- function( N, n, x,
        null.value=0, alternative="greater", test="z", correct=TRUE,
        succ.crit = "trial", Z.crit.final = 1.96, alpha.final = 0.025,  clin.succ.threshold = NULL, 
        a = 1, b=1) 
{
par(mar=c(1,1,1,1))
plot.new()
dev.off()

if(N<2 | n<1 | x<0 | N-n<1 | n<x) stop("Please correct the values of N, n, x.\n")

if(succ.crit=="clinical" & is.null(clin.succ.threshold)) stop("For clinical success, clin.succ.threshold must be specified.\n")
if(succ.crit=="trial" & is.null(Z.crit.final) & is.null(alpha.final)) stop("For trial success, Z.crit.final or alpha.final must be specified.\n")

#--- updated on Jan 17, 2021
Z1.crit=ifelse(!is.na(Z.crit.final), Z.crit.final, abs(qnorm(alpha.final)))
alpha<- 1-pnorm(Z1.crit)

# Beta Binomial Predictive distribution function
BetaBinom <- Vectorize(function(y, N, n, x, a,b){
	log.val <- lchoose(N-n, y) + lbeta(a+x+y,b+N-x-y) - lbeta(a+x,b+n-x)
	return(exp(log.val))
	})

#--- Treatment arm
y=c(0: (N - n))
y.predprob<- BetaBinom(y, N, n, x, a, b) 

#--- updated on Jan 17, 2021
pval<- NULL
if(succ.crit == "trial" & test=="binom"){
	for(y.i in 0:(N - n)) pval<- c(pval, binom.test(x=x+y.i, n=N, p=null.value, alternative = alternative)$p.value)
	pval<- cbind(y, y.predprob, pval)
	pval<- subset(pval, pval[,3]<alpha)
	} else if(succ.crit == "trial" & test=="z"){
	for(y.i in 0:(N - n)) pval<- c(pval, prop.test(x=x+y.i, n=N, p=null.value, alternative = alternative, correct=correct)$p.value)
	pval<- cbind(y, y.predprob, pval)
	pval<- subset(pval, pval[,3]<alpha)
	} else if(succ.crit == "clinical"){
	for(y.i in 0:(N - n)) pval<- c(pval, (x+y.i)/N)
	pval<- cbind(y, y.predprob, pval)
	pval<- subset(pval, pval[,3]>=clin.succ.threshold)
	}
ppos<- sum((pval[,2]))

#----- Summary: Start
alt.sign<- ifelse(alternative=="greater", ">", "<")
parm<- bquote(pi)
parm.nam="proportion"   
line1<- bquote("Test of hypothesis." ~ H[0] ~ ":" ~ .(parm) == .(null.value) ~ 
                               "vs." ~ H[1] ~ ":" ~ .(parm) ~ .(alt.sign) ~ .(null.value) )

test.nam<- ifelse(test=="exact", "Exact binomial", "Z")
line1.1<- paste0("   (Statistical test: ", test.nam, " test)")

ia.prop<- round(x/n, digits=3)
line2<- paste0("Interim results: Estimated proportion = ", ia.prop) 

prior.dist<- paste0("Beta(", a, ",", b, ")")
line3<- bquote("Prior distribution: " ~ pi  %~%  ~ .(prior.dist) ) 

succ.nam<- ifelse(succ.crit == "trial", 
                  paste0("trial success (i.e., achieving statistical significance at 1-sided ", round(1-pnorm(Z1.crit), digits=4)," level)"),
                  paste0("clinical success (i.e., achieving threshold of ", clin.succ.threshold, ")"))
line4<- bquote("Evaluating for " ~ .(succ.nam))

#----- Summary: End

#---- Print result
cat("Incoroprating prior information to the interim data, the predictive probability of success is", round(ppos, digits=3), "\n") 

#---- Plot of predictive distribution
par(mar=par()$mar + c(3,0,3,0), xpd=TRUE)
ymax=max(y.predprob)
ymin=min(y.predprob)
xmax=max(y)
xmin=min(y)
plot(y, y.predprob, type="o", lwd=1, col=1, xlab="value", ylab="Density")
text(xmin , ymax +0.6*(ymax-ymin), bquote(.(line1) ~ .(line1.1)), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.5*(ymax-ymin), bquote(.(line2) ~ ''), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.4*(ymax-ymin), line3, cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.3*(ymax-ymin), line4, cex = 0.8, adj=c(0,0))

text(xmin, ymax +0.1*(ymax-ymin),"Predictive distribution of number of post-interim response", 
     cex =0.9, adj=c(0,0))
}

