
succ_ia_betabinom_two<- function( N.trt, N.con, 
        n.trt, x.trt, n.con, x.con,
        alternative="greater", test="z", 
        succ.crit = "trial", Z.crit.final = 1.96, alpha.final = 0.025,  clin.succ.threshold = NULL, 
        a.trt = 1, b.trt=1, a.con=1, b.con=1) 
{
par(mar=c(1,1,1,1))
plot.new()
dev.off()

if(N.trt<2 | N.con<2 | n.trt<1 | n.con<1 | x.trt<0 | x.con<0 | N.trt-n.trt<1 | N.con-n.con<1 | n.trt<x.trt| n.con<x.con) stop("Please correct the values of N.trt, N.con, n.trt, x.trt, n.con, x.con.\n")
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
y.trt=c(0: (N.trt - n.trt))
y.predprob.trt<- BetaBinom(y.trt, N.trt, n.trt, x.trt, a.trt, b.trt) 

#--- Control arm
y.con=c(0:(N.con - n.con))
y.predprob.con<- BetaBinom(y.con, N.con, n.con, x.con, a.con, b.con) 

pval<- NULL

if(succ.crit == "trial" & test=="fisher"){
      #print("Performing Fisher exact test")
	for(y.trt.i in 0:(N.trt - n.trt)){
		for(y.con.i in 0:(N.con - n.con)){
			dat<- matrix(c(x.trt+y.trt.i, N.trt-x.trt-y.trt.i, x.con+y.con.i, N.con-x.con-y.con.i), nrow=2)
			pval<- rbind(pval, c(y.trt.i, y.con.i, fisher.test(dat, alternative = alternative)$p.value))
			}
   		}
	pval<- subset(pval, pval[,3]<alpha)
	} else if(succ.crit == "trial" & test=="z"){
      #print("Performing z test")
	for(y.trt.i in 0:(N.trt - n.trt)){
		for(y.con.i in 0:(N.con - n.con)){
			pval<- rbind(pval, c(y.trt.i, y.con.i, prop.test(x=c((x.trt+y.trt.i), (x.con+y.con.i)), 
                                                                   n=c(N.trt, N.con), 
                                                                   alternative = alternative)$p.value))
			}
   		}
	pval<- subset(pval, pval[,3]<alpha)
	} else if(succ.crit == "clinical"){
	for(y.trt.i in 0:(N.trt - n.trt)){
		for(y.con.i in 0:(N.con - n.con)){
                  diff.prop<- (x.trt+y.trt.i)/N.trt - (x.con+y.con.i)/N.con
			pval<- rbind(pval, c(y.trt.i, y.con.i, diff.prop))
			}
   		}
	pval<- subset(pval, pval[,3]>=clin.succ.threshold)
	}

if(nrow(pval)==0) ppos=0
if(nrow(pval)>=1) ppos<- sum(BetaBinom(pval[,1], N.trt, n.trt, x.trt, a.trt, b.trt)*
           BetaBinom(pval[,2], N.con, n.con, x.con, a.con, b.con))


#----- Summary: Start
alt.sign<- ifelse(alternative=="greater", ">", "<")
parm<- bquote(pi)
parm.nam="difference between two proportions"   
parm.ext<- bquote(.(parm)[T] ~ "-" ~.(parm)[C])
line1<- bquote("Test of hypothesis." ~ H[0] ~ ":" ~ .(parm.ext) == 0 ~ 
                                     "vs." ~ H[1] ~ ":" ~ .(parm.ext) ~ .(alt.sign) ~ 0 )

test.nam<- ifelse(test=="fisher", "Fisher's exact", "Z")
line1.1<- paste0("   (Statistical test: ", test.nam, " test)")

ia.prop.trt<- round(x.trt/n.trt, digits=3)
ia.prop.con<- round(x.con/n.con, digits=3)
line2<- paste0("Interim results: Estimated proportion (treatment arm vs. control arm) = ", ia.prop.trt, " vs. ", ia.prop.con) 

prior.dist.trt<- paste0("Beta(", a.trt, ",", b.trt, ")")
prior.dist.con<- paste0("Beta(", a.con, ",", b.con, ")")
line3<- bquote("Prior distribution: " ~ pi[T]  %~%  ~ .(prior.dist.trt) ~ ";  " ~ pi[C]  %~%  ~ .(prior.dist.con)) 

succ.nam<- ifelse(succ.crit == "trial", 
                  paste0("trial success (i.e., achieving statistical significance at 1-sided ", round(1-pnorm(Z1.crit), digits=4)," level)"),
                  paste0("clinical success (i.e., achieving threshold of ", clin.succ.threshold, ")"))
line4<- bquote("Evaluating for " ~ .(succ.nam))

#----- Summary: End

#---- Print result
cat("Incoroprating prior information to the interim data, the predictive probability of success is", round(ppos, digits=3), "\n") 

#---- Plot of predictive distribution
par(mar=par()$mar + c(3,0,3,0), xpd=TRUE)
ymax=max(y.predprob.trt, y.predprob.con)
ymin=min(y.predprob.trt, y.predprob.con)
xmax=max(y.trt, y.con)
xmin=min(y.trt, y.con)
plot(y.trt, y.predprob.trt, type="o", lwd=2, col=1, lty=1, 
     xlim=c(xmin, xmax), ylim=c(ymin, ymax), xlab="value", ylab="Density")
points(y.con, y.predprob.con, type="o", lwd=2, col=2, lty=2)
text(xmin , ymax +0.6*(ymax-ymin), bquote(.(line1) ~ .(line1.1)), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.5*(ymax-ymin), bquote(.(line2) ~ ''), cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.4*(ymax-ymin), line3, cex = 0.8, adj=c(0,0))
text(xmin , ymax +0.3*(ymax-ymin), line4, cex = 0.8, adj=c(0,0))

text(xmin, ymax +0.1*(ymax-ymin),"Predictive distribution of number of post-interim response", 
     cex =0.9, adj=c(0,0))

legend(x="bottom", horiz=TRUE, inset=-0.5, bty = "n", cex=0.9,
       legend=c('Treatment arm', 'Control arm'), lty=1:2, lwd=2, col=1:2)
}


