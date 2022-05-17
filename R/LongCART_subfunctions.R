#-------------------------------------------------------------------------
#  LongCART algorithm to construct regression tree with Longitudinal data
#  Date: 28-Sep-2016
#  Developed by: Madan G Kundu (Email: madan_g.kundu@novartis.com)
#-------------------------------------------------------------------------
##### Following checks need to be done
#1. Check: No missing value for patid and response variable
#2. Check: Grouping varaible should have only observations for each patid
#3. Check: Length of gvars and tgvars should be equal
#4. Check: The variables should exist
############## End of Checks


#--------------------------------------------------------#
#   Stability test for categorical grouping variable     #
#--------------------------------------------------------#

StabCat<- function(data, patid, fixed, splitvar)
{
  #------------Checks
  if(!exists(as.character(substitute(data)), envir=sys.frame(-1L))) stop("Dataset does not exist\n")
  if(!is.data.frame(data)) stop("Dataset does not exist\n")

  #-------- Checks for patid variable
  if(!patid %in% colnames(data))  stop("The column ", patid, " containing subjects id is missing in dataset.\n")
  data$patid<- data[[patid]]
  #-------------Drop observations with missing patid
  data<- data[!is.na(data[[patid]]),]

  #-------- Checks for Y variable
  Y.name<- as.character(attr(as.Formula(fixed), "lhs"))
  if(!Y.name %in% colnames(data))
    stop("The column ", Y.name, " containing subjects id is missing in dataset.\n")
  #-------------Drop observations with missing Y and patid
  data<- data[!is.na(data[[Y.name]]),]

  #--------- Check whether splitvar exist or not
  if(!splitvar %in% colnames(data))
    stop("The column ", splitvar, " is missing in dataset.\n")


  sig.stab<- NA
  message("Stability Test for Categorical grouping variable \n", appendLF=FALSE)
  data.t<- data[!is.na(splitvar),]
  splitvar<- data[[splitvar]]
  splitvar<- splitvar[!is.na(splitvar)]

  temp<- gsub(pattern=" ", replacement="", fixed, fixed=TRUE)[3]
  noint=grepl(pattern="-1", temp, fixed=TRUE)

  form<-Formula(fixed)
  Y.name<- as.character(attr(form, "lhs"))
  X.name1<- as.character(attr(form, "rhs"))
  X.name1<- gsub(pattern=" ", replacement="", x=X.name1)
  X.name1<- gsub(pattern="-1", replacement="", x=X.name1)
  X.name<- unlist(strsplit(X.name1, "+", fixed = TRUE))

  Y<- data.t[[Y.name]]
  p<- length(X.name)

  if(X.name1=='1') X<- rep(1, length(Y))
  if(X.name1!='1'){
    for(p.i in 1:p){
      if(p.i==1) X<-data.t[[X.name[p.i]]]
      if(p.i>1)  X<-cbind(X, data.t[[X.name[p.i]]])
    }
  }

  id<- data.t$patid
  Group<- tapply(splitvar, id, unique)
  Group<- Group[!is.na(Group)]
  G<- length(unique(Group))

  ## Fit mixed model
  fit<- NULL
  if (X.name1=='1') {
    p<- 1
    fit <- try(fit <- lme(fixed=Y~1, random=~1|id) , silent=TRUE)
  } else if (noint) {
    p<- p
    fit <- try(fit <- lme(fixed=Y~X-1, random=~1|id) , silent=TRUE)
  } else if (!noint) {
    p<- p+1
    fit <- try(fit <- lme(fixed=Y~X, random=~1|id) , silent=TRUE)
  }

  if(!is.null(fit)){
    #Get sd_int and sigma.e
    sd_int<- as.numeric(VarCorr(fit)[1,2])
    sd_err<- as.numeric(VarCorr(fit)[2,2])

    #Get beta.hat
    beta <- matrix(fit$coeff$fixed,ncol=1)

    #Get J
    idvec<- unique(id)
    N<- length(idvec)
    J=solve(vcov(fit))/N

    u<- matrix(0, N, p)

    for(i in 1:N)	{
      sub=idvec[i]
      Y.i<- subset(Y, id==sub)
      X.i<- subset(X, id==sub)
      if (noint|X.name1=='1') X.i<- model.matrix(~X.i-1)
      if (!noint) X.i<- model.matrix(~X.i)
      n.i<- length(Y.i)
      if(p==1) X.i<- matrix(X.i, nrow=n.i, ncol=1)
      if(n.i==1) V.i<- matrix(sd_int^2+sd_err^2, nrow=n.i, ncol=n.i)
      else V.i<- matrix(sd_int^2, nrow=n.i, ncol=n.i) + diag(rep(sd_err^2, n.i))
      u[i,]<- t(t(X.i)%*%solve(V.i)%*%(Y.i-X.i%*%beta))
    }

    u.grp<- W<- W.thetahat<- matrix(NA, nrow=G, ncol=p)
    n.grp<- table(splitvar)

    for(pi in 1:p) u.grp[,pi]<- tapply(u[,pi], Group, sum)

    #Test statistic
    Test.stat=0
    for(Gi in 1:G){
      txx<- matrix(u.grp[Gi,], ncol=1)
      Test.stat<- Test.stat+t(txx)%*%solve(n.grp[Gi]*J)%*%txx
    }

    #p value
    sig.stab<- 1- pchisq(Test.stat, (G-1)*p)

    message('Test.statistic=', round(c(Test.stat), digits=3), ',    p-value=', round(c(sig.stab),digits=3), ' \n', appendLF=FALSE)
  }
  ret<- list(pval=c(sig.stab))
  ret
}

#--------------------------------------------------------#
#   Stability test for continuous grouping variable      #
#--------------------------------------------------------#


StabCont<- function(data, patid, fixed, splitvar)
{
  CDF.D <- function(x)
  {
    k <- seq(1:20)
    1 + 2 * sum((-1)^k * exp(- 2 * k^2 * x^2))
  }

  #------------Checks
  if(!exists(as.character(substitute(data)), envir=sys.frame(-1L))) stop("Dataset does not exist\n")
  if(!is.data.frame(data)) stop("Dataset does not exist\n")

  #-------- Checks for patid variable
  if(!patid %in% colnames(data))  stop("The column ", patid, " containing subjects id is missing in dataset.\n")
  data$patid<- data[[patid]]
  #-------------Drop observations with missing patid
  data<- data[!is.na(data[[patid]]),]

  #-------- Checks for Y variable
  Y.name<- as.character(attr(as.Formula(fixed), "lhs"))
  if(!Y.name %in% colnames(data))
    stop("The column ", Y.name, " containing subjects id is missing in dataset.\n")
  #-------------Drop observations with missing Y and patid
  data<- data[!is.na(data[[Y.name]]),]

  #--------- Check whether splitvar exist or not
  if(!splitvar %in% colnames(data))
    stop("The column ", splitvar, " is missing in dataset.\n")

  sig.stab<- NA
  message("Stability Test for Continuous grouping variable \n", appendLF=FALSE)
  data.t<- data[!is.na(splitvar),]
  splitvar<- data[[splitvar]]
  splitvar<- splitvar[!is.na(splitvar)]

  temp<- gsub(pattern=" ", replacement="", fixed, fixed=TRUE)[3]
  noint=grepl(pattern="-1", temp, fixed=TRUE)

  form<-Formula(fixed)
  Y.name<- as.character(attr(form, "lhs"))
  X.name1<- as.character(attr(form, "rhs"))
  X.name1<- gsub(pattern=" ", replacement="", x=X.name1)
  X.name1<- gsub(pattern="-1", replacement="", x=X.name1)
  X.name<- unlist(strsplit(X.name1, "+", fixed = TRUE))

  Y<- data.t[[Y.name]]
  p<- length(X.name)

  if(X.name1=='1') X<- rep(1, length(Y))
  if(X.name1!='1'){
    for(p.i in 1:p){
      if(p.i==1) X<-data.t[[X.name[p.i]]]
      if(p.i>1)  X<-cbind(X, data.t[[X.name[p.i]]])
    }
  }

  id<- data.t$patid
  Group<- tapply(splitvar, id, unique)
  Group<- Group[!is.na(Group)]
  G<- length(unique(Group))

  ## Fit mixed model
  fit<- NULL
  if (X.name1=='1') {
    p<- 1
    fit <- try(fit <- lme(fixed=Y~1, random=~1|id) , silent=TRUE)
  } else if (noint) {
    p<- p
    fit <- try(fit <- lme(fixed=Y~X-1, random=~1|id) , silent=TRUE)
  } else if (!noint) {
    p<- p+1
    fit <- try(fit <- lme(fixed=Y~X, random=~1|id) , silent=TRUE)
  }

  if(!is.null(fit)){
    #Get sd_int and sigma.e
    sd_int<- as.numeric(VarCorr(fit)[1,2])
    sd_err<- as.numeric(VarCorr(fit)[2,2])

    #Get beta.hat
    beta <- matrix(fit$coeff$fixed,ncol=1)

    #Get J
    idvec<- unique(id)
    N<- length(idvec)
    J=solve(vcov(fit))/N

    if(p==1){
      J.sqrt<- sqrt(J)
    } else {
      J.eig <- eigen(J)
      J.sqrt <- J.eig$vectors %*% diag(sqrt(J.eig$values)) %*% solve(J.eig$vectors)
    }

    u<- matrix(0, N, p)

    for(i in 1:N)	{
      sub=idvec[i]
      Y.i<- subset(Y, id==sub)
      n.i<- length(Y.i)
      X.i<- subset(X, id==sub)
      if (noint|X.name1=='1') X.i<- model.matrix(~X.i-1)
      if (!noint) X.i<- model.matrix(~X.i)
      if(p==1) X.i<- matrix(X.i, nrow=n.i, ncol=1)
      if(n.i==1) V.i<- matrix(sd_int^2+sd_err^2, nrow=n.i, ncol=n.i)
      else V.i<- matrix(sd_int^2, nrow=n.i, ncol=n.i) + diag(rep(sd_err^2, n.i))
      u[i,]<- t(t(X.i)%*%solve(V.i)%*%(Y.i-X.i%*%beta))
    }

    u.grp<- W<- W.thetahat<- matrix(NA, nrow=G, ncol=p)
    n.grp<- table(splitvar)
    t.grp<- cumsum(n.grp)/N

    for(pi in 1:p)	{
      u.grp[,pi]<- tapply(u[,pi], Group, sum)
      W.thetahat[,pi]<- cumsum(u.grp[,pi])/sqrt(N)
    }

    M.thetahat<- solve(J.sqrt)%*%t(W.thetahat)

    #Test statistic & p-value
    Test.stat<- apply(abs(M.thetahat),1,max)
    out.p<- 1- unlist(lapply(Test.stat, FUN=CDF.D))
    out.p.adj<- p.adjust(out.p, method = "hochberg")
    message('Test.statistic=', paste(" ", round(Test.stat, digits=3), sep=''), ',    Adj. p-value=', round(out.p.adj, digits=3), ' \n', appendLF=FALSE)
    sig.stab<- min(out.p.adj)
  }
  ret<- list(pval=c(sig.stab))
  ret
}


plot.LongCART<- function(x, uniform = FALSE, branch = 1, compress = FALSE,
                         nspace=branch, margin = 0, minbranch = 0.3, ...){
  if(!inherits(x, "LongCART")) stop("Need a LongCART object\n")
  class(x)<- "rpart"
  plot(x=x, uniform = uniform, branch = branch, compress = compress,
       nspace=nspace, margin = margin, minbranch = minbranch, ...)
}




text.LongCART<- function(x, splits = TRUE, all = FALSE, use.n = FALSE,
                         minlength = 1L, ...){
  if(!inherits(x, "LongCART")) stop("Need a LongCART object\n")
  class(x)<- "rpart"
  text(x=x, splits = splits , all = all, use.n = use.n, minlength = minlength, ...)
}









