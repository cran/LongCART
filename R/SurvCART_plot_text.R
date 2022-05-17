
plot.SurvCART<- function(x, uniform = FALSE, branch = 1, compress = FALSE, 
                         nspace=branch, margin = 0, minbranch = 0.3, ...){
      if(!inherits(x, "SurvCART")) stop("Need a SurvCART object\n")
      class(x)<- "rpart"
	plot(x=x, uniform = uniform, branch = branch, compress = compress,
           nspace=nspace, margin = margin, minbranch = minbranch, ...)
}




text.SurvCART<- function(x, splits = TRUE, all = FALSE, use.n = FALSE, 
                minlength = 1L, ...){
	if(!inherits(x, "SurvCART")) stop("Need a SurvCART object\n")
      class(x)<- "rpart"
	text(x=x, splits = splits , all = all, use.n = use.n, minlength = minlength, ...)
}












