

#    Add node.labels() function within SurvCART() and LongCART()
#    create a object nodelab=node.labels(Treeout)
#    Add this code 	
#      varnam<- unique(frame$var)
#	 varnam<- varnam[varnam!="<leaf>" ]
#    Update "ret" object with varname and nodelab in LongCART() and SurvCART()
#    Replace "nodes<- node.labels(x$Treeout)" by "nodes<- x$nodelab"
#    Update help files for LongCART() and SurvCART() and KMSurvPlot(): Input, output and example
#    Update the ProfilePlot.LongCART if necessary
  

KMPlot<- function (x, type = 1, overlay=TRUE, conf.type="log-log", mfrow=NULL, ...) 
{
    SubGroup<- NULL
    if (!inherits(x, "SurvCART")) 
        stop("Need a SurvCART object\n")
    ds <- predict.SurvCART(object=x, newdata=x$ds)
    ds2<- data.frame(time=ds[[x$timevar]], cens=ds[[x$censorvar]], SubGroup=ds$node)
    ds2$event <- I(ds2$cens == x$event.ind) * 1
    ds2$censor <- 1 - ds2$event
    if(overlay){
    	if (type == 1) km.out <- survfit(Surv(time, event) ~ SubGroup, data = ds2)
    	if (type == 2) km.out <- survfit(Surv(time, censor) ~ SubGroup, data = ds2)
    	ggsurvplot(km.out, data = ds2, ...)
    } else {
    	nodes<- x$nodelab
    	n.nodes<- nrow(nodes)
    	if(is.null(mfrow)) mfrow=c( ceiling(sqrt(n.nodes)), ceiling(sqrt(n.nodes))) 
    	par(mfrow=mfrow)
    	for(k in 1:n.nodes){
		if (type == 1) km.obj <- survfit(Surv(time, event) ~ 1, data = ds2, 
                              subset=(SubGroup==nodes$node[k]), conf.type = conf.type)
		if (type == 2) km.obj <- survfit(Surv(time, censor) ~ 1, data = ds2, 
                              subset=(SubGroup==nodes$node[k]), conf.type = conf.type)
		#title(main=paste0(nodes[k,1], ": ", nodes[k,2]))conf
            plot(km.obj, 
                 main=paste0(nodes[k,1], ": ", nodes[k,2]), ...)
                 #xlab="Time", ylab="Survival Probability", 
		rm(km.obj)
    	}
    }
}


#KMPlot.SurvCART(out, xscale=365.25, type=1)
#KMPlot.SurvCART(out, xscale=365.25, type=2, overlay=FALSE, mfrow=c(2,2), xlab="abc", ylab="xyz")

