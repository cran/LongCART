
predict.LongCART<- function(object, newdata, patid, ...){ 
                            #--- fit: a SurvCART obj
                            #--- data: dataset to be predicted
        if (!inherits(object, "LongCART")) 
        stop("Not a legitimate \"LongCART\" object")

       Terminal<- NULL

        fit<- object
        data<- newdata
	frame<- fit$frame
	varnam<- fit$varnam
	vnum <- match(varnam, colnames(data))   
	if (any(is.na(vnum))) 
        stop("Tree has variables not found in new data")

	treeobj<- fit$Treeout
	tree.nt<- subset(treeobj, !Terminal)
	tree.t<- subset(treeobj, Terminal)
	n.group<- nrow(tree.t)
	id.nt<- tree.nt$ID

	nd<- data
	nd$temp.lc.patid<- nd[,patid]   #add a temporary id variable to associate with sub-group info
	N<- unique(nd$temp.lc.patid)
	pred.group<- data.frame(temp.lc.patid=numeric(),
                 			node=integer(),
                 			profile=character(),
                 			stringsAsFactors=FALSE)

	for(k in 1:n.group){
		node.t<- tree.t$ID[k]
		#cat(node.t, "\n")
		text=NULL
		while(node.t>1){
			if((node.t/2)>floor(node.t/2)){
				text<- paste0(tree.nt$var[id.nt==floor(node.t/2)], ">=", tree.nt$index[id.nt==floor(node.t/2)], 
                          ifelse(!is.null(text), paste0(" & ", text), ""))
			} else if((node.t/2)==floor(node.t/2)){
				text<- paste0(tree.nt$var[id.nt==node.t/2], "<", tree.nt$index[id.nt==node.t/2],  
                          ifelse(!is.null(text), paste0(" & ", text), ""))
			} 
		#cat(text, "\n")
			node.t=floor(node.t/2)
		} #--- while loop ends
		selsubj<- subset(nd, eval(parse(text=text)))$temp.lc.patid
      	if(length(selsubj)>0)
			pred.group<- rbind(pred.group, 
                               	 data.frame(temp.lc.patid=selsubj,
                 						node=tree.t$ID[k],
                 						profile=tree.t$yval[k],
                 						stringsAsFactors=FALSE))
	} #--- for loop ends

pred.obj<- merge(nd, pred.group, by="temp.lc.patid")
pred.obj<- transform(pred.obj, predval=eval(parse(text=profile)))
pred.obj[,-c(1)]
}





