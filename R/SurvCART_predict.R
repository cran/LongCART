
#--- Replace the following line in SurvCART()
#         Treeout1$yval = paste("Med=", Treeout1$median.T, sep = "")
#    by
#        Treeout1$yval = paste("q", quantile, "=", Treeout1$yval, sep = "")  #--- replace median.T by yval at the very outset
#  
#---  Add Q1.T, Q1.C, Q3.T, Q3.C to the subj.class and Treeout obj
#
#---  Add quantile=0.5 as a input parameter
#
#--- Revise the node.summary function
#
#--- Add quantile=quantile when calling node.summary function
#
#--- Delete env$idlist
#
#--- Specify in the help file that AIC is actually -2*loglikelihood + 2


predict.SurvCART<- function(object, newdata, ...){ 
                            #--- fit: a SurvCART obj
                            #--- data: dataset to be predicted
        if (!inherits(object, "SurvCART")) 
        stop("Not a legitimate \"SurvCART\" object")

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
	N<- nrow(nd)
	nd$temp.sc.patid<- 1:N   #add a temporary id variable to associate with sub-group info
	pred.group<- NULL

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
		selsubj<- subset(nd, eval(parse(text=text)))$temp.sc.patid
      	if(length(selsubj)>0)
			pred.group<- rbind(pred.group, 
                               cbind(temp.sc.patid=selsubj, 
                                     node=tree.t$ID[k], 
                                     median.T=tree.t$median.T[k],
                                     median.C=tree.t$median.C[k],
                                     Q1.T=tree.t$Q1.T[k],  Q3.T=tree.t$Q3.T[k], 
                                     Q1.C=tree.t$Q1.C[k],  Q3.C=tree.t$Q3.C[k]))
	} #--- for loop ends

pred.obj<- merge(nd, pred.group, by="temp.sc.patid")
pred.obj[,-c(1)]
}

