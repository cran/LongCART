SurvCART<- function (data, patid, timevar, censorvar, gvars, tgvars, time.dist = "exponential", 
    cens.dist = "NA", event.ind = 1, alpha = 0.05, minsplit = 40, 
    minbucket = 20, quantile=0.50, print = FALSE) 
{

    node.labels<- function(x){   #the x should be Treeout
        Terminal<- NULL
	tree.nt<- subset(x, !Terminal)
	tree.t<- subset(x, Terminal)
	n.group<- nrow(tree.t)
	id.nt<- tree.nt$ID

	node.labels<- data.frame(node=numeric(), label=character())
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
	node.labels[nrow(node.labels) + 1,]<- c(tree.t$ID[k], text)
	} #--- for loop ends
	node.labels
    }

    nodesummary <- function(data, censorvar, timevar, event.ind, 
        time.dist, quantile=0.5, print = FALSE) {
        Time <- data[[timevar]]
        Delta <- data[[censorvar]]
        D <- sum(I(Delta == event.ind))
        N = nrow(data)
        km.obj1 = survfit(Surv(Time, Delta == event.ind) ~ 1)
        km.obj2 = survfit(Surv(Time, Delta == !event.ind) ~ 1)

	  #--- Summary for Treeout, subject class & prediction
        summary.T = quantile(km.obj1, probs = c(0.25, 0.5, 0.75), conf.int = FALSE)  
        summary.C = quantile(km.obj2, probs = c(0.25, 0.5, 0.75), conf.int = FALSE)  

	  #--- Summary for display
	  yval = quantile(km.obj1, probs = quantile, conf.int = FALSE) 

        if (print) {
            message("\n", appendLF = FALSE)
            message("#of subjects at risk  :", N, "\n", 
                appendLF = FALSE)
            message("#of events            :", D, "\n", 
                appendLF = FALSE)
            message("Event (%)             :", round(D * 
                100/N, 2), "\n", appendLF = FALSE)
            message("Summary (event time)    :", summary.T, 
                "\n", appendLF = FALSE)
            message("Summary (censoring time):", summary.C, 
                "\n", appendLF = FALSE)
            message("\n", appendLF = FALSE)
        }
	  fit.model<- survreg(Surv(Time, I(Delta==event.ind)) ~ 1, dist=time.dist)
	  loglik<- c(logLik(fit.model)) #--- log-likelihood based on fitted parametric model
	  AIC<- c(AIC(fit.model)) #--- log-likelihood based on fitted parametric model

        if(cens.dist!="NA"){
		fit.model.cens<- survreg(Surv(Time, I(Delta!=event.ind)) ~ 1, dist=cens.dist)
		loglik<- loglik + c(logLik(fit.model.cens)) #--- log-likelihood based on fitted parametric model
		AIC<- AIC + c(AIC(fit.model.cens)) #--- log-likelihood based on fitted parametric model
	  }
        ret <- c(N = N, D = D, summary.T, summary.C, yval, loglik, AIC)
        ret
    }
    single.group <- function(data, timevar, censorvar, event.ind, 
        time.dist, splitvar, type, minbucket) {
        spvar <- data[[splitvar]]
        data.t <- data[!is.na(spvar), ]
        Time <- data.t[[timevar]]
        Delta <- data.t[[censorvar]]
        Group <- data.t[[splitvar]]
        fit.model <- loglik.full <- c(logLik(survreg(Surv(Time, 
            I(Delta == event.ind)) ~ 1, dist = time.dist)))
        if (type == 1) {
            event <- I(Delta == event.ind) * 1
            message("Finding cut-off based on heterogeneity in event distribution\n", 
                appendLF = FALSE)
        }
        else if (type == 2) {
            event <- I(Delta != event.ind) * 1
            message("Finding cut-off based on heterogeneity in censoring distribution\n", 
                appendLF = FALSE)
        }
        xcut <- NA
        improve <- NA
        vals <- sort(unlist(unique(Group)))
        vals <- vals[-c(1)]
        message("Evaluations of cutoffs for maximum improvements (Maximum cutoff value = ", 
            max(vals), ")\n", appendLF = FALSE)
        nvals <- length(vals)
        logrank.vec <- NULL
        for (i in 1:nvals) {
            cutoff <- vals[i]
            message(cutoff, appendLF = FALSE)
            group.ind <- ifelse(Group < cutoff, 0, 1)
            if (min(table(group.ind)) >= minbucket) {
                fit <- survdiff(Surv(Time, event) ~ group.ind)
                logrank.vec <- c(logrank.vec, fit$chisq)
            }
            else logrank.vec <- c(logrank.vec, NA)
            message(" ", appendLF = FALSE)
        }
        message(".\n", appendLF = FALSE)
        if (any(!is.na(logrank.vec))) {
            xcut <- vals[which.max(logrank.vec)]
            loglik.left <- c(logLik(survreg(Surv(Time, I(Delta == 
                event.ind)) ~ 1, subset = (Group < xcut), dist = time.dist)))
            loglik.right <- c(logLik(survreg(Surv(Time, I(Delta == 
                event.ind)) ~ 1, subset = (Group >= xcut), dist = time.dist)))
            improve <- loglik.left + loglik.right - loglik.full
        }
        ret <- list(xcut = xcut, improve = improve)
        ret
    }
    bestsplit <- function(data, censorvar, timevar, time.dist, 
        cens.dist, event.ind, gvars, tgvars, node.name, minbucket, 
        alpha) {
        ngvars <- length(gvars)
        best.gvar <- NA
        best.cutoff <- NA
        improve <- NA
        min.pval.adj <- NA
        vi = 0
        stab.res <- matrix(NA, nrow = ngvars, ncol = 2)
        for (v in 1:ngvars) {
            message("\nTesting splitting variable: ", " ", 
                gvars[v], "\n", appendLF = FALSE)
            splitvar <- gvars[v]
            splitvar1 <- data[[gvars[v]]]
            splitvar1 <- splitvar1[!is.na(splitvar1)]
            data1 <- data[!is.na(splitvar1), ]
            G <- length(unique(splitvar1))
            if (G > 1) {
                if (tgvars[v] == 0) 
                  testres <- StabCat.surv(data1, timevar, censorvar, 
                    splitvar, time.dist, cens.dist, event.ind)
                if (tgvars[v] == 1) 
                  testres <- StabCont.surv(data1, timevar, censorvar, 
                    splitvar, time.dist, cens.dist, event.ind)
                stab.res[v, ] <- c(testres$pval, testres$type)
            }
            else message("\nInstability test is NOT performed. \n", 
                appendLF = FALSE)
        }
        stab.pval <- stab.res[, 1]
        if (any(!is.na(stab.pval))) {
            stab.pval.adj <- p.adjust(stab.pval, method = "hochberg")
            sel.v <- which.min(stab.pval)
            min.pval.adj <- stab.pval.adj[sel.v]
            best.gvar <- gvars[sel.v]
            best.gvar.type <- stab.res[sel.v, 2]
            message("\n***Selected Splitting variable:", 
                best.gvar, "***\n", appendLF = FALSE)
            if (min.pval.adj < alpha) {
                best.cut <- single.group(data, timevar, censorvar, 
                  event.ind, time.dist, best.gvar, best.gvar.type, 
                  minbucket)
                best.cutoff <- best.cut$xcut
                improve <- best.cut$improve
                if (is.na(improve)) {
                  best.gvar <- NA
                  best.cutoff <- NA
                  improve <- NA
                  min.pval.adj <- NA
                }
            }
        }
        return(list(node = node.name, gvar = best.gvar, cutoff = best.cutoff, 
            improve = improve, pval = min.pval.adj))
    }
    rsplit <- function(data, censorvar, timevar, time.dist, cens.dist, 
        event.ind, gvars, tgvars, id, split, alpha, minsplit, 
        minbucket, print, nodeout, env = parent.frame()) {
        s.var <- unlist(split[2])
        s.cut <- unlist(split[3])
        s.improve <- unlist(split[4])
        s.pval <- unlist(split[5])
        if (id == 1) {
   
            #--- Updated the env$Treeout
            env$Treeout <- data.frame(id = id, N = nodeout[1], D = nodeout[2], 
                Q1.T=nodeout[3], median.T = nodeout[4], Q3.T = nodeout[5],
                Q1.C=nodeout[6], median.C = nodeout[7], Q3.C = nodeout[8],
                yval=nodeout[9],
                loglik = nodeout[10], AIC = nodeout[11], splitvar = s.var, 
                cutoff = s.cut, pstab = s.pval, improve = s.improve, 
                stringsAsFactors = FALSE)
        }
        else {
            env$Treeout <- rbind(env$Treeout, c(id, nodeout, 
                s.var, s.cut, s.pval, s.improve))
        }
        if (!is.na(s.var) && !is.na(s.cut) && s.pval < alpha) {
            data <- data[!is.na(data[[s.var]]), ]
            id_l = id * 2
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            message("NODE ", id_l, "- Rule:", s.var, 
                " <", s.cut, "\n", appendLF = FALSE)
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            data_l <- subset(data, data[[s.var]] < s.cut)
            left.subj <- nrow(data_l)
            message("No. of subjects in left node: ", left.subj, 
                " \n", appendLF = FALSE)
            nodeout <- nodesummary(data = data_l, censorvar = censorvar, 
                timevar = timevar, event.ind = event.ind, time.dist = time.dist, 
                quantile=quantile, print = print)    #add quantile
            if (left.subj >= minsplit) {
                message("\nDECISION: Go to the next level \n", 
                  appendLF = FALSE)
                split_l <- bestsplit(data_l, censorvar, timevar, 
                  time.dist, cens.dist, event.ind, gvars, tgvars, 
                  id_l, minbucket, alpha)
                rsplit(data_l, censorvar, timevar, time.dist, 
                  cens.dist, event.ind, gvars, tgvars, id_l, 
                  split_l, alpha, minsplit, minbucket, print, 
                  nodeout, env = env)
            }
            else {
                env$Treeout <- rbind(env$Treeout, c(id_l, nodeout, 
                  NA, NA, NA, NA))
            }
            id_r = id * 2 + 1
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            message("NODE ", id_r, "- Rule:", s.var, 
                " >=", s.cut, "\n", appendLF = FALSE)
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            data_r <- subset(data, data[[s.var]] >= s.cut)
            right.subj <- nrow(data_r)
            message("No. of subjects in right node: ", 
                right.subj, " \n", appendLF = FALSE)
            nodeout <- nodesummary(data = data_r, censorvar = censorvar, 
                timevar = timevar, event.ind = event.ind, time.dist = time.dist, 
                quantile=quantile, print = print)    # add quantile
            if (right.subj >= minsplit) {
                message("\nDECISION: Go to the next level \n", 
                  appendLF = FALSE)
                split_r <- bestsplit(data_r, censorvar, timevar, 
                  time.dist, cens.dist, event.ind, gvars, tgvars, 
                  id_r, minbucket, alpha)
                rsplit(data_r, censorvar, timevar, time.dist, 
                  cens.dist, event.ind, gvars, tgvars, id_r, 
                  split_r, alpha, minsplit, minbucket, print, 
                  nodeout, env = env)
            }
            else {
                env$Treeout <- rbind(env$Treeout, c(id_r, nodeout, 
                  NA, NA, NA, NA))
            }
        }
        else message("\nDECISION: NO more splitting required \n", 
            appendLF = FALSE)
    }
    if (!exists(as.character(substitute(data)))) 
        stop("Dataset does not exist\n")
    if (!is.data.frame(data)) 
        stop("Dataset does not exist\n")
    if (!timevar %in% colnames(data)) 
        stop("The column ", timevar, " containing follow-up time is missing in dataset.\n")
    data$timevar <- data[[timevar]]
    if (!censorvar %in% colnames(data)) 
        stop("The column ", censorvar, " containing censoring status is missing in dataset.\n")
    data$censorvar <- data[[censorvar]]
    if (I(time.dist != "exponential") & I(time.dist != 
        "weibull") & I(time.dist != "lognormal") & 
        I(time.dist != "gaussian")) 
        stop("Please check value for time.dist\n")
    if (I(cens.dist != "exponential") & I(cens.dist != 
        "weibull") & I(cens.dist != "NA") & I(cens.dist != 
        "lognormal") & I(cens.dist != "gaussian")) 
        stop("Please check value for cens.dist\n")
    if (!patid %in% colnames(data)) 
        stop("The column ", patid, " containing subject id is missing in dataset.\n")
    data$patid <- data[[patid]]
    data <- data[!is.na(data[[timevar]]) & !is.na(data[[censorvar]]) & 
        !is.na(data[[patid]]), ]
    if (length(gvars) != length(tgvars)) 
        stop("gvars and tgvars are not of equal length. \n")
    if (any(is.na(tgvars))) 
        stop("tgvars cannot have NA value. \n")
    for (var.i in 1:length(gvars)) {
        if (!gvars[var.i] %in% colnames(data)) 
            stop("The column ", gvars[var.i], " is missing in dataset.\n")
    }
    SurvCART.env <- new.env()
    ds<- data
    message("------------------------------------------ \n", 
        appendLF = FALSE)
    message("           ROOT NODE: NODE 1               \n", 
        appendLF = FALSE)
    message("------------------------------------------ \n", 
        appendLF = FALSE)
    message("No. of subjects in root node: ", nrow(data), 
        " \n", appendLF = FALSE)
    rootout <- nodesummary(data, censorvar, timevar, event.ind, 
        time.dist = time.dist, quantile=quantile, print)
    split <- bestsplit(data, censorvar, timevar, time.dist, cens.dist, 
        event.ind, gvars, tgvars, 1, minbucket, alpha)
    rsplit(data, censorvar, timevar, time.dist, cens.dist, event.ind, 
        gvars, tgvars, 1, split, alpha, minsplit, minbucket, 
        print, rootout, env = SurvCART.env)
    Treeout <- SurvCART.env$Treeout
    colnames(Treeout) <- c("ID", "n", "D", 
        "Q1.T", "median.T", "Q3.T", 
        "Q1.C", "median.C", "Q3.C",
        "yval", 
        "loglik", "AIC", 
        "var", "index", "p (Instability)", 
        "improve")
    row.names(Treeout) <- NULL
    Treeout[, 1] <- as.numeric(Treeout[, 1])
    Treeout[, 2] <- as.numeric(Treeout[, 2])
    Treeout[, 3] <- as.numeric(Treeout[, 3])
    Treeout[, 4] <- round(as.numeric(Treeout[, 4]), digits = 3)
    Treeout[, 5] <- round(as.numeric(Treeout[, 5]), digits = 3)
    Treeout[, 6] <- round(as.numeric(Treeout[, 6]), digits = 3)
    Treeout[, 7] <- round(as.numeric(Treeout[, 7]), digits = 3)
    Treeout[, 8] <- round(as.numeric(Treeout[, 8]), digits = 3)
    Treeout[, 9] <- round(as.numeric(Treeout[, 9]), digits = 3)
    Treeout[, 10] <- round(as.numeric(Treeout[, 10]), digits = 3)
    Treeout[, 11] <- round(as.numeric(Treeout[, 11]), digits = 1)
    Treeout[, 12] <- round(as.numeric(Treeout[, 12]), digits = 1)
    Treeout[, 14] <- round(as.numeric(Treeout[, 14]))
    Treeout[, 15] <- round(as.numeric(Treeout[, 15]), digits = 3)
    Treeout[, 16] <- round(as.numeric(Treeout[, 16]), digits = 1)
    Treeout$Terminal <- ifelse(is.na(Treeout[, 14]), TRUE, FALSE)
    print(Treeout)

    logLik.tree <- sum(Treeout$loglik * Treeout$Terminal)
    logLik.root <- Treeout[1, "loglik"]
    improve.loglik <- logLik.tree - logLik.root
    AIC.tree <- sum(Treeout$AIC * Treeout$Terminal)
    AIC.root <- Treeout[1, "AIC"]
    message("logLikelihood (root)=", logLik.root, "   logLikelihood (tree)=", 
        logLik.tree, "\n", appendLF = FALSE)
    message("AIC (root)=", AIC.root, "    AIC (tree)=", 
        AIC.tree, "\n", appendLF = FALSE)
    Treeout1 <- as.data.frame(Treeout[, c("ID", "var", 
        "n", "D", "yval",        #--- updated 
        "loglik", "AIC", "Terminal", "index", 
        "improve")])
    row.names(Treeout1) <- Treeout1$ID
    Treeout1$var = ifelse(Treeout1$Terminal, "<leaf>", 
        Treeout1$var)
    Treeout1$yval = paste("q", quantile, "=", Treeout1$yval, sep = "")  #--- updated
    Treeout1$dev = -Treeout1$loglik
    Treeout1$wt <- Treeout1$count <- Treeout1$n
    Treeout1$ncompete <- Treeout1$nsurrogate <- 0
    Treeout1$ncat = -1
    frame <- Treeout1[c("var", "wt", "dev", 
        "yval", "ncompete", "nsurrogate")]
    frame$n <- paste(Treeout1$n, ",D=", Treeout1$D, sep = "")
    splits <- Treeout1[!Treeout1$Terminal, ]
    splits <- splits[c("count", "ncat", "improve", 
        "index")]
    splits <- as.matrix(splits)
    temp <- frame[frame$var == "<leaf>", ]
    cptable <- 0:(dim(temp)[1] - 1)
    functions <- NULL
    functions$text <- function(yval, dev, wt, ylevel, digits, 
        n, use.n) {
        if (use.n) 
            paste(yval, "\nN=", n, "")
        else yval
    }

    nodelab<- node.labels(Treeout)

    varnam<- unique(frame$var)
    varnam<- varnam[varnam!="<leaf>" ]

    ret <- list(Treeout, frame, splits, 
                cptable, functions, logLik.tree, 
                logLik.root, AIC.tree, AIC.root, 
                event.ind, nodelab, varnam,
                ds, timevar, censorvar)

    names(ret) <- c("Treeout", "frame", "splits", 
        "cptable", "functions", "logLik.tree", 
        "logLik.root", "AIC.tree", "AIC.root", 
        "event.ind", "nodelab", "varnam",
        "ds", "timevar", "censorvar")   
    class(ret) <- c("SurvCART")
    ret
}
