LongCART<- function (data, patid, fixed, gvars, tgvars, minsplit = 40, minbucket = 20, 
    alpha = 0.05, coef.digits = 2, print.lme = FALSE) 
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

    llike <- function(data, fixed) {
        if (exists("txx.out")) 
            rm(txx.out)
        loglik <- NA
        try(txx.out <- lme(fixed, data = data, method = "ML", 
            random = ~1 | patid, na.action = na.omit), silent = TRUE)
        if (exists("txx.out")) 
            loglik <- logLik(txx.out)[1]
        loglik
    }
    test4split <- function(cutoff, data, fixed, splitvar) {
        LRTS <- NA
        message(cutoff, appendLF = FALSE)
        left.data <- data[data[[splitvar]] < cutoff, ]
        right.data <- data[data[[splitvar]] >= cutoff, ]
        rn.ll <- llike(data = data, fixed = fixed)
        left.ll <- llike(data = left.data, fixed = fixed)
        right.ll <- llike(data = right.data, fixed = fixed)
        LRTS <- 2 * (left.ll + right.ll - rn.ll)
        message(" ", appendLF = FALSE)
        LRTS
    }
    single.group <- function(data, fixed, splitvar, minbucket) {
        xcut <- NA
        improve <- NA
        Y.name <- as.character(attr(as.Formula(fixed), "lhs"))
        data <- data[!is.na(data[[Y.name]]), ]
        data <- data[!is.na(data[[splitvar]]), ]
        temp <- unique(data[, c("patid", splitvar)])
        N <- tapply(temp$patid, temp[[splitvar]], length)
        temp1 <- as.data.frame(list(cuts = as.numeric(names(N)), 
            N = N))
        temp1$left <- cumsum(temp1$N) - temp1$N
        temp1$right <- sum(temp1$N) - temp1$left
        temp2 <- temp1[temp1$left >= minbucket & temp1$right >= 
            minbucket, ]
        vals <- temp2$cuts
        message("Evaluations of cutoffs for maximum improvements (Maximum cutoff value = ", 
            max(vals), ")\n", appendLF = FALSE)
        LRT.vec <- lapply(vals, FUN = test4split, data = data, 
            fixed = fixed, splitvar = splitvar)
        message(".\n", appendLF = FALSE)
        LRT.vec <- unlist(LRT.vec)
        if (any(!is.na(LRT.vec))) {
            index <- which.max(LRT.vec)
            xcut = vals[index]
            improve = LRT.vec[index]
        }
        ret <- list(xcut = xcut, improve = improve)
        ret
    }
    bestsplit <- function(data, fixed, gvars, tgvars, node.name, 
        minbucket, alpha) {
        ngvars <- length(gvars)
        best.gvar <- NA
        best.cutoff <- NA
        improve <- NA
        min.pval.adj <- NA
        stab.pval <- numeric(length = ngvars)
        for (v in 1:ngvars) {
            stab.pval[v] <- NA
            splitvar <- gvars[v]
            message("\nSplitting variable: ", splitvar, 
                "\n", appendLF = FALSE)
            data1 <- data[!is.na(splitvar), ]
            G <- length(unique(data1[[splitvar]]))
            message("G=", G, "\n", appendLF = FALSE)
            if (G > 1) {
                if (tgvars[v] == 0) 
                  stab.pval[v] <- StabCat(data = data1, patid = "patid", 
                    fixed = fixed, splitvar = splitvar)$pval
                if (tgvars[v] == 1) 
                  stab.pval[v] <- StabCont(data = data1, patid = "patid", 
                    fixed = fixed, splitvar = splitvar)$pval
            }
            else message("Instability test was NOT performed. \n", 
                appendLF = FALSE)
        }
        message("\n stab.pval=", stab.pval, "\n", 
            appendLF = FALSE)
        if (any(!is.na(stab.pval))) {
            stab.pval.adj <- p.adjust(stab.pval, method = "hochberg")
            sel.v <- which.min(stab.pval)
            min.pval.adj <- stab.pval.adj[sel.v]
            best.gvar <- gvars[sel.v]
            message("\n stab.pval.adj=", stab.pval.adj, 
                "\n", appendLF = FALSE)
            message("\n alpha=", alpha, "\n", appendLF = FALSE)
            if (min.pval.adj < alpha) {
                bestcut <- single.group(data, fixed, best.gvar, 
                  minbucket)
                best.cutoff <- bestcut$xcut
                improve <- bestcut$improve
            }
        }
        return(list(node = node.name, gvar = best.gvar, cutoff = best.cutoff, 
            improve = improve, pval = min.pval.adj))
    }
    coeff.Estimate <- function(LMEobj, coef.digits = 2) {
        ct <- LMEobj$tTable
        p <- nrow(ct)
        ct.names <- row.names(ct)
        yval <- " "
        for (ct.i in 1:p) yval <- paste0(yval, "+", round(ct[ct.i, 
            1], coef.digits), "*", ct.names[ct.i])
        yval <- gsub(pattern = "(Intercept)", replacement = "", 
            x = yval, fixed = TRUE)
        yval <- gsub(pattern = "+-", replacement = "-", 
            yval, fixed = TRUE)
        yval <- gsub(pattern = " +", replacement = "", 
            yval, fixed = TRUE)
        yval <- gsub(pattern = " ", replacement = "", 
            yval, fixed = TRUE)
        yval <- gsub(pattern = "*+", replacement = "+", 
            yval, fixed = TRUE)
        yval <- gsub(pattern = "*-", replacement = "-", 
            yval, fixed = TRUE)
        yval
    }
    rsplit <- function(data, fixed, gvars, tgvars, id, split, 
        alpha, minsplit, minbucket, Rate, loglik, env = parent.frame(), 
        coef.digits = 2) {
        s.var <- unlist(split[2])
        s.cut <- unlist(split[3])
        s.improve <- unlist(split[4])
        s.pval <- unlist(split[5])
        N <- length(unique(data$patid))
        if (id == 1) {
            env$Treeout <- data.frame(id = id, N = N, yval = Rate, 
                splitvar = s.var, cutoff = s.cut, pstab = s.pval, 
                loglik = loglik, improve = s.improve, stringsAsFactors = FALSE)
        }
        else {
            env$Treeout <- rbind(env$Treeout, c(id, N, Rate, 
                s.var, s.cut, s.pval, loglik, s.improve))
        }
        if (!is.na(s.var) && !is.na(s.cut)) {
            data <- data[!is.na(data[[s.var]]), ]
            id_l = id * 2
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            message("NODE ", id_l, "- Rule:", s.var, 
                " <", s.cut, "\n", appendLF = FALSE)
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            data_l <- subset(data, data[[s.var]] < s.cut)
            left.subj <- length(unique(data_l$patid))
            message("No. of individual in left node: ", 
                left.subj, " \n", appendLF = FALSE)
            tjj <- NULL
            tjj <- summary(lme(fixed, data = data_l, method = "REML", 
                random = ~1 | patid))
            if (print.lme) 
                cat("NODE ", id_l, "- Rule:", s.var, 
                  " <", s.cut, "\n")
            if (print.lme) 
                print(tjj$tTable)
            Rate <- coeff.Estimate(tjj, coef.digits = coef.digits)
            loglik <- logLik(lme(fixed, data = data_l, method = "REML", 
                random = ~1 | patid))
            if (left.subj >= minsplit) {
                message("\nDECISION: Go to the next level \n", 
                  appendLF = FALSE)
                split_l <- bestsplit(data = data_l, fixed = fixed, 
                  gvars = gvars, tgvars = tgvars, node.name = id_l, 
                  minbucket = minbucket, alpha = alpha)
                rsplit(data = data_l, fixed = fixed, gvars = gvars, 
                  tgvars = tgvars, id = id_l, split = split_l, 
                  alpha = alpha, minsplit = minsplit, minbucket = minbucket, 
                  Rate = Rate, loglik = loglik, env = env, coef.digits = coef.digits)
            }
            else {
                env$Treeout <- rbind(env$Treeout, c(id_l, left.subj, 
                  Rate, NA, NA, NA, loglik, NA))
            }
            id_r = id * 2 + 1
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            message("NODE ", id_r, "- Rule:", s.var, 
                " >=", s.cut, "\n", appendLF = FALSE)
            message("---------------------------------------- \n", 
                appendLF = FALSE)
            data_r <- subset(data, data[[s.var]] >= s.cut)
            right.subj <- length(unique(data_r$patid))
            message("No. of individual in right node: ", 
                right.subj, " \n", appendLF = FALSE)
            tjj <- NULL
            tjj <- summary(lme(fixed, data = data_r, method = "REML", 
                random = ~1 | patid))
            if (print.lme) 
                cat("NODE ", id_r, "- Rule:", s.var, 
                  " >=", s.cut, "\n")
            if (print.lme) 
                print(tjj$tTable)
            Rate <- coeff.Estimate(tjj, coef.digits = coef.digits)
            loglik <- logLik(lme(fixed, data = data_r, method = "REML", 
                random = ~1 | patid))
            if (right.subj >= minsplit) {
                message("\nDECISION: Go to the next level \n", 
                  appendLF = FALSE)
                split_r <- bestsplit(data = data_r, fixed = fixed, 
                  gvars = gvars, tgvars = tgvars, node.name = id_r, 
                  minbucket = minbucket, alpha = alpha)
                rsplit(data = data_r, fixed = fixed, gvars = gvars, 
                  tgvars = tgvars, id = id_r, split = split_r, 
                  alpha = alpha, minsplit = minsplit, minbucket = minbucket, 
                  Rate = Rate, loglik = loglik, env = env, coef.digits = coef.digits)
            }
            else {
                env$Treeout <- rbind(env$Treeout, c(id_r, right.subj, 
                  Rate, NA, NA, NA, loglik, NA))
            }
        }
        else message("\nDECISION: NO more splitting required \n", 
            appendLF = FALSE)
    }
    if (!exists(as.character(substitute(data)))) 
        stop("Dataset does not exist\n")
    if (!is.data.frame(data)) 
        stop("Dataset does not exist\n")
    if (!patid %in% colnames(data)) 
        stop("The column ", patid, " containing subjects id is missing in dataset.\n")
    data$patid <- data[[patid]]
    data <- data[!is.na(data[["patid"]]), ]
    Y.name <- as.character(attr(as.Formula(fixed), "lhs"))
    if (!Y.name %in% colnames(data)) 
        stop("The column ", Y.name, " containing subjects id is missing in dataset.\n")
    data <- data[!is.na(data[[Y.name]]), ]
    if (length(gvars) != length(tgvars)) 
        stop("gvars and tgvars are not of equal length. \n")
    if (any(is.na(tgvars))) 
        stop("tgvars cannot have NA value. \n")
    for (var.i in 1:length(gvars)) {
        if (!gvars[var.i] %in% colnames(data)) 
            stop("The column ", gvars[var.i], " is missing in dataset.\n")
    }
    for (var.i in 1:length(gvars)) {
        temp <- as.data.frame(list(patid = cbind(data["patid"], 
            gvar = data[gvars[var.i]])))
        colnames(temp) = c("patid", "gvar")
        temp <- unique(temp)
        temp <- temp[!is.na("gvar"), ]
        n.gvar.pat <- tapply(temp$gvar, temp$patid, function(x) length(unique(x)))
        if (any(n.gvar.pat != 1)) 
            stop("One subject is associated with more than one value of ", 
                gvars[var.i], ".\n")
        rm(temp, n.gvar.pat)
    }
    LongCART.env <- new.env()
    message("------------------------------------------ \n", 
        appendLF = FALSE)
    message("           ROOT NODE: NODE 1               \n", 
        appendLF = FALSE)
    message("------------------------------------------ \n", 
        appendLF = FALSE)
    tjj <- NULL
    tjj <- summary(lme(fixed, data = data, method = "REML", 
        random = ~1 | patid))
    if (print.lme) 
        cat("ROOT NODE: NODE 1\n")
    if (print.lme) 
        print(tjj$tTable)
    Rate <- coeff.Estimate(tjj, coef.digits = coef.digits)
    loglik <- logLik(lme(fixed, data = data, method = "REML", 
        random = ~1 | patid))
    split <- bestsplit(data = data, fixed = fixed, gvars = gvars, 
        tgvars = tgvars, node.name = 1, minbucket = minbucket, 
        alpha = alpha)
    rsplit(data = data, fixed = fixed, gvars = gvars, tgvars = tgvars, 
        id = 1, split = split, alpha = alpha, minsplit = minsplit, 
        minbucket = minbucket, Rate = Rate, loglik = loglik, 
        env = LongCART.env, coef.digits = coef.digits)
    Treeout <- LongCART.env$Treeout
    colnames(Treeout) <- c("ID", "n", "yval", 
        "var", "index", "p (Instability)", 
        "loglik", "improve")
    row.names(Treeout) <- NULL
    Treeout[, 1] <- as.numeric(Treeout[, 1])
    Treeout[, 2] <- as.numeric(Treeout[, 2])
    Treeout[, 5] <- round(as.numeric(Treeout[, 5]), digits = 2)
    Treeout[, 6] <- round(as.numeric(Treeout[, 6]), digits = 3)
    Treeout[, 7] <- round(as.numeric(Treeout[, 7]), digits = 0)
    Treeout[, 8] <- round(as.numeric(Treeout[, 8]), digits = 0)
    Treeout$Terminal <- ifelse(is.na(Treeout[, 5]), TRUE, FALSE)
    print(Treeout)
    form <- Formula(fixed)
    Y.name <- as.character(attr(form, "lhs"))
    X.name <- as.character(attr(form, "rhs"))
    X.name <- gsub(pattern = " ", replacement = "", 
        x = X.name)
    X.name <- unlist(strsplit(X.name, "+", fixed = TRUE))
    p <- length(X.name) + 1
    AIC.tree <- 2 * sum(Treeout[, 7] * Treeout[, 9]) - 2 * (p + 
        2) * sum(Treeout[, 9])
    AIC.root <- 2 * Treeout[1, 7] - 2 * (p + 2)
    improve.AIC <- AIC.tree - AIC.root
    logLik.tree <- sum(Treeout[, 7] * Treeout[, 9])
    logLik.root <- Treeout[1, 7]
    Deviance <- 2 * (logLik.tree - logLik.root)
    LRT.df <- p * sum(Treeout[, 9]) - p
    LRT.p <- 1 - pchisq(Deviance, LRT.df)
    message("AIC(tree)=", AIC.tree, "   AIC(Root)=", 
        AIC.root, "\n", appendLF = FALSE)
    message("logLikelihood (tree)=", logLik.tree, "   logLikelihood (Root)=", 
        logLik.root, "\n", appendLF = FALSE)
    message("Deviance=", Deviance, " (df=", LRT.df, 
        ", p-val=", LRT.p, ") \n", appendLF = FALSE)
    sel <- as.numeric(rownames(Treeout[Treeout$Terminal, ]))
    Treeout1 <- as.data.frame(Treeout[, c("ID", "var", 
        "n", "yval", "loglik", "Terminal", 
        "index", "improve")])
    row.names(Treeout1) <- Treeout1$ID
    Treeout1$var = ifelse(Treeout1$Terminal, "<leaf>", 
        Treeout1$var)
    Treeout1$dev = -Treeout1$loglik
    Treeout1$wt <- Treeout1$count <- Treeout1$n
    Treeout1$ncompete <- Treeout1$nsurrogate <- 0
    Treeout1$ncat = -1
    frame <- Treeout1[c("var", "n", "wt", "dev", 
        "yval", "ncompete", "nsurrogate")]
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
            paste0(yval, "\nn=", n)
        else yval
    }

    nodelab<- node.labels(Treeout)

    varnam<- unique(frame$var)
    varnam<- varnam[varnam!="<leaf>" ]

    ret <- list(Treeout, frame, splits, cptable, functions, p, 
        AIC.tree, AIC.root, improve.AIC, logLik.tree, logLik.root, 
        Deviance, LRT.df, LRT.p, data, patid, fixed,
        nodelab, varnam)
    names(ret) <- c("Treeout", "frame", "splits", 
        "cptable", "functions", "p", "AIC.tree", 
        "AIC.root", "improve.AIC", "logLik.tree", 
        "logLik.root", "Deviance", "LRT.df", 
        "LRT.pval", "data", "patid", "fixed",
        "nodelab", "varnam")
    class(ret) <- "LongCART"
    ret
}
