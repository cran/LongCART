ProfilePlot<- function (x, timevar, timevar.power = NULL, covariate.val = NULL, 
    xlab = NULL, ylab = NULL, sg.title = 4, mfrow=NULL, ...) 
{
    if (!inherits(x, "LongCART")) 
        stop("Need a LongCART object\n")
    time.range <- range(x$data[[timevar]])
    profile.data <- data.frame(var1 = seq(time.range[1], time.range[2], 
        length = 100))
    names(profile.data) <- timevar
    form <- Formula(x$fixed)
    Y.name <- as.character(attr(form, "lhs"))
    X.name1 <- as.character(attr(form, "rhs"))
    X.name1 <- gsub(pattern = " ", replacement = "", 
        x = X.name1)
    X.name1 <- gsub(pattern = "-1", replacement = "", 
        x = X.name1)
    X.name <- unlist(strsplit(X.name1, "+", fixed = TRUE))
    p <- length(X.name)
    for (p.i in 1:p) {
        var.i <- X.name[p.i]
        if (is.null(timevar.power)) 
            indic = 0
        else if (is.na(timevar.power[p.i])) 
            indic = 1
        else indic = 2
        if (var.i != timevar) {
            if (indic == 2) 
                profile.data[[var.i]] <- (profile.data[[timevar]])^timevar.power[p.i]
            else if (indic < 2) {
                if (is.null(covariate.val)) 
                  profile.data[[var.i]] <- median(x$data[[var.i]], 
                    na.rm = TRUE)
                else if (is.na(covariate.val[p.i])) 
                  profile.data[[var.i]] <- median(x$data[[var.i]], 
                    na.rm = TRUE)
                else profile.data[[var.i]] <- covariate.val[p.i]
            }
        }
    }
    Tout <- x$Treeout
    Tout <- Tout[Tout$Terminal, ]
    n.sg <- nrow(Tout)
    sg.i <- 1
    crit <- Tout$yval[sg.i]
    fit.val <- with(profile.data, eval(parse(text = crit)))
    if (is.null(xlab)) 
        xlab = timevar
    if (is.null(ylab)) 
        ylab = Y.name

    nodes<- x$nodelab

    if (sg.title == 1) 
        Tout$main = paste0("Node=", Tout$ID)
    if (sg.title == 2) 
        Tout$main = paste0("Sub-group=", 1:n.sg)
    if (sg.title == 3) 
        Tout$main = paste0("Sub-group=", 1:n.sg, " (Node=", 
            Tout$ID, ")")
    if (sg.title == 4) 
        Tout$main = paste0(nodes[,1], ": ", nodes[,2])

    y.range <- NULL
    for (sg.i in 1:n.sg) {
        crit <- Tout$yval[sg.i]
        fit.val <- with(profile.data, eval(parse(text = crit)))
        y.range <- range(y.range, fit.val)
    }

    if(is.null(mfrow)) mfrow=c( ceiling(sqrt(n.sg)), ceiling(sqrt(n.sg))) 
    par(mfrow=mfrow)
    for (sg.i in 1:n.sg) {
        crit <- Tout$yval[sg.i]
        fit.val <- with(profile.data, eval(parse(text = crit)))
        plot(profile.data$time, fit.val, type = "l", xlab = xlab, 
            ylab = ylab, main = Tout$main[sg.i], ylim = y.range, 
            ...)
    }
}
