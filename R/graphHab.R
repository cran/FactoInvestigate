graphHab <-
function(res, file = "", dim = 1:2, hab = NULL, ellipse = TRUE, Iselec = "contrib", Rselec = "cos2", Cselec = "contrib", Icoef = 1, Rcoef = 1, Ccoef = 1, figure.title = "Figure", graph = TRUE, cex = 0.7, options=NULL) {

test.de.Wilks <- function(x,grouping){
  if (any(summary(grouping)<2)){
    notok=grouping%in%(levels(grouping)[which(summary(grouping)<2)])
    aux <- Wilks.test(x[!notok,,drop=FALSE],grouping[!notok])$parameter
  } else aux <- Wilks.test(x,grouping)$parameter
  pchisq(aux[1],aux[2],lower.tail=FALSE)
}

    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(Iselec) & !is.character(Iselec)) {return(warning("the argument 'Iselec' should be a numeric or character vector"))}
    if(!is.numeric(Rselec) & !is.character(Rselec)) {return(warning("the argument 'Rselec' should be a numeric or character vector"))}
    if(!is.numeric(Cselec) & !is.character(Cselec)) {return(warning("the argument 'Cselec' should be a numeric or character vector"))}
    
    if(!is.numeric(Icoef)) {return(warning("the argument 'Icoef' must be numeric"))}
    if(!is.numeric(Rcoef)) {return(warning("the argument 'Rcoef' must be numeric"))}
    if(!is.numeric(Ccoef)) {return(warning("the argument 'Ccoef' must be numeric"))}
    
    if(Icoef < 0) {return(warning("the argument 'Icoef' must be positive"))}
    if(Rcoef < 0) {return(warning("the argument 'Rcoef' must be positive"))}
    if(Ccoef < 0) {return(warning("the argument 'Ccoef' must be positive"))}
    
    if(!is.numeric(hab) & !is.character(hab) & !is.null(hab)) {return(warning("the argument 'hab' should be the name or the index of the variable used to color the individuals"))}
    
    if(!is.numeric(cex)) {return(warning("the argument 'cex' must be numeric"))}
    if(cex < 0) {return(warning("the argument 'cex' must be positive"))}
    
    if(!is.logical(ellipse)) {return(warning("the argument 'ellipse' must be logical"))}
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) | length(dim) != 2) {return(warning("the argument 'dim' has to be a 2 dimensional numeric vector"))}
    if(any(dim < 0)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    data = param$data
    hab.param = hab
    switch(analyse,
           PCA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Iselec, coef = Icoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             
             factors = sapply(rownames(res$quali.sup$eta2), function(x) {length(levels(data[, x]))})
             
             quali = names(factors)[factors > 1]
             reject = names(factors)[factors <= 1]
             
             wilks.p = sapply(quali, function(x, res, dim) {test.de.Wilks(rbind(res$ind$coord[, dim[1]:dim[2]], res$ind.sup$coord[, dim[1]:dim[2]]), res$call$X[, x])}, res = res, dim = dim)
             names(wilks.p) = quali
             wilks.p = sort(wilks.p)
             
             hab = names(which(wilks.p == min(wilks.p))) # maybe more than 1
             p.value = min(wilks.p)
             
             writeRmd(gettext("The Wilks test p-value indicates which variable factors are the best separated on the plane",domain="R-FactoInvestigate"), " (",
                      gettext("i.e. which one explain the best the distance between individuals",domain="R-FactoInvestigate"), ")", sep = "", file = file, end = ".\n")
             
              # wilks.s = NULL
             # if(length(hab) > 1) {
               # wilks.s = sapply(names(wilks.p[wilks.p == 0]), function(x, res, dim) {test.de.Wilks(rbind(res$ind$coord[, dim[1]:dim[2]], res$ind.sup$coord[, dim[1]:dim[2]]), res$call$X[, x])$statistic}, res = res, dim = dim)
               # names(wilks.s) = quali[wilks.p == 0]
               # wilks.s = sort(wilks.s, decreasing = TRUE)
               
               # if(length(wilks.s) > 12) {wilks.s = wilks.s[1:12]}
               
               # hab = names(which.max(wilks.s))
             # }
             
             if(length(wilks.p) > 12) {wilks.p = wilks.p[1:12]}
             
             
             if(graph) {
               show(wilks.p)
             }
             writeRmd(start = TRUE, end = "", options = options, file = file)
             dump("wilks.p", file = file, append = TRUE)
             writeRmd("wilks.p", stop = TRUE, sep = "", file = file)
             
             
             if(length(quali) == 1) {
               writeRmd(gettext("There only is one possible qualitative variable to illustrate the distance between individuals",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
             } else {
               if(length(names(which(wilks.p == min(wilks.p)))) == 1) {
                 writeRmd(gettext("The best qualitative variable to illustrate the distance between individuals on this plane is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               } else {
                 writeRmd(gettext("Many qualitative variables have a Wilks p-value equal to zero",domain="R-FactoInvestigate"), ". ", gettext("To arbitrate which one to select, we need to compare their statistic value",domain="R-FactoInvestigate"), 
                          " : *", paste(names(which(wilks.p == min(wilks.p))), collapse = "*, *"), file = file, sep = "", end = "*.\n")
                 
                 # if(graph) {
                   # show(wilks.s)
                 # }
                 # writeRmd(start = TRUE, end = "", options = options, file = file)
                 # dump("wilks.s", file = file, append = TRUE)
                 # writeRmd("wilks.s", stop = TRUE, sep = "", file = file)
                 
                 # writeRmd(gettext("The best qualitative variable to illustrate the distance between individuals on this plane is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               }
             }
             
             if(length(reject) != 0) {
               writeRmd(gettext("The qualitative variables",domain="R-FactoInvestigate"), " *", paste(reject, collapse = "*, *"), "* ",
                        gettext("cannot separate the individuals on the plane, cause they are unimodal",domain="R-FactoInvestigate"), end = ".\n", file = file, sep = "")
             }
             
             if(!is.null(hab.param)) {
               if(hab.param %in% quali) {
                 hab = hab.param
                 writeRmd(gettext("Here, the qualitative variable selected is",domain="R-FactoInvestigate"), " : *", hab, file = file, end = "*.\n", sep = "")
               } else {
                 writeRmd(gettext("The variable",domain="R-FactoInvestigate"), " *", hab.param, "* ", gettext("cannot be selected to illustrate the plane",domain="R-FactoInvestigate"), file = file, sep = "", end = ".\n")
               }
             }
             
             if(graph) {
               sample = sample(rownames(data), length(rownames(data)))
               res$call$X = res$call$X[sample,]
               
               #shuffle
               res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
               res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],] # works even if ind.sup = NULL
               
               if(ellipse) {
                 plotellipses(res, axes = dim[1]:dim[2], invisible = 'quali', select = drawn, keepvar = hab, title = gettext("Individuals factor map (PCA)",domain="R-FactoInvestigate"), cex = cex)
               } else {
                 plot.PCA(res, axes = dim[1]:dim[2], choix = 'ind', invisible = 'quali', select = drawn, habillage = hab, title = gettext("Individuals factor map (PCA)",domain="R-FactoInvestigate"), cex = cex)
               }
             }
             
             writeRmd(file = file)
             writeRmd("sample = sample(rownames(res$call$X), length(rownames(res$call$X)))", file = file, start = TRUE, options = options)
             writeRmd("res$call$X = res$call$X[sample,]", file = file)
             writeRmd("res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]", file = file)
             writeRmd("res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]", file = file)
             dump("drawn", file = file, append = TRUE)
             dump("hab", file = file, append = TRUE)
             
             if(ellipse){
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplotellipses(res, axes = ", dim[1], ":", dim[2], ", invisible = 'quali', select = drawn, keepvar = hab, title = '', cex = cex)", file = file, sep = "", stop = TRUE, end = "\n\n")
             } else {
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res, axes = ", dim[1], ":", dim[2], ", choix = 'ind', invisible = 'quali', select = drawn, habillage = hab, title = '', cex = cex)", file = file, sep = "", stop = TRUE, end = "\n\n")
             }
             
             writeRmd("**", figure.title, " - ", gettext("Individuals factor map (PCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             writeRmd(what.drawn, file = file, sep = "")
             if(!is.null(param$ind.sup)) {
               writeRmd("*", gettext("The italic individuals represented by an empty circle are the illustrative ones, those represented by a point are the active ones",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
             }
             writeRmd("*", gettext("The individuals are coloured after their category for the variable",domain="R-FactoInvestigate"), "* ", hab, ".", sep = "", file = file)
           },
           
           CA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Rselec, coef = Rcoef)
             r.drawn = selec.res[[1]]
             r.what.drawn = selec.res[[2]]
             
             selec.res = selection(res, dim = dim, margin = 2, selec = Cselec, coef = Ccoef)
             c.drawn = selec.res[[1]]
             c.what.drawn = selec.res[[2]]
             
             
             factors = sapply(rownames(res$quali.sup$eta2), function(x) {length(levels(data[, x]))})
             
             quali = names(factors)[factors > 1]
             reject = names(factors)[factors <= 1]
             
             wilks.p = sapply(quali, function(x, res, dim) {test.de.Wilks(rbind(res$row$coord[, dim[1]:dim[2]], res$row.sup$coord[, dim[1]:dim[2]]), res$call$Xtot[, x])}, res = res, dim = dim)
             names(wilks.p) = quali
             wilks.p = sort(wilks.p)
             
             hab = names(which(wilks.p == min(wilks.p))) # maybe more than 1
             p.value = min(wilks.p)
             
             writeRmd(gettext("The Wilks test p-value indicates which variable factors are the best separated on the plane",domain="R-FactoInvestigate"), " (",
                      gettext("i.e. which one explain the best the distance between individuals",domain="R-FactoInvestigate"), ")", sep = "", file = file, end = ".\n")
             
             # wilks.s = NULL
             # if(length(hab) > 1) {
               # wilks.s = sapply(names(wilks.p[wilks.p == 0]), function(x, res, dim) {test.de.Wilks(rbind(res$row$coord[, dim[1]:dim[2]], res$row.sup$coord[, dim[1]:dim[2]]), res$call$Xtot[, x])$statistic}, res = res, dim = dim)
               # names(wilks.s) = quali[wilks.p == 0]
               # wilks.s = sort(wilks.s, decreasing = TRUE)
               
               # if(length(wilks.s) > 12) {wilks.s = wilks.s[1:12]}
               
               # hab = names(which.max(wilks.s))
             # }
             
             if(length(wilks.p) > 12) {wilks.p = wilks.p[1:12]}
             
             
             if(graph) {
               show(wilks.p)
             }
             writeRmd(start = TRUE, end = "", options = options, file = file)
             dump("wilks.p", file = file, append = TRUE)
             writeRmd("wilks.p", stop = TRUE, sep = "", file = file)
             
             
             if(length(quali) == 1) {
               writeRmd(gettext("There only is one possible qualitative variable to illustrate the distance between individuals",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
             } else {
               if(length(names(which(wilks.p == min(wilks.p)))) == 1) {
                 writeRmd(gettext("The best qualitative variable to illustrate the distance between individuals on this plane is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               } else {
                 writeRmd(gettext("Many qualitative variables have a Wilks p-value equal to zero",domain="R-FactoInvestigate"), ". ", gettext("To arbitrate which one to select, we need to compare their statistic value",domain="R-FactoInvestigate"), 
                          " : *", paste(names(which(wilks.p == min(wilks.p))), collapse = "*, *"), file = file, sep = "", end = "*.\n")
                 
                 # if(graph) {
                   # show(wilks.s)
                 # }
                 # writeRmd(start = TRUE, end = "", options = options, file = file)
                 # dump("wilks.s", file = file, append = TRUE)
                 # writeRmd("wilks.s", stop = TRUE, sep = "", file = file)
                 
                 # writeRmd(gettext("The best qualitative variable to illustrate the distance between individuals on this plane is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               }
             }
             
             if(length(reject) != 0) {
               writeRmd(gettext("The qualitative variables",domain="R-FactoInvestigate"), " *", paste(reject, collapse = ", "), "* ",
                        gettext("cannot separate the individuals on the plane, cause they are unimodal",domain="R-FactoInvestigate"), end = ".\n", file = file, sep = "")
             }
             
             if(!is.null(hab.param)) {
               if(hab.param %in% quali) {
                 hab = hab.param
                 writeRmd(gettext("Here, the qualitative variable selected is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               } else {
                 writeRmd(gettext("The variable",domain="R-FactoInvestigate"), " *", hab.param, "* ", gettext("cannot be selected to illustrate the plane",domain="R-FactoInvestigate"), file = file, sep = "", end = ".\n")
               }
             }
             
             if(graph) {
               sample = sample(rownames(data), length(rownames(data)))
               res$call$X = res$call$X[sample,]
               
               #shuffle
               res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
               res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],] # works even if ind.sup = NULL
               
               plot.CA(res, axes = dim[1]:dim[2], choix = 'CA', invisible = c('var', 'quali'), selectRow = r.drawn, selectCol = c.drawn, habillage = hab, title = gettext("Overlayed factor map (CA)",domain="R-FactoInvestigate"), cex = cex)
             }
             
             writeRmd(file = file)
             writeRmd("sample = sample(rownames(res$call$Xtot), length(rownames(res$call$Xtot)))", file = file, start = TRUE, options = options)
             writeRmd("res$call$Xtot = res$call$Xtot[sample,]", file = file)
             writeRmd("res$row$coord = res$ind$coord[sample[!sample %in% rownames(res$row.sup$coord)],]", file = file)
             writeRmd("res$row.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$row.sup$coord)],]", file = file)
             dump("r.drawn", file = file, append = TRUE)
             dump("c.drawn", file = file, append = TRUE)
             dump("hab", file = file, append = TRUE)
             
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.CA(res, axes = ", dim[1], ":", dim[2], ", choix = 'CA', invisible = c('var', 'quali'), selectRow = r.drawn, selectCol = c.drawn, habillage = hab, title = '', cex = cex)", file = file, sep = "", stop = TRUE, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Overlayed factor map (CA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             writeRmd(r.what.drawn, file = file, sep = "")
             writeRmd(c.what.drawn, file = file, sep = "")
             if(!is.null(param$ind.sup)) {
               writeRmd("*", gettext("The italic individuals represented by an empty circle are the illustrative ones, those represented by a point are the active ones",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
             }
             writeRmd("*", gettext("The individuals are coloured after their category for the variable",domain="R-FactoInvestigate"), "* ", hab, ".", sep = "", file = file)
           },
           
           CaGalt = {},
           
           MCA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Iselec, coef = Icoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             
             factors = sapply(rownames(res$quali.sup$eta2), function(x) {length(levels(data[, x]))})
             
             quali = names(factors)[factors > 1]
             reject = names(factors)[factors <= 1]
             
             wilks.p = sapply(quali, function(x, res, dim) {test.de.Wilks(rbind(res$ind$coord[, dim[1]:dim[2]], res$ind.sup$coord[, dim[1]:dim[2]]), res$call$X[, x])}, res = res, dim = dim)
             names(wilks.p) = quali
             wilks.p = sort(wilks.p)
             
             hab = names(which(wilks.p == min(wilks.p))) # maybe more than 1
             p.value = min(wilks.p)
             
             writeRmd(gettext("The Wilks test p-value indicates which variable factors are the best separated on the plane",domain="R-FactoInvestigate"), " (",
                      gettext("i.e. which one explain the best the distance between individuals",domain="R-FactoInvestigate"), ")", sep = "", file = file, end = ".\n")
             
             # wilks.s = NULL
             # if(length(hab) > 1) {
               # wilks.s = sapply(names(wilks.p[wilks.p == 0]), function(x, res, dim) {test.de.Wilks(rbind(res$ind$coord[, dim[1]:dim[2]], res$ind.sup$coord[, dim[1]:dim[2]]), res$call$X[, x])$statistic}, res = res, dim = dim)
               # names(wilks.s) = quali[wilks.p == 0]
               # wilks.s = sort(wilks.s, decreasing = TRUE)
               
               # if(length(wilks.s) > 12) {wilks.s = wilks.s[1:12]}
               
               # hab = names(which.max(wilks.s))
             # }
             
             if(length(wilks.p) > 12) {wilks.p = wilks.p[1:12]}
             
             
             if(graph) {
               show(wilks.p)
             }
             writeRmd(start = TRUE, end = "", options = options, file = file)
             dump("wilks.p", file = file, append = TRUE)
             writeRmd("wilks.p", stop = TRUE, sep = "", file = file)
             
             
             if(length(quali) == 1) {
               writeRmd(gettext("There only is one possible qualitative variable to illustrate the distance between individuals",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
             } else {
               if(length(names(which(wilks.p == min(wilks.p)))) == 1) {
                 writeRmd(gettext("The best qualitative variable to illustrate the distance between individuals on this plane is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               } else {
                 writeRmd(gettext("Many qualitative variables have a Wilks p-value equal to zero",domain="R-FactoInvestigate"), ". ", gettext("To arbitrate which one to select, we need to compare their statistic value",domain="R-FactoInvestigate"), 
                          " : *", paste(names(which(wilks.p == min(wilks.p))), collapse = "*, *"), file = file, sep = "", end = "*.\n")
                 
                 # if(graph) {
                   # show(wilks.s)
                 # }
                 # writeRmd(start = TRUE, end = "", options = options, file = file)
                 # dump("wilks.s", file = file, append = TRUE)
                 # writeRmd("wilks.s", stop = TRUE, sep = "", file = file)
                 
                 # writeRmd(gettext("The best qualitative variable to illustrate the distance between individuals on this plane is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               }
             }
             
             if(length(reject) != 0) {
               writeRmd(gettext("The qualitative variables",domain="R-FactoInvestigate"), " *", paste(reject, collapse = ", "), "* ",
                        gettext("cannot separate the individuals on the plane, cause they are unimodal",domain="R-FactoInvestigate"), end = ".\n", file = file, sep = "")
             }
             
             if(!is.null(hab.param)) {
               if(hab.param %in% quali) {
                 hab = hab.param
                 writeRmd(gettext("Here, the qualitative variable selected is",domain="R-FactoInvestigate"), " : *", hab, sep = "", file = file, end = "*.\n")
               } else {
                 writeRmd(gettext("The variable",domain="R-FactoInvestigate"), " *", hab.param, "* ", gettext("cannot be selected to illustrate the plane",domain="R-FactoInvestigate"), file = file, sep = "", end = ".\n")
               }
             }
             
             if(graph) {
               sample = sample(rownames(data), length(rownames(data)))
               res$call$X = res$call$X[sample,]
               
               #shuffle
               res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]
               res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],] # works even if ind.sup = NULL
               
               if(ellipse) {
                 plotellipses(res, axes = dim[1]:dim[2], invisible = c('var', 'quali'), select = drawn, keepvar = hab, title = gettext("Individuals factor map (MCA)",domain="R-FactoInvestigate"), cex = cex)
               } else {
                 plot.MCA(res, axes = dim[1]:dim[2], choix = 'ind', invisible = c('var', 'quali'), select = drawn, habillage = hab, title = gettext("Individuals factor map (MCA)",domain="R-FactoInvestigate"), cex = cex)
               }
             }
             
             writeRmd(file = file)
             writeRmd("sample = sample(rownames(res$call$X), length(rownames(res$call$X)))", file = file, start = TRUE, options = options)
             writeRmd("res$call$X = res$call$X[sample,]", file = file)
             writeRmd("res$ind$coord = res$ind$coord[sample[!sample %in% rownames(res$ind.sup$coord)],]", file = file)
             writeRmd("res$ind.sup$coord = res$ind.sup$coord[sample[sample %in% rownames(res$ind.sup$coord)],]", file = file)
             dump("drawn", file = file, append = TRUE)
             dump("hab", file = file, append = TRUE)
             
             if(ellipse){
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplotellipses(res, axes = ", dim[1], ":", dim[2], ", invisible = c('var', 'quali'), select = drawn, keepvar = hab, title = '', cex = cex)", file = file, sep = "", stop = TRUE, end = "\n\n")
             } else {
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res, axes = ", dim[1], ":", dim[2], ", choix = 'ind', invisible = c('var', 'quali'), select = drawn, habillage = hab, title = '', cex = cex)", file = file, sep = "", stop = TRUE, end = "\n\n")
             }
             
             writeRmd("**", figure.title, " - ", gettext("Individuals factor map (MCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             writeRmd(what.drawn, file = file, sep = "")
             if(!is.null(param$ind.sup)) {
               writeRmd("*", gettext("The italic individuals represented by an empty circle are the illustrative ones, those represented by a point are the active ones",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
             }
             writeRmd("*", gettext("The individuals are coloured after their category for the variable",domain="R-FactoInvestigate"), "* ", hab, ".", sep = "", file = file)
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
  }
