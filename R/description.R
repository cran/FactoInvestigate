description <-
function(res, file = "", dim = 1:2, desc = dim, Iselec = "contrib", Vselec = "cos2", Rselec = "cos2", Cselec = "cos2", Icoef = 1, Vcoef = 1, Rcoef = 1, Ccoef = 1, mmax = 10, nmax = 10) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(any(!desc %in% dim)) {return(warning("all the dimensions to describe must belong to the plane"))}
    
    if(!is.numeric(Iselec) & !is.character(Iselec)) {return(warning("the argument 'Iselec' should be a numeric or character vector"))}
    if(!is.numeric(Vselec) & !is.character(Vselec)) {return(warning("the argument 'Vselec' should be a numeric or character vector"))}
    if(!is.numeric(Rselec) & !is.character(Rselec)) {return(warning("the argument 'Rselec' should be a numeric or character vector"))}
    if(!is.numeric(Cselec) & !is.character(Cselec)) {return(warning("the argument 'Cselec' should be a numeric or character vector"))}
    
    if(!is.numeric(Icoef)) {return(warning("the argument 'Icoef' must be numeric"))}
    if(!is.numeric(Vcoef)) {return(warning("the argument 'Vcoef' must be numeric"))}
    if(!is.numeric(Rcoef)) {return(warning("the argument 'Rcoef' must be numeric"))}
    if(!is.numeric(Ccoef)) {return(warning("the argument 'Ccoef' must be numeric"))}
    
    if(Icoef < 0) {return(warning("the argument 'Icoef' must be positive"))}
    if(Vcoef < 0) {return(warning("the argument 'Vcoef' must be positive"))}
    if(Rcoef < 0) {return(warning("the argument 'Rcoef' must be positive"))}
    if(Ccoef < 0) {return(warning("the argument 'Ccoef' must be positive"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) | length(dim) != 2) {return(warning("the argument 'dim' has to be a 2 dimensional numeric vector"))}
    if(any(dim < 0)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    data = param$data
    switch(analyse,
           PCA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Iselec, coef = Icoef)
             drawn = selec.res[[1]]
             
             
             ind.cos2 = apply(rbind(res$ind$cos2[,dim], res$ind.sup$cos2[,dim]), 1, sum)
             selec.ind = names(ind.cos2[ind.cos2 > mean(ind.cos2)])
             
             # if(length(selec.ind) > 100) {
             #   selec.ind = names(sort(ind.cos2, decreasing = TRUE))[1:100]
             # }
             Idata = data.frame(rbind(res$ind$coord[,dim], res$ind.sup$coord[,dim]))
             if(length(selec.ind) > 100) {
               ind.hcpc = HCPC(data.frame(Idata[selec.ind,]), kk = 100, consol = FALSE, graph = FALSE)
               Idata$clust = NA
               Idata[selec.ind, "clust"] = ind.hcpc$data.clust$clust
               last.clust = length(levels(ind.hcpc$data.clust$clust)) + 1
               Idata[!rownames(Idata) %in% selec.ind, "clust"] = last.clust
             } else if(length(selec.ind) < 10) {
               ind.hcpc = HCPC(data.frame(Idata), graph = FALSE)
               Idata$clust = ind.hcpc$data.clust$clust
               last.clust = length(levels(ind.hcpc$data.clust$clust)) + 1
             } else {
               ind.hcpc = HCPC(data.frame(Idata[selec.ind,]), graph = FALSE)
               Idata$clust = NA
               Idata[selec.ind, "clust"] = ind.hcpc$data.clust$clust
               last.clust = length(levels(ind.hcpc$data.clust$clust)) + 1
               Idata[!rownames(Idata) %in% selec.ind, "clust"] = last.clust
             }
             
             
             Idata$clust = as.factor(Idata$clust)
             CD.dim = catdes(Idata[Idata$clust != last.clust,], 3, proba = 0.15)
             Itest = sapply(CD.dim$quanti, is.null)
             
             if(any(Itest)) { # identification des clusters non-caracteristiques du plan
               false.clust = names(Itest)[Itest]
               Idata$clust[Idata$clust %in% false.clust] = last.clust # fusion avec le cluster moyen
               Idata$clust = factor(Idata$clust, labels = 1:(last.clust - length(false.clust)))
               last.clust = which.max(levels(Idata$clust))
               CD.dim = catdes(Idata[Idata$clust != last.clust,], 3, proba = 0.15)
             }
             
             CD.var = catdes(cbind(data, Idata$clust)[Idata$clust != last.clust,], ncol(data) + 1)
             
             Idata = Idata[order(ind.cos2, decreasing = TRUE),] # tri par qualite de projection decroissante sur le plan (utile pour selection des 5 meilleurs individus ensuite)
             
             
             
             for(d in desc) {
               writeRmd("\n* * *", file = file, end = "\n\n")
               
               pos.groups = which(sapply(sapply(CD.dim$quanti, function(x, d) {x[rownames(x) %in% paste("Dim.", d, sep = ""), 1]}, d = d), function(x) x %dim0% 0) > 0) # groupes significativements positifs
               neg.groups = which(sapply(sapply(CD.dim$quanti, function(x, d) {x[rownames(x) %in% paste("Dim.", d, sep = ""), 1]}, d = d), function(x) x %dim0% 0) < 0) # groupes significativements negatifs
               
               ind.pos = rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]
               ind.neg = rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]
               
               if(length(ind.pos) > mmax) {
                 ind.pos = paste(paste(ind.pos[1:(mmax - 1)], collapse = "*, *"), ind.pos[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                 ind.neg = paste(paste(ind.neg[1:(mmax - 1)], collapse = "*, *"), ind.neg[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
               } else if(length(ind.pos) > 1) {
                 ind.pos = paste(paste(ind.pos[- length(ind.pos)], collapse = "*, *"), ind.pos[length(ind.pos)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                 ind.neg = paste(paste(ind.neg[- length(ind.neg)], collapse = "*, *"), ind.neg[length(ind.neg)], sep = gettext("* and *",domain="R-FactoInvestigate"))
               }
               
               if(d %% 2 == 1) {
                 if(length(rownames(Idata)[Idata$clust %in% pos.groups]) != 0) {
                   if(length(rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   }
                 } else {
                   if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("does not sufficiently discriminate individuals",domain="R-FactoInvestigate"), end = ".\n\n", file = file, sep = "")
                   }
                 }
               } else {
                 if(length(rownames(Idata)[Idata$clust %in% pos.groups]) != 0) {
                   if(length(rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   }
                 } else {
                   if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("does not sufficiently discriminate individuals",domain="R-FactoInvestigate"), end = ".\n\n", file = file, sep = "")
                   }
                 }
               }
               
               compteur = 0
               for(i in c(pos.groups, neg.groups)) { # pour chaque cluster dans le groupe
                 compteur = compteur + 1
                 if(length(c(pos.groups, neg.groups)) == 1) {
                   writeRmd(gettext("These individuals form a group sharing",domain="R-FactoInvestigate"), end = " :\n", file = file)
                 } else {
                   if(length(rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]) == 0) {
                     writeRmd(gettext("The group",domain="R-FactoInvestigate"), compteur, end = " ", file = file)
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   } else if(length(rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]) == 1) {
                     writeRmd(gettext("The group in which the individual",domain="R-FactoInvestigate"), " *", rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn],"* ", gettext("stands",domain="R-FactoInvestigate"), end = " ", file = file, sep = "")
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   } else {
                     ind.clust = rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]
                     if(length(ind.clust) > mmax) {
                       ind.clust = paste(paste(ind.clust[1:(mmax - 1)], collapse = "*, *"), ind.clust[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                     } else {
                       ind.clust = paste(paste(ind.clust[- length(ind.clust)], collapse = "*, *"), ind.clust[length(ind.clust)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                     }
                     writeRmd(gettext("The group in which the individuals",domain="R-FactoInvestigate"), " *", ind.clust, "* ", gettext("stand",domain="R-FactoInvestigate"), end = " ", file = file, sep = "")
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   }
                 }
                 
                 if(!is.null(CD.var$quanti[[i]])) {
                   if(any(CD.var$quanti[[i]][, 1] > 0)) {
                     if(nrow(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, 1] > 0,]) == 1) {
                       writeRmd("- ", gettext("high values for the variable",domain="R-FactoInvestigate"), " *", rownames(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, 1] > 0,]), end = "*.\n", file = file, sep = "")
                     } else {
                       high.var = rownames(CD.var$quanti[[i]][CD.var$quanti[[i]][, 1] > 0,])
                       if(length(high.var) > nmax) {
                         high.var = paste(paste(high.var[1:(nmax - 1)], collapse = "*, *"), high.var[nmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("high values for variables like",domain="R-FactoInvestigate"), " *", high.var, "* (", gettext("variables are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       } else {
                         high.var = paste(paste(high.var[- length(high.var)], collapse = "*, *"), high.var[length(high.var)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("high values for the variables",domain="R-FactoInvestigate"), " *", high.var, "* (", gettext("variables are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       }
                     }
                   }
                   if(any(CD.var$quanti[[i]][, 1] < 0)) {
                     if(nrow(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, 1] < 0,]) == 1) {
                       writeRmd("- ", gettext("low values for the variable",domain="R-FactoInvestigate"), " *", rownames(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, 1] < 0,]), end = "*.\n", file = file, sep = "") 
                     } else {
                       low.var = rownames(CD.var$quanti[[i]][CD.var$quanti[[i]][, 1] < 0,])
                       low.var = low.var[length(low.var):1]
                       if(length(low.var) > nmax) {
                         low.var = paste(paste(low.var[1:(nmax - 1)], collapse = "*, *"), low.var[nmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("low values for variables like",domain="R-FactoInvestigate"), " *", low.var, "* (", gettext("variables are sorted from the weakest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       } else {
                         low.var = paste(paste(low.var[- length(low.var)], collapse = "*, *"), low.var[length(low.var)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("low values for the variables",domain="R-FactoInvestigate"), " *", low.var, "* (", gettext("variables are sorted from the weakest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       }
                     }
                   }
                 } else {
                   writeRmd("- ", gettext("variables whose values do not differ significantly from the mean",domain="R-FactoInvestigate"), end = ".\n", file = file, sep = "")
                 }
                 writeRmd(file = file)
               }
               var.cos = rbind(res$var$cos2, res$quanti.sup$cos2, res$quali.sup$cos2)
               if(any(var.cos[, d] >= 0.9)) {
                 var.dim = names(which(var.cos[, d] >= 0.9))
                 if(length(which(var.cos[, d] >= 0.9)) == 1) {
                   writeRmd(gettext("Note that the variable",domain="R-FactoInvestigate"), " *", var.dim,"* ", gettext("is highly correlated with this dimension (correlation of",domain="R-FactoInvestigate"), " ", 
                            round(var.cos[var.dim, 1], 2), "). ", gettext("This variable could therefore summarize itself the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 } else if(length(which(var.cos[, d] >= 0.9)) >= nmax) {
                   writeRmd(gettext("Note that the variables",domain="R-FactoInvestigate"), " *", paste(paste(var.dim[1:(nmax - 1)], collapse = "*, *"), var.dim[nmax], sep = gettext("* and *",domain="R-FactoInvestigate")), "* ", 
                            gettext("are highly correlated with this dimension (respective correlation of",domain="R-FactoInvestigate"), " ", paste(round(var.cos[var.dim, 1], 2), collapse = ", "), "). ", 
                            gettext("These variables could therefore summarize themselve the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 } else {
                   writeRmd(gettext("Note that the variables",domain="R-FactoInvestigate"), " *", paste(paste(var.dim[1:(length(var.dim) - 1)], collapse = "*, *"), var.dim[length(var.dim)], sep = gettext("* and *",domain="R-FactoInvestigate")), "* ", 
                            gettext("are highly correlated with this dimension (respective correlation of",domain="R-FactoInvestigate"), " ", paste(round(var.cos[var.dim, 1], 2), collapse = ", "), "). ", 
                            gettext("These variables could therefore summarize themselve the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 }
               }
             }
           },
           
           CA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Rselec, coef = Rcoef)
             r.drawn = selec.res[[1]]
             
             selec.res = selection(res, dim = dim, margin = 2, selec = Cselec, coef = Ccoef)
             c.drawn = selec.res[[1]]
             
             drawn = c(r.drawn, c.drawn)
             
             
             row.cos2 = apply(rbind(res$row$cos2[,dim], res$row.sup$cos2[,dim]), 1, sum)
             selec.row = names(row.cos2[row.cos2 > mean(row.cos2)])
             
             # if(length(selec.row) > 100) {
             #   selec.row = names(sort(row.cos2, decreasing = TRUE))[1:100]
             # }
             Idata = data.frame(rbind(res$row$coord[,dim], res$row.sup$coord[,dim]))
             if(length(selec.row) > 100) {
               row.hcpc = HCPC(data.frame(Idata[selec.row,]), kk = 100, consol = FALSE, graph = FALSE)
               Idata$clust = NA
               Idata[selec.row, "clust"] = row.hcpc$data.clust$clust
               last.clust = length(levels(row.hcpc$data.clust$clust)) + 1
               Idata[!rownames(Idata) %in% selec.row, "clust"] = last.clust
             } else if(length(selec.row) < 10) {
               row.hcpc = HCPC(data.frame(Idata[selec.row,]), graph = FALSE)
               Idata$clust = row.hcpc$data.clust$clust
               last.clust = length(levels(row.hcpc$data.clust$clust)) + 1
             } else {
               row.hcpc = HCPC(data.frame(Idata[selec.row,]), graph = FALSE)
               Idata$clust = NA
               Idata[selec.row, "clust"] = row.hcpc$data.clust$clust
               last.clust = length(levels(row.hcpc$data.clust$clust)) + 1
               Idata[!rownames(Idata) %in% selec.row, "clust"] = last.clust
             }
             
             
             Idata$clust = as.factor(Idata$clust)
             CD.dim = catdes(Idata[Idata$clust != last.clust,], 3, proba = 0.15)
             Itest = sapply(CD.dim$quanti, is.null)
             
             if(any(Itest)) { # identification des clusters non-caracteristiques du plan
               false.clust = names(Itest)[Itest]
               Idata$clust[Idata$clust %in% false.clust] = last.clust # fusion avec le cluster moyen
               Idata$clust = factor(Idata$clust, labels = 1:(last.clust - length(false.clust)))
               last.clust = which.max(levels(Idata$clust))
               CD.dim = catdes(Idata[Idata$clust != last.clust,], 3, proba = 0.15)
             }
             
             CD.var = catdes(cbind(data, Idata$clust)[Idata$clust != last.clust,], ncol(data) + 1)
             
             Idata = Idata[order(row.cos2, decreasing = TRUE),] # tri par qualite de projection decroissante sur le plan (utile pour selection des 5 meilleurs individus ensuite)
             
             
             
             for(d in desc) {
               writeRmd("\n* * *", file = file, end = "\n\n")
               
               pos.groups = which(sapply(sapply(CD.dim$quanti, function(x, d) {x[rownames(x) %in% paste("Dim.", d, sep = ""), 1]}, d = d), function(x) x %dim0% 0) > 0) # groupes significativements positifs
               neg.groups = which(sapply(sapply(CD.dim$quanti, function(x, d) {x[rownames(x) %in% paste("Dim.", d, sep = ""), 1]}, d = d), function(x) x %dim0% 0) < 0) # groupes significativements negatifs
               
               ind.pos = rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]
               ind.neg = rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]
               
               if(length(ind.pos) > mmax) {
                 ind.pos = paste(paste(ind.pos[1:(mmax - 1)], collapse = "*, *"), ind.pos[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                 ind.neg = paste(paste(ind.neg[1:(mmax - 1)], collapse = "*, *"), ind.neg[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
               } else if(length(ind.pos) > 1) {
                 ind.pos = paste(paste(ind.pos[- length(ind.pos)], collapse = "*, *"), ind.pos[length(ind.pos)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                 ind.neg = paste(paste(ind.neg[- length(ind.neg)], collapse = "*, *"), ind.neg[length(ind.neg)], sep = gettext("* and *",domain="R-FactoInvestigate"))
               }
               
               if(d %% 2 == 1) {
                 if(length(rownames(Idata)[Idata$clust %in% pos.groups]) != 0) {
                   if(length(rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   }
                 } else {
                   if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("does not sufficiently discriminate factors",domain="R-FactoInvestigate"), end = ".\n\n", file = file, sep = "")
                   }
                 }
               } else {
                 if(length(rownames(Idata)[Idata$clust %in% pos.groups]) != 0) {
                   if(length(rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes factors characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to factors characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   }
                 } else {
                   if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes factors characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("does not sufficiently discriminate factors",domain="R-FactoInvestigate"), end = ".\n\n", file = file, sep = "")
                   }
                 }
               }
               
               compteur = 0
               for(i in c(pos.groups, neg.groups)) { # pour chaque cluster dans le groupe
                 compteur = compteur + 1
                 if(length(c(pos.groups, neg.groups)) == 1) {
                   writeRmd(gettext("These factors form a group sharing",domain="R-FactoInvestigate"), end = " :\n", file = file)
                 } else {
                   if(length(rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]) == 0) {
                     writeRmd(gettext("The group",domain="R-FactoInvestigate"), compteur, end = " ", file = file)
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   } else if(length(rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]) == 1) {
                     writeRmd(gettext("The group in which the factor",domain="R-FactoInvestigate"), " *", rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn],"* ", gettext("stands",domain="R-FactoInvestigate"), end = " ", file = file, sep = "")
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   } else {
                     ind.clust = rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]
                     if(length(ind.clust) > mmax) {
                       ind.clust = paste(paste(ind.clust[1:(mmax - 1)], collapse = "*, *"), ind.clust[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                     } else {
                       ind.clust = paste(paste(ind.clust[- length(ind.clust)], collapse = "*, *"), ind.clust[length(ind.clust)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                     }
                     writeRmd(gettext("The group in which the factors",domain="R-FactoInvestigate"), " *", ind.clust, "* ", gettext("stand",domain="R-FactoInvestigate"), end = " ", file = file, sep = "")
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   }
                 }
                 
                 if(!is.null(CD.var$quanti[[i]])) {
                   if(any(CD.var$quanti[[i]][, "v.test"] > 0)) {
                     if(nrow(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, "v.test"] > 0,]) == 1) {
                       writeRmd("- ", gettext("high frequency for the variable",domain="R-FactoInvestigate"), " *", rownames(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, "v.test"] > 0,]), end = "*.\n", file = file, sep = "")
                     } else {
                       high.var = rownames(CD.var$quanti[[i]][CD.var$quanti[[i]][, "v.test"] > 0,])
                       if(length(high.var) > nmax) {
                         high.var = paste(paste(high.var[1:(nmax - 1)], collapse = "*, *"), high.var[nmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("high frequency for factors like",domain="R-FactoInvestigate"), " *", high.var, "* (", gettext("factors are sorted from the most common",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       } else {
                         high.var = paste(paste(high.var[- length(high.var)], collapse = "*, *"), high.var[length(high.var)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("high frequency for the factors",domain="R-FactoInvestigate"), " *", high.var, "* (", gettext("factors are sorted from the most common",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       }
                     }
                   }
                   if(any(CD.var$quanti[[i]][, "v.test"] < 0)) {
                     if(nrow(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, "v.test"] < 0,]) == 1) {
                       writeRmd("- ", gettext("low frequency for the variable",domain="R-FactoInvestigate"), " *", rownames(data.frame(CD.var$quanti[[i]])[CD.var$quanti[[i]][, "v.test"] < 0,]), end = "*.\n", file = file, sep = "") 
                     } else {
                       low.var = rownames(CD.var$quanti[[i]][CD.var$quanti[[i]][, "v.test"] < 0,])
                       low.var = low.var[length(low.var):1] 
                       if(length(low.var) > nmax) {
                         low.var = paste(paste(low.var[1:(nmax - 1)], collapse = "*, *"), low.var[nmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("low frequency for factors like",domain="R-FactoInvestigate"), " *", low.var, "* (", gettext("factors are sorted from the rarest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       } else {
                         low.var = paste(paste(low.var[- length(low.var)], collapse = "*, *"), low.var[length(low.var)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("low frequency for the factors",domain="R-FactoInvestigate"), " *", low.var, "* (", gettext("factors are sorted from the rarest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       }
                     }
                   }
                 } else {
                   writeRmd("- ", gettext("factors whose frequency does not differ significantly from the mean",domain="R-FactoInvestigate"), end = ".\n", file = file, sep = "")
                 }
                 writeRmd(file = file)
               }
               var.cos = rbind(res$row$cos2, res$col$cos2, res$quanti.sup$cos2, res$quali.sup$cos2)
               if(any(var.cos[, d] >= 0.9)) {
                 var.dim = names(which(var.cos[, d] >= 0.9))
                 if(length(which(var.cos[, d] >= 0.9)) == 1) {
                   writeRmd(gettext("Note that the factor",domain="R-FactoInvestigate"), " *", var.dim,"* ", gettext("is highly correlated with the dimension (correlation of",domain="R-FactoInvestigate"), " ", 
                            round(var.cos[var.dim, 1], 2), "). ", gettext("This factor could therefore summarize itself the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 } else if(length(which(var.cos[, d] >= 0.9)) >= nmax) {
                   writeRmd(gettext("Note that the factors",domain="R-FactoInvestigate"), " *", paste(paste(var.dim[1:(nmax - 1)], collapse = "*, *"), var.dim[nmax], sep = gettext("* and *",domain="R-FactoInvestigate")), "* ", 
                            gettext("are highly correlated with the dimension (respective correlation of",domain="R-FactoInvestigate"), " ", paste(round(var.cos[var.dim, 1], 2), collapse = ", "), "). ", 
                            gettext("These factors could therefore summarize themselve the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 } else {
                   writeRmd(gettext("Note that the factors",domain="R-FactoInvestigate"), " *", paste(paste(var.dim[1:(length(var.dim) - 1)], collapse = "*, *"), var.dim[length(var.dim)], sep = gettext("* and *",domain="R-FactoInvestigate")), "* ", 
                            gettext("are highly correlated with the dimension (respective correlation of",domain="R-FactoInvestigate"), " ", paste(round(var.cos[var.dim, 1], 2), collapse = ", "), "). ", 
                            gettext("These factors could therefore summarize themselve the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 }
               }
             }
           },
           
           CaGalt = {},
           
           MCA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Iselec, coef = Icoef)
             drawn = selec.res[[1]]
             
             ind.cos2 = apply(rbind(res$ind$cos2[,dim], res$ind.sup$cos2[,dim]), 1, sum)
             selec.ind = names(ind.cos2[ind.cos2 > mean(ind.cos2)])
             
             # if(length(selec.ind) > 100) {
             #   selec.ind = names(sort(ind.cos2, decreasing = TRUE))[1:100]
             # }
             Idata = data.frame(rbind(res$ind$coord[,dim], res$ind.sup$coord[,dim]))
             if(length(selec.ind) > 100) {
               ind.hcpc = HCPC(data.frame(Idata[selec.ind,]), kk = 100, consol = FALSE, graph = FALSE)  ## pb need that more than 100 points are distinct
               Idata$clust = NA
               Idata[selec.ind, "clust"] = ind.hcpc$data.clust$clust
               last.clust = length(levels(ind.hcpc$data.clust$clust)) + 1
               Idata[!rownames(Idata) %in% selec.ind, "clust"] = last.clust
             } else if(length(selec.ind) < 10) {
               ind.hcpc = HCPC(data.frame(Idata[selec.ind,]), graph = FALSE)
               Idata$clust = ind.hcpc$data.clust$clust
               last.clust = length(levels(ind.hcpc$data.clust$clust)) + 1
             } else {
               ind.hcpc = HCPC(data.frame(Idata[selec.ind,]), graph = FALSE)
               Idata$clust = NA
               Idata[selec.ind, "clust"] = ind.hcpc$data.clust$clust
               last.clust = length(levels(ind.hcpc$data.clust$clust)) + 1
               Idata[!rownames(Idata) %in% selec.ind, "clust"] = last.clust
             }
             
             Idata$clust = as.factor(Idata$clust)
             CD.dim = catdes(Idata[Idata$clust != last.clust,], 3, proba = 0.15)
             Itest = sapply(CD.dim$quanti, is.null)
             
             if(any(Itest)) { # identification des clusters non-caracteristiques du plan
               false.clust = names(Itest)[Itest]
               Idata$clust[Idata$clust %in% false.clust] = last.clust # fusion avec le cluster moyen
               Idata$clust = factor(Idata$clust, labels = 1:(last.clust - length(false.clust)))
               last.clust = which.max(levels(Idata$clust))
               CD.dim = catdes(Idata[Idata$clust != last.clust,], 3, proba = 0.15)
             }
             
             CD.var = catdes(cbind(data, Idata$clust)[Idata$clust != last.clust,], ncol(data) + 1)
             
             Idata = Idata[order(ind.cos2, decreasing = TRUE),] # tri par qualite de projection decroissante sur le plan (utile pour selection des 5 meilleurs individus ensuite)
             
             
             
             for(d in desc) {
               writeRmd("\n* * *", file = file, end = "\n\n")
               
               pos.groups = which(sapply(sapply(CD.dim$quanti, function(x, d) {x[rownames(x) %in% paste("Dim.", d, sep = ""), 1]}, d = d), function(x) x %dim0% 0) > 0) # groupes significativements positifs
               neg.groups = which(sapply(sapply(CD.dim$quanti, function(x, d) {x[rownames(x) %in% paste("Dim.", d, sep = ""), 1]}, d = d), function(x) x %dim0% 0) < 0) # groupes significativements negatifs
               
               ind.pos = rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]
               ind.neg = rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]
               
               if(length(ind.pos) > mmax) {
                 ind.pos = paste(paste(ind.pos[1:(mmax - 1)], collapse = "*, *"), ind.pos[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                 ind.neg = paste(paste(ind.neg[1:(mmax - 1)], collapse = "*, *"), ind.neg[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
               } else if(length(ind.pos) > 1) {
                 ind.pos = paste(paste(ind.pos[- length(ind.pos)], collapse = "*, *"), ind.pos[length(ind.pos)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                 ind.neg = paste(paste(ind.neg[- length(ind.neg)], collapse = "*, *"), ind.neg[length(ind.neg)], sep = gettext("* and *",domain="R-FactoInvestigate"))
               }
               
               if(d %% 2 == 1) {
                 if(length(rownames(Idata)[Idata$clust %in% pos.groups]) != 0) {
                   if(length(rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                "* (", gettext("to the right of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the right of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   }
                 } else {
                   if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                "* (", gettext("to the left of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the left of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("does not sufficiently discriminate individuals",domain="R-FactoInvestigate"), end = ".\n\n", file = file, sep = "")
                   }
                 }
               } else {
                 if(length(rownames(Idata)[Idata$clust %in% pos.groups]) != 0) {
                   if(length(rownames(Idata)[Idata$clust %in% pos.groups & rownames(Idata) %in% drawn]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                  "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.pos, 
                                "* (", gettext("to the top of the graph, characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                       if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                  "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       } else {
                         writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("opposes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                  " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ")\n", file = file, sep = "")
                         writeRmd(gettext("to individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                       }
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly positive coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the top of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   }
                 } else {
                   if(length(rownames(Idata)[Idata$clust %in% neg.groups]) != 0) {
                     if(length(rownames(Idata)[Idata$clust %in% neg.groups & rownames(Idata) %in% drawn]) != 0) {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals such as",domain="R-FactoInvestigate"), " *", ind.neg, 
                                "* (", gettext("to the bottom of the graph, characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     } else {
                       writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("particularly distinguishes individuals characterized by a strongly negative coordinate on the axis",domain="R-FactoInvestigate"), 
                                " (", gettext("to the bottom of the graph",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                     }
                   } else {
                     writeRmd(gettext("The **dimension",domain="R-FactoInvestigate"), " ", d, "** ", gettext("does not sufficiently discriminate individuals",domain="R-FactoInvestigate"), end = ".\n\n", file = file, sep = "")
                   }
                 }
               }
               
               compteur = 0
               for(i in c(pos.groups, neg.groups)) { # pour chaque cluster dans le groupe
                 compteur = compteur + 1
                 if(length(c(pos.groups, neg.groups)) == 1) {
                   writeRmd(gettext("These individuals form a group sharing",domain="R-FactoInvestigate"), end = " :\n", file = file)
                 } else {
                   if(length(rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]) == 0) {
                     writeRmd(gettext("The group",domain="R-FactoInvestigate"), compteur, end = " ", file = file)
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   } else if(length(rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]) == 1) {
                     writeRmd(gettext("The group in which the individual",domain="R-FactoInvestigate"), " *", rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn],"* ", gettext("stands",domain="R-FactoInvestigate"), end = " ", file = file, sep = "")
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   } else {
                     ind.clust = rownames(Idata)[Idata$clust == i & rownames(Idata) %in% drawn]
                     if(length(ind.clust) > mmax) {
                       ind.clust = paste(paste(ind.clust[1:(mmax - 1)], collapse = "*, *"), ind.clust[mmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                     } else {
                       ind.clust = paste(paste(ind.clust[- length(ind.clust)], collapse = "*, *"), ind.clust[length(ind.clust)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                     }
                     writeRmd(gettext("The group in which the individuals",domain="R-FactoInvestigate"), " *", ind.clust, "* ", gettext("stand",domain="R-FactoInvestigate"), end = " ", file = file, sep = "")
                     if(i %in% pos.groups) {
                       writeRmd("(", gettext("characterized by a positive coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     } else {
                       writeRmd("(", gettext("characterized by a negative coordinate on the axis",domain="R-FactoInvestigate"), ")", end = " ", file = file, sep = "")
                     }
                     writeRmd(gettext("is sharing",domain="R-FactoInvestigate"), end = " :\n\n", file = file, sep = "")
                   }
                 }
                 
                 if(!is.null(CD.var$category[[i]])) {
                   if(any(CD.var$category[[i]][, "v.test"] > 0)) {
                     if(nrow(data.frame(CD.var$category[[i]])[CD.var$category[[i]][, "v.test"] > 0,]) == 1) {
                       writeRmd("- ", gettext("high frequency for the variable",domain="R-FactoInvestigate"), " *", rownames(data.frame(CD.var$category[[i]])[CD.var$category[[i]][, "v.test"] > 0,]), end = "*.\n", file = file, sep = "")
                     } else {
                       high.var = rownames(CD.var$category[[i]][CD.var$category[[i]][, "v.test"] > 0,])
                       if(length(high.var) > nmax) {
                         high.var = paste(paste(high.var[1:(nmax - 1)], collapse = "*, *"), high.var[nmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("high frequency for factors like",domain="R-FactoInvestigate"), " *", high.var, "* (", gettext("factors are sorted from the most common",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       } else {
                         high.var = paste(paste(high.var[- length(high.var)], collapse = "*, *"), high.var[length(high.var)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("high frequency for the factors",domain="R-FactoInvestigate"), " *", high.var, "* (", gettext("factors are sorted from the most common",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       }
                     }
                   }
                   if(any(CD.var$category[[i]][, "v.test"] < 0)) {
                     if(nrow(data.frame(CD.var$category[[i]])[CD.var$category[[i]][, "v.test"] < 0,]) == 1) {
                       writeRmd("- ", gettext("low frequency for the variable",domain="R-FactoInvestigate"), " *", rownames(data.frame(CD.var$category[[i]])[CD.var$category[[i]][, "v.test"] < 0,]), end = "*.\n", file = file, sep = "") 
                     } else {
                       low.var = rownames(CD.var$category[[i]][CD.var$category[[i]][, "v.test"] < 0,])
                       low.var = low.var[length(low.var):1] 
                       if(length(low.var) > nmax) {
                         low.var = paste(paste(low.var[1:(nmax - 1)], collapse = "*, *"), low.var[nmax], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("low frequency for factors like",domain="R-FactoInvestigate"), " *", low.var, "* (", gettext("factors are sorted from the rarest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       } else {
                         low.var = paste(paste(low.var[- length(low.var)], collapse = "*, *"), low.var[length(low.var)], sep = gettext("* and *",domain="R-FactoInvestigate"))
                         writeRmd("- ", gettext("low frequency for the factors",domain="R-FactoInvestigate"), " *", low.var, "* (", gettext("factors are sorted from the rarest",domain="R-FactoInvestigate"), end = ").\n", file = file, sep = "")
                       }
                     }
                   }
                 } else {
                   writeRmd("- ", gettext("factors whose frequency does not differ significantly from the mean",domain="R-FactoInvestigate"), end = ".\n", file = file, sep = "")
                 }
                 writeRmd(file = file)
               }
               var.cos = rbind(res$var$cos2, res$quanti.sup$cos2, res$quali.sup$cos2)
               if(any(var.cos[, d] >= 0.9)) {
                 var.dim = names(which(var.cos[, d] >= 0.9))
                 if(length(which(var.cos[, d] >= 0.9)) == 1) {
                   writeRmd(gettext("Note that the factor",domain="R-FactoInvestigate"), " *", var.dim,"* ", gettext("is highly correlated with the dimension (correlation of",domain="R-FactoInvestigate"), " ", 
                            round(var.cos[var.dim, 1], 2), "). ", gettext("This factor could therefore summarize itself the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 } else if(length(which(var.cos[, d] >= 0.9)) >= nmax) {
                   writeRmd(gettext("Note that the factors",domain="R-FactoInvestigate"), " *", paste(paste(var.dim[1:(nmax - 1)], collapse = "*, *"), var.dim[nmax], sep = gettext("* and *",domain="R-FactoInvestigate")), "* ", 
                            gettext("are highly correlated with the dimension (respective correlation of",domain="R-FactoInvestigate"), " ", paste(round(var.cos[var.dim, 1], 2), collapse = ", "), "). ", 
                            gettext("These factors could therefore summarize themselve the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 } else {
                   writeRmd(gettext("Note that the factors",domain="R-FactoInvestigate"), " *", paste(paste(var.dim[1:(length(var.dim) - 1)], collapse = "*, *"), var.dim[length(var.dim)], sep = gettext("* and *",domain="R-FactoInvestigate")), "* ", 
                            gettext("are highly correlated with the dimension (respective correlation of",domain="R-FactoInvestigate"), " ", paste(round(var.cos[var.dim, 1], 2), collapse = ", "), "). ", 
                            gettext("These factors could therefore summarize themselve the dimension",domain="R-FactoInvestigate"), " ", d, end = ".\n", file = file, sep = "")
                 }
               }
             }
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
  }
