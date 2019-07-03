outliers <-
function(res, file = "", Vselec = "cos2", Vcoef = 1, nmax = 10, figure.title = "Figure", graph = TRUE, cex = 0.7, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(nmax)) {return(warning("the argument 'nmax' must be numeric"))}
    if(nmax < 0) {return(warning("the argument 'nmax' must be positive"))}
    if(!is.numeric(Vcoef)) {return(warning("the argument 'Vcoef' must be numeric"))}
    if(Vcoef < 0) {return(warning("the argument 'Vcoef' must be positive"))}
    if(!is.numeric(cex)) {return(warning("the argument 'cex' must be numeric"))}
    if(cex < 0) {return(warning("the argument 'cex' must be positive"))}
    
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC")) {
      return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))
    }
    else if(analyse == 'CA') {
      writeRmd(gettext("The detection of outliers does not apply to CA results",domain="R-FactoInvestigate"), file = file, end = ".\n\n")
      return(list(new.res = res, res.out = res, memory = NULL))
    } else if(analyse == "CaGalt") {
      writeRmd(gettext("The detection of outliers does not apply to CaGalt results",domain="R-FactoInvestigate"), file = file, end = ".\n\n")
      return(list(new.res = res, res.out = res, memory = NULL))
    }
    param = getParam(res)
    
    memory = res
    data = param$data
    ind.names = rownames(param$data)
    ind.sup = param$ind.sup
    extrem = NULL
    continue = TRUE
    
    while(any(scale(apply(res$ind$contrib[, 1:2], 1, function(x, coeff) {weighted.mean(x, coeff)}, 
                          coeff = res$eig[1:2, 1])) > 3) & continue) {
res.temp = res
      test.ind = scale(apply(res$ind$contrib[, 1:2], 1, function(x, coeff) {
        weighted.mean(x, coeff)
      }, coeff = res$eig[1:2, 1]))
      names = rownames(test.ind)[which(test.ind > 3)] # noms des individus identifies exceptionnels
      actual.names = rownames(res$ind$coord) # noms des individus de l'ACP actuelle
      new.sup = which(ind.names %in% names) # numeros de ces individus dans l'ACP d'origine
      actual.new.sup = which(actual.names %in% names)
      if (!is.null(res$call$row.w.init)) row.w = res$call$row.w.init[-actual.new.sup] # correction de la longueur du vecteur poids pour la nouvelle ACP
      else row.w = res$call$row.w[-actual.new.sup] # correction de la longueur du vecteur poids pour la nouvelle ACP
      ind.sup = c(ind.sup, new.sup)
      
      switch(analyse,
             # PCA = {
               # res = PCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, 
                         # graph = FALSE, scale.unit = param$scale, row.w = row.w, col.w = param$col.w, ncp = param$ncp.mod)
             # },
             
             # MCA = {
               # res = MCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, 
                         # graph = FALSE, row.w = row.w, ncp = param$ncp.mod)
             # },
             
             PCA = {
               if (is.null(c(param$quanti.sup, param$quali.sup))) res = PCA(param$data[-ind.sup,], graph = FALSE, scale.unit = param$scale, row.w = row.w, col.w = param$col.w, ncp = param$ncp.mod)
               else res = PCA(param$data[-ind.sup,-c(param$quanti.sup, param$quali.sup)], graph = FALSE, scale.unit = param$scale, row.w = row.w, col.w = param$col.w, ncp = param$ncp.mod)
             },
             
             MCA = {
               if (is.null(c(param$quanti.sup, param$quali.sup))) res = MCA(droplevels(param$data[-ind.sup,]), graph = FALSE, row.w = row.w, ncp = param$ncp.mod)
               else res = MCA(droplevels(param$data[-ind.sup,-c(param$quanti.sup, param$quali.sup)]), graph = FALSE, row.w = row.w, ncp = param$ncp.mod)
             },

             MFA = {},
             
             HMFA = {},
             
             DMFA = {},
             
             FAMD = {},
             
             GPA = {},
             
             HCPC = {})
      
      COR = cor(res.temp$ind$coord[-actual.new.sup, ], res$ind$coord)
      if(weighted.mean(abs(c(COR[1,1], COR[2,2])), res.temp$eig[1:2, 1]) > 0.9){ # si on considere que l'individu est tres particulier, mais pas anormal
        res = res.temp
        ind.sup = res$call$ind.sup %dim0% NULL # on recupere la derniere selection validee
        continue = FALSE
      } else {
        extrem = c(extrem, new.sup)
      }
    }
  
    if(!is.null(extrem)) {
      if(length(extrem) == 1){
        writeRmd(gettext("The analysis of the graphs leads to detect an outlier that strongly influences the results",,domain="R-FactoInvestigate"),
                 gettext("First we will describe this outlier and then we will suppress it from the analysis",domain="R-FactoInvestigate"), file = file, sep = ". ", end = ".\n")
        writeRmd(gettext("Looking at the graph, we can note that a particular individual strongly contributes to the construction of the plane",domain="R-FactoInvestigate"), ". ",
                 gettext("Its contribution to the construction of the plane equals",domain="R-FactoInvestigate"), " **", round(weighted.mean(memory$ind$contrib[extrem, 1:2], memory$eig[1:2,1]), 1), 
                 end = "%**.\n", file = file, sep = "")
      } else {
        writeRmd(gettext("The analysis of the graphs leads to detect outliers that strongly influence the results",domain="R-FactoInvestigate"),
                 gettext("First we will describe these outliers and then we will suppress them from the analysis",domain="R-FactoInvestigate"), file = file, sep = ". ", end = ".\n")
        writeRmd(gettext("Looking at the graph, we can note that",domain="R-FactoInvestigate"), " ", length(extrem), " ", gettext("particular individuals strongly contribute to the construction of the plane",domain="R-FactoInvestigate"), ". ",
                 gettext("The cumulative contribution of these individuals to the construction of the plane equals",domain="R-FactoInvestigate"), " **", round(weighted.mean(apply(memory$ind$contrib[extrem, 1:2], 2, sum), memory$eig[1:2,1]), 1), 
                 end = "%**.\n", file = file, sep = "")
      }
      
#      if(length(which(unlist(lapply(param$data[- extrem, which(!sapply(param$data, is.numeric)), drop=FALSE], summary)) == 0)) > 0) {  # pb un individu estreme est seul a prendre une modalite
        param$data <- param$data[- extrem, ]
        if (!is.null(memory$call$row.w.init)) row.w <- memory$call$row.w.init[- extrem]
		    else row.w <- memory$call$row.w[- extrem]
        if(!is.null(param$ind.sup)){
          ind.sup = which(rownames(param$data) %in% ind.names[extrem])
        } else {
          ind.sup = NULL
        }
        drawn = integer(0)
#       } else {
#         ind.sup = c(param$ind.sup, extrem)
#         if (!is.null(param$row.w.init)) row.w = param$row.w.init[-which(ind.names %in% ind.names[extrem])]
# 		    else row.w = param$row.w[-which(ind.names %in% ind.names[extrem])]
#         drawn = ind.names[extrem]
#       }
      
		switch(analyse,
             PCA = {
               res = PCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, graph = FALSE, scale.unit = param$scale, row.w = row.w, col.w = param$col.w, ncp = param$ncp.mod)
               
               if(graph) {
                 plot.PCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = ind.names[extrem], title = gettext("Individuals factor map (PCA) before correction",domain="R-FactoInvestigate"), cex = cex)
                 plot.PCA(res, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = gettext("Individuals factor map (PCA) after correction",domain="R-FactoInvestigate"), cex = cex)
               }
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 1, sep = "."), " - ", gettext("Individuals factor map (PCA) before correction",domain="R-FactoInvestigate"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of",domain="R-FactoInvestigate"), " ", length(extrem), " ", gettext("outliers",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               }
               
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 2, sep = "."), " - ", gettext("Individuals factor map (PCA) after correction",domain="R-FactoInvestigate"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of",domain="R-FactoInvestigate"), " ", length(extrem), " ", gettext("outliers",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               }
               
               selec.res = selection(res, dim = 1:2, margin = 2, selec = Vselec, coef = Vcoef)
               drawn = selec.res[[1]]
               what.drawn = selec.res[[2]]
               
               if(graph) {
                 plot.PCA(memory, choix = 'var', select = drawn, title = gettext("Variables factor map (PCA) before correction",domain="R-FactoInvestigate"), cex = cex)
                 plot.PCA(res, choix = 'var', select = drawn, title = gettext("Variables factor map (PCA) after correction",domain="R-FactoInvestigate"), cex = cex)
               }
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(memory, choix = 'var', select = drawn, title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
               writeRmd("**", paste(figure.title, 3, sep = "."), " - ", gettext("Variables factor map (PCA) before correction",domain="R-FactoInvestigate"), "**", file = file, sep = "")
               if(!is.null(param$quali.sup)) {
                 writeRmd("*", gettext("The variables in black are considered as active whereas those in blue are illustrative",domain="R-FactoInvestigate"), ".*", file = file, sep = "")
               }
               writeRmd(what.drawn, file = file, sep = "")
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res, choix = 'var', select = drawn, title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
               writeRmd("**", paste(figure.title, 3, sep = "."), " - ", gettext("Variables factor map (PCA) after correction",domain="R-FactoInvestigate"), "**", file = file, sep = "")
               
               
               writeRmd("\n- - -", file = file, end = "\n\n")
               
               if(is.null(c(param$quali.sup, param$quanti.sup))){
                 if(is.null(param$ind.sup)){
                   centred = scale(data)
                   rownames(centred)=rownames(data)
                 } else {
                   centred = scale(data[-param$ind.sup,])
                   rownames(centred)=rownames(data[-param$ind.sup,])
                 }
			   } else {
                 if(is.null(param$ind.sup)){
                   centred = scale(data[, -c(param$quali.sup, param$quanti.sup)])
                   rownames(centred)=rownames(data)
                 } else {
                   centred = scale(data[-param$ind.sup, -c(param$quali.sup, param$quanti.sup)])
                   rownames(centred)=rownames(data[-param$ind.sup,])
                 }
               }
               
               for(i in extrem){
                 writeRmd("**", gettext("The individual",domain="R-FactoInvestigate"), " ", i, end = "** :\n\n", file = file, sep = "")
                 writeRmd("-", gettext("takes very high values for the variable(s)",domain="R-FactoInvestigate"), end = " :\n", file = file)
                 liste.pos = names(which(sort(centred[rownames(data)[i],], decreasing = TRUE) > 2))
                 if(length(liste.pos) == 1) {
                   writeRmd("*", liste.pos, end = "*.\n\n", file = file, sep = "")
                 } else if(length(liste.pos) > 1 & length(liste.pos) < nmax) {
                   writeRmd("*", paste(paste(liste.pos[- length(liste.pos)], collapse = "*, *"), liste.pos[length(liste.pos)], sep = gettext("* and *",domain="R-FactoInvestigate")), "* (", gettext("variables are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                 } else if(length(liste.pos) >= nmax) {
                   writeRmd("*", paste(liste.pos[1:(nmax - 1)], collapse = "*, *"), gettext("* and *",domain="R-FactoInvestigate"), liste.pos[nmax], "* (", gettext("variables are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                 }
                 
                 writeRmd("-", gettext("takes very low values for the variable(s)",domain="R-FactoInvestigate"), end = " :\n", file = file)
                 liste.neg = names(which(sort(centred[rownames(data)[i],], decreasing = FALSE) < -2))
                 if(length(liste.neg) == 1) {
                   writeRmd("*", liste.neg, end = "*.\n\n", file = file, sep = "")
                 } else if(length(liste.neg) > 1 & length(liste.neg) < nmax) {
                   writeRmd("*", paste(paste(liste.neg[- length(liste.neg)], collapse = "*, *"), liste.neg[length(liste.neg)], sep = gettext("* and *",domain="R-FactoInvestigate")), "* (", gettext("variables are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                 } else if(length(liste.neg) >= nmax) {
                   writeRmd("*", paste(liste.neg[1:(nmax - 1)], collapse = "*, *"), gettext("* and *",domain="R-FactoInvestigate"), liste.neg[nmax], "* (", gettext("variables are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                 }
               }
               
               if(length(extrem) == 1) {
                 writeRmd(gettext("This outlier is suppressed from the analysis and a second one is performed on the rest of the individuals",domain="R-FactoInvestigate"), end = ".\n", file = file)
               } else {
                 writeRmd(gettext("These outliers are suppressed from the analysis and a second one is performed on the rest of the individuals",domain="R-FactoInvestigate"), end = ".\n", file = file)
               }
             },
             
             MCA = {
               res = MCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, graph = FALSE, row.w = row.w, ncp = param$ncp.mod)
               
               if(graph) {
                 plot.MCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = ind.names[extrem], title = gettext("Individuals factor map (MCA) before correction",domain="R-FactoInvestigate"), cex = cex)
                 plot.MCA(res, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = gettext("Individuals factor map (MCA) after correction",domain="R-FactoInvestigate"), cex = cex)
               }
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 1, sep = "."), " - ", gettext("Individuals factor map (MCA) before correction",domain="R-FactoInvestigate"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of",domain="R-FactoInvestigate"), " ", length(extrem), " ", gettext("outliers",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               }

               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 2, sep = "."), " - ", gettext("Individuals factor map (MCA) after correction",domain="R-FactoInvestigate"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of",domain="R-FactoInvestigate"), " ", length(extrem), " ", gettext("outliers",domain="R-FactoInvestigate"), end = ".* \n", file = file, sep = "")
               }
               
               selec.res = selection(res, dim = 1:2, margin = 2, selec = Vselec, coef = Vcoef)
               drawn = selec.res[[1]]
               what.drawn = selec.res[[2]]
               
               if(graph) {
                 plot.MCA(memory, choix = 'ind', invisible = 'ind', selectMod = drawn, title = gettext("Variables factor map (MCA) before correction",domain="R-FactoInvestigate"), cex = cex)
                 plot.MCA(res, choix = 'ind', invisible = 'ind', selectMod = drawn, title = gettext("Variables factor map (MCA) after correction",domain="R-FactoInvestigate"), cex = cex)
               }
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(memory, choix = 'ind', invisible = 'ind', selectMod = drawn, title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
               writeRmd("**", paste(figure.title, 3, sep = "."), " - ", gettext("Variables factor map (MCA) before correction",domain="R-FactoInvestigate"), "**", file = file, sep = "")
               if(!is.null(param$quali.sup)) {
                 writeRmd("*", gettext("The factors in red are considered as active whereas those in green are illustrative",domain="R-FactoInvestigate"), ".*", file = file, sep = "")
               }
               writeRmd(what.drawn, file = file, sep = "")
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = options, end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res, choix = 'ind', invisible = 'ind', selectMod = drawn, title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
               writeRmd("**", paste(figure.title, 3, sep = "."), " - ", gettext("Variables factor map (MCA) after correction",domain="R-FactoInvestigate"), "**", file = file, sep = "")
               
               
               writeRmd("\n- - -", file = file, end = "\n\n")
               
               if(is.null(c(param$quali.sup, param$quanti.sup))){
                 if(is.null(param$ind.sup)){
                   centred = scale(tab.disjonctif(data))
                   rownames(centred)=rownames(data)
                 } else {
                   centred = scale(tab.disjonctif(data[-param$ind.sup,]))
                   rownames(centred)=rownames(data[-param$ind.sup,])
                 }
               } else {
                 if(is.null(param$ind.sup)){
                   centred = scale(tab.disjonctif(data[, -c(param$quali.sup, param$quanti.sup)]))
                   rownames(centred)=rownames(data)
                 } else {
                   centred = scale(tab.disjonctif(data[-param$ind.sup, -c(param$quali.sup, param$quanti.sup)]))
                   rownames(centred)=rownames(data[-param$ind.sup,])
                 }
               }
               
               for(i in extrem){
                 writeRmd("**", gettext("The individual",domain="R-FactoInvestigate"), " ", i, end = "** :\n\n", file = file, sep = "")
                 
                 writeRmd("-", gettext("is characterized by the factor(s)",domain="R-FactoInvestigate"), end = " :\n", file = file)
                 liste.pos = names(which(sort(centred[rownames(data)[i],], decreasing = TRUE) > 2))
                 if(length(liste.pos) == 1) {
                   writeRmd("*", liste.pos, end = "*.\n\n", file = file, sep = "")
                 } else if(length(liste.pos) > 1 & length(liste.pos) < nmax) {
                   writeRmd("*", paste(liste.pos, collapse = "*, *"), "* (", gettext("factors are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                 } else if(length(liste.pos) >= nmax) {
                   writeRmd("*", paste(liste.pos[1:(nmax - 1)], collapse = "*, *"), gettext("* and *",domain="R-FactoInvestigate"), liste.pos[nmax], "* (", gettext("factors are sorted from the strongest",domain="R-FactoInvestigate"), end = ").\n\n", file = file, sep = "")
                 }
               }
               
               if(length(extrem) == 1) {
                 writeRmd(gettext("This outlier is suppressed from the analysis and a second one is performed on the rest of the individuals",domain="R-FactoInvestigate"), end = ".\n", file = file)
               } else {
                 writeRmd(gettext("These outliers are suppressed from the analysis and a second one is performed on the rest of the individuals",domain="R-FactoInvestigate"), end = ".\n", file = file)
               }
             },
             
             MFA = {},
             
             HMFA = {},
             
             DMFA = {},
             
             FAMD = {},
             
             GPA = {},
             
             HCPC = {})
      
    } else {
      writeRmd(gettext("The analysis of the graphs does not detect any outlier",domain="R-FactoInvestigate"), file = file, end = ".\n")
      res.out = res
      memory = NULL
    }
    
    writeRmd("\n- - -", file = file, end = "\n\n")
    
    list(new.res = res, memory = memory, N = length(extrem), ID = extrem)
  }
