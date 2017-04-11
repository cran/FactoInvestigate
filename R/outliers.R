outliers <-
function(res, file = "", Vselec = "cos2", Vcoef = 1, nmax = 10, figure.title = "Figure", graph = TRUE, cex = 0.7) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(nmax)) {return(warning("the argument 'nmax' must be numeric"))}
    if(nmax < 0) {return(warning("the argument 'nmax' must be positive"))}
    if(!is.numeric(Vcoef)) {return(warning("the argument 'Vcoef' must be numeric"))}
    if(Vcoef < 0) {return(warning("the argument 'Vcoef' must be positive"))}
    if(!is.numeric(cex)) {return(warning("the argument 'cex' must be numeric"))}
    if(cex < 0) {return(warning("the argument 'cex' must be positive"))}
    
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    if(analyse == 'CA') {
      writeRmd(gettext("The detection of outliers does not apply to CA results"), file = file, end = ".\n\n")
      return(list(new.res = res, res.out = res, memory = NULL))
    } else if(analyse == "CaGalt") {
      writeRmd(gettext("The detection of outliers does not apply to CaGalt results"), file = file, end = ".\n\n")
      return(list(new.res = res, res.out = res, memory = NULL))
    }
    
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
      names = rownames(test.ind)[which(test.ind > 3)] # noms des individus identifiés exceptionnels
      actual.names = rownames(res$ind$coord) # noms des individus de l'ACP actuelle
      new.sup = which(ind.names %in% names) # numéros de ces individus dans l'ACP d'origine
      actual.new.sup = which(actual.names %in% names)
      row.w = res$call$row.w[-actual.new.sup] # correction de la longueur du vecteur poids pour la nouvelle ACP
      ind.sup = c(ind.sup, new.sup)
      
      switch(analyse,
             PCA = {
               res = PCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, 
                         graph = FALSE, scale.unit = param$scale, row.w = row.w, col.w = param$col.w, ncp = param$ncp.mod)
             },
             
             MCA = {
               res = MCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, 
                         graph = FALSE, row.w = row.w, ncp = param$ncp.mod)
             },
             
             MFA = {},
             
             HMFA = {},
             
             DMFA = {},
             
             FAMD = {},
             
             GPA = {},
             
             HCPC = {})
      
      COR = cor(res.temp$ind$coord[-actual.new.sup, ], res$ind$coord)
      if(weighted.mean(abs(c(COR[1,1], COR[2,2])), res.temp$eig[1:2, 1]) > 0.9){ # si on considère que l'individu est très particulier, mais pas anormal
        res = res.temp
        ind.sup = res$call$ind.sup %dim0% NULL # on récupère la dernière sélection validée
        continue = FALSE
      } else {
        extrem = c(extrem, new.sup)
      }
    }
    
    
    if(!is.null(extrem)) {
      if(length(extrem) == 1){
        writeRmd(gettext("The analysis of the graphs leads to detect an outlier that strongly influences the results"),
                 gettext("First we will describe this outlier and then we will suppress it from the analysis"), file = file, sep = ". ", end = ".\n")
        writeRmd(gettext("Looking at the graph, we can note that a particular individual strongly contributes to the construction of the plane"), ". ",
                 gettext("Its contribution to the construction of the plane equals"), " **", round(weighted.mean(memory$ind$contrib[extrem, 1:2], memory$eig$eigenvalue[1:2]), 1), 
                 end = "%**.\n", file = file, sep = "")
      } else {
        writeRmd(gettext("The analysis of the graphs leads to detect outliers that strongly influence the results"),
                 gettext("First we will describe these outliers and then we will suppress them from the analysis"), file = file, sep = ". ", end = ".\n")
        writeRmd(gettext("Looking at the graph, we can note that"), " ", length(extrem), " ", gettext("particular individuals strongly contribute to the construction of the plane"), ". ",
                 gettext("The cumulative contribution of these individuals to the construction of the plane equals"), " **", round(weighted.mean(apply(memory$ind$contrib[extrem, 1:2], 2, sum), memory$eig$eigenvalue[1:2]), 1), 
                 end = "%**.\n", file = file, sep = "")
      }
      
      
      switch(analyse,
             PCA = {
               if(graph) {
                 plot.PCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = ind.names[extrem], title = gettext("Individuals factor map (PCA) before correction"), cex = cex)
                 plot.PCA(res, choix = 'ind', invisible = c('var', 'quali'), select = ind.names[extrem], title = gettext("Individuals factor map (PCA) after correction"), cex = cex)
               }
               
               drawn = ind.names[extrem]
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 1, sep = "."), " - ", gettext("Individuals factor map (PCA) before correction"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of"), " ", length(extrem), " ", gettext("outliers"), end = ".* \n", file = file, sep = "")
               }
               
               ind.sup = c(param$ind.sup, extrem)
               row.w = param$row.w[- which(ind.names %in% ind.names[extrem])]
               res.out = PCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, graph = FALSE, scale.unit = param$scale, row.w = row.w, col.w = param$col.w, ncp = param$ncp.mod)
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res.out, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 2, sep = "."), " - ", gettext("Individuals factor map (PCA) after correction"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of"), " ", length(extrem), " ", gettext("outliers"), end = ".* \n", file = file, sep = "")
               }
               
               selec.res = selection(res, dim = 1:2, margin = 2, selec = Vselec, coef = Vcoef)
               drawn = selec.res[[1]]
               what.drawn = selec.res[[2]]
               
               if(graph) {plot.PCA(memory, choix = 'var', select = drawn, title = gettext("Variables factor map (PCA) before correction"), cex = cex)}
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(memory, choix = 'var', select = drawn, title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
               writeRmd("**", paste(figure.title, 3, sep = "."), " - ", gettext("Variables factor map (PCA) before correction"), "**", file = file, sep = "")
               if(!is.null(param$quali.sup)) {
                 writeRmd("*", gettext("The variables in black are considered as active whereas those in blue are illustrative"), ".*", file = file, sep = "")
               }
               writeRmd(what.drawn, file = file, sep = "")
               
               writeRmd("\n- - -", file = file, end = "\n\n")
               
               if(is.null(c(param$quali.sup, param$quanti.sup))){
                 if(is.null(param$ind.sup)){
                   centred = scale(data)
                 } else {
                   centred = scale(data[-param$ind.sup,])
                 }
               } else {
                 if(is.null(param$ind.sup)){
                   centred = scale(data[, -c(param$quali.sup, param$quanti.sup)])
                 } else {
                   centred = scale(data[-param$ind.sup, -c(param$quali.sup, param$quanti.sup)])
                 }
               }
               
               for(i in extrem){
                 writeRmd("**", gettext("The individual"), " ", i, end = "** :\n\n", file = file, sep = "")
                 
                 writeRmd("-", gettext("takes very high values for the variable(s)"), end = " :\n", file = file)
                 liste.pos = names(which(sort(centred[rownames(data)[i],], decreasing = TRUE) > 2))
                 if(length(liste.pos) == 1) {
                   writeRmd("*", liste.pos, end = "*.\n\n", file = file, sep = "")
                 } else if(length(liste.pos) > 1 & length(liste.pos) < nmax) {
                   writeRmd("*", paste(paste(liste.pos[- length(liste.pos)], collapse = "*, *"), liste.pos[length(liste.pos)], sep = gettext("* and *")), "* (", gettext("variables are sorted from the strongest"), end = ").\n\n", file = file, sep = "")
                 } else if(length(liste.pos) >= nmax) {
                   writeRmd("*", paste(liste.pos[1:(nmax - 1)], collapse = "*, *"), gettext("* and *"), liste.pos[nmax], "* (", gettext("variables are sorted from the strongest"), end = ").\n\n", file = file, sep = "")
                 }
                 
                 writeRmd("-", gettext("takes very low values for the variable(s)"), end = " :\n", file = file)
                 liste.neg = names(which(sort(centred[rownames(data)[i],], decreasing = FALSE) < -2))
                 if(length(liste.neg) == 1) {
                   writeRmd("*", liste.neg, end = "*.\n\n", file = file, sep = "")
                 } else if(length(liste.neg) > 1 & length(liste.neg) < nmax) {
                   writeRmd("*", paste(paste(liste.neg[- length(liste.neg)], collapse = "*, *"), liste.neg[length(liste.neg)], sep = gettext("* and *")), "* (", gettext("variables are sorted from the strongest"), end = ").\n\n", file = file, sep = "")
                 } else if(length(liste.neg) >= nmax) {
                   writeRmd("*", paste(liste.neg[1:(nmax - 1)], collapse = "*, *"), gettext("* and *"), liste.neg[nmax], "* (", gettext("variables are sorted from the strongest"), end = ").\n\n", file = file, sep = "")
                 }
               }
               
               if(length(extrem) == 1) {
                 writeRmd(gettext("This outlier is suppressed from the analysis and a second one is performed on the rest of the individuals"), end = ".\n", file = file)
               } else {
                 writeRmd(gettext("These outliers are suppressed from the analysis and a second one is performed on the rest of the individuals"), end = ".\n", file = file)
               }
               suppr = which(ind.names %in% ind.names[extrem])
               ind.sup = which(ind.names %in% ind.names[param$ind.sup]) %dim0% NULL
               res = PCA(param$data[-suppr, ], quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, 
                         graph = FALSE, scale.unit = param$scale, row.w = param$row.w[-suppr], col.w = param$col.w, ncp = param$ncp.mod)
             },
             
             MCA = {
               if(graph) {
                 plot.MCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = ind.names[extrem], title = gettext("Individuals factor map (MCA) before correction"), cex = cex)
                 plot.MCA(res, choix = 'ind', invisible = c('var', 'quali'), select = ind.names[extrem], title = gettext("Individuals factor map (MCA) after correction"), cex = cex)
               }
               
               drawn = ind.names[extrem]
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(memory, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 1, sep = "."), " - ", gettext("Individuals factor map (MCA) before correction"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of"), " ", length(extrem), " ", gettext("outliers"), end = ".* \n", file = file, sep = "")
               }
               
               ind.sup = c(param$ind.sup, extrem)
               row.w = param$row.w[- which(ind.names %in% ind.names[extrem])]
               res.out = MCA(param$data, quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, graph = FALSE, row.w = row.w, ncp = param$ncp.mod)
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res.out, choix = 'ind', invisible = c('var', 'quali'), select = drawn, title = '', cex = cex)", file = file, stop = TRUE, end = "\n\n")
               writeRmd("**", paste(figure.title, 2, sep = "."), " - ", gettext("Individuals factor map (MCA) after correction"), end = ".** \n", file = file, sep = "")
               if(length(extrem) == 1){
                 writeRmd("*", gettext("Highlighting of an outlier"), end = ".* \n", file = file, sep = "")
               } else {
                 writeRmd("*", gettext("Highlighting of"), " ", length(extrem), " ", gettext("outliers"), end = ".* \n", file = file, sep = "")
               }
               
               selec.res = selection(res, dim = 1:2, margin = 2, selec = Vselec, coef = Vcoef)
               drawn = selec.res[[1]]
               what.drawn = selec.res[[2]]
               
               if(graph) {plot.MCA(memory, choix = 'ind', invisible = 'ind', selectMod = drawn, title = gettext("Variables factor map (MCA) before correction"), cex = cex)}
               
               writeRmd(file = file)
               writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
               dump("drawn", file = file, append = TRUE)
               writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(memory, choix = 'ind', invisible = 'ind', selectMod = drawn, title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
               writeRmd("**", paste(figure.title, 3, sep = "."), " - ", gettext("Variables factor map (MCA) before correction"), "**", file = file, sep = "")
               if(!is.null(param$quali.sup)) {
                 writeRmd("*", gettext("The factors in red are considered as active whereas those in green are illustrative"), ".*", file = file, sep = "")
               }
               writeRmd(what.drawn, file = file, sep = "")
               
               writeRmd("\n- - -", file = file, end = "\n\n")
               
               if(is.null(c(param$quali.sup, param$quanti.sup))){
                 if(is.null(param$ind.sup)){
                   centred = scale(tab.disjonctif(data))
                 } else {
                   centred = scale(tab.disjonctif(data[-param$ind.sup,]))
                 }
               } else {
                 if(is.null(param$ind.sup)){
                   centred = scale(tab.disjonctif(data[, -c(param$quali.sup, param$quanti.sup)]))
                 } else {
                   centred = scale(tab.disjonctif(data[-param$ind.sup, -c(param$quali.sup, param$quanti.sup)]))
                 }
               }
               
               for(i in extrem){
                 writeRmd("**", gettext("The individual"), " ", i, end = "** :\n\n", file = file, sep = "")
                 
                 writeRmd("-", gettext("is characterized by the factor(s)"), end = " :\n", file = file)
                 liste.pos = names(which(sort(centred[rownames(data)[i],], decreasing = TRUE) > 2))
                 if(length(liste.pos) == 1) {
                   writeRmd("*", liste.pos, end = "*.\n\n", file = file, sep = "")
                 } else if(length(liste.pos) > 1 & length(liste.pos) < nmax) {
                   writeRmd("*", paste(liste.pos, collapse = "*, *"), "* (", gettext("factors are sorted from the strongest"), end = ").\n\n", file = file, sep = "")
                 } else if(length(liste.pos) >= nmax) {
                   writeRmd("*", paste(liste.pos[1:(nmax - 1)], collapse = "*, *"), gettext("* and *"), liste.pos[nmax], "* (", gettext("factors are sorted from the strongest"), end = ").\n\n", file = file, sep = "")
                 }
               }
               
               if(length(extrem) == 1) {
                 writeRmd(gettext("This outlier is suppressed from the analysis and a second one is performed on the rest of the individuals"), end = ".\n", file = file)
               } else {
                 writeRmd(gettext("These outliers are suppressed from the analysis and a second one is performed on the rest of the individuals"), end = ".\n", file = file)
               }
               suppr = which(ind.names %in% ind.names[extrem])
               ind.sup = which(ind.names %in% ind.names[param$ind.sup]) %dim0% NULL
               res = MCA(param$data[-suppr, ], quanti.sup = param$quanti.sup, quali.sup = param$quali.sup, ind.sup = ind.sup, 
                         graph = FALSE, row.w = param$row.w[-suppr], ncp = param$ncp.mod)
             },
             
             MFA = {},
             
             HMFA = {},
             
             DMFA = {},
             
             FAMD = {},
             
             GPA = {},
             
             HCPC = {})
      
    } else {
      writeRmd(gettext("The analysis of the graphs does not detect any outlier"), file = file, end = ".\n")
      res.out = res
      memory = NULL
    }
    
    writeRmd("\n- - -", file = file, end = "\n\n")
    
    list(new.res = res, res.out = res.out, memory = memory, N = length(extrem), ID = extrem)
  }
