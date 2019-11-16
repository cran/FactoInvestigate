graphInd <-
function(res, file = "", dim = 1:2, Iselec = "contrib", Icoef = 1, figure.title = "Figure", graph = TRUE, cex = 0.7, codeGraphInd = NULL, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(Iselec) & !is.character(Iselec)) {return(warning("the argument 'Iselec' should be a numeric or character vector"))}
    if(!is.numeric(Icoef)) {return(warning("the argument 'Icoef' must be numeric"))}
    if(Icoef < 0) {return(warning("the argument 'Icoef' must be positive"))}
    
    if(!is.numeric(cex)) {return(warning("the argument 'cex' must be numeric"))}
    if(cex < 0) {return(warning("the argument 'cex' must be positive"))}
    
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) | length(dim) != 2) {return(warning("the argument 'dim' has to be a 2 dimensionnal numeric vector"))}
    if(any(dim < 0)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    switch(analyse,
         PCA = {
           if (is.null(codeGraphInd)){
			 selec.res = selection(res, dim = dim, margin = 1, selec = Iselec, coef = Icoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph) plot.PCA(res, select = drawn, axes = c(dim[1],dim[2]), choix = 'ind', invisible = 'quali', title = gettext("Individuals factor map (PCA)",domain="R-FactoInvestigate"), cex = cex)
           } else{
			   eval(parse(text=codeGraphInd))
           }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             if (is.null(codeGraphInd)) writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res, select = drawn, axes = c(", dim[1], ",", dim[2],
                      "), choix = 'ind', invisible = 'quali', title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             else writeRmd(paste0("par(mar = c(4.1, 4.1, 1.1, 2.1))\n",codeGraphInd), stop = TRUE, sep = "", file = file, end = "\n\n")
			 
             writeRmd("**", figure.title, " - ", gettext("Individuals factor map (PCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             if (is.null(codeGraphInd)) {
			   if(!is.null(param$ind.sup))  writeRmd("*", gettext("The individuals in black are considered as active whereas the individuals in blue are illustrative",domain="R-FactoInvestigate"), ".*", file = file, sep = "")
               writeRmd(what.drawn, file = file, sep = "")
			 }
           },
           
           MCA = {
           if (is.null(codeGraphInd)){
             selec.res = selection(res, dim = dim, margin = 1, selec = Iselec, coef = Icoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph)  plot.MCA(res, select = drawn, axes = c(dim[1],dim[2]), choix = 'ind', invisible = c('var', 'quali'), title = gettext("Individuals factor map (MCA)",domain="R-FactoInvestigate"), cex = cex)
           } else{
			   eval(parse(text=codeGraphInd))
           }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             if (is.null(codeGraphInd)) writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res, select = drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'ind', invisible = c('var', 'quali'), title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             else writeRmd(paste0("par(mar = c(4.1, 4.1, 1.1, 2.1))\n",codeGraphInd), stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Individuals factor map (MCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             if (is.null(codeGraphInd)) {
               if(!is.null(param$ind.sup)) writeRmd("*", gettext("The individuals in light blue are considered as active whereas the individuals in dark blue are illustrative",domain="R-FactoInvestigate"), ".*", file = file, sep = "")
               writeRmd(what.drawn, file = file, sep = "")
             }
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
  }
