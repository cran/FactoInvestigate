graphVar <-
function(res, file = "", dim = 1:2, Vselec = "cos2", Vcoef = 1, figure.title = "Figure", graph = TRUE, cex = 0.7, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(Vselec) & !is.character(Vselec)) {return(warning("the argument 'Vselec' should be a numeric or character vector"))}
    if(!is.numeric(Vcoef)) {return(warning("the argument 'Vcoef' must be numeric"))}
    if(Vcoef < 0) {return(warning("the argument 'Vcoef' must be positive"))}
    
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
             selec.res = selection(res, dim = dim, margin = 2, selec = Vselec, coef = Vcoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph) {
               plot.PCA(res, select = drawn, axes = dim[1]:dim[2], choix = 'var', title = gettext("Variables factor map (PCA)",domain="R-FactoInvestigate"), cex = cex)
             }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res, select = drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'var', title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Variables factor map (PCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             if(!is.null(param$quanti.sup)) {
               writeRmd("*", gettext("The variables in black are considered as active whereas those in blue are illustrative",domain="R-FactoInvestigate"), ".*", file = file, sep = "")
             }
             writeRmd(what.drawn, file = file, sep = "")
           },
           
           MCA = {
             selec.res = selection(res, dim = dim, margin = 2, selec = Vselec, coef = Vcoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph) {
               plot.MCA(res, selectMod = drawn, axes = dim[1]:dim[2], choix = 'ind', invisible = 'ind', title = gettext("Variables factor map (MCA)",domain="R-FactoInvestigate"), cex = cex)
             }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res, selectMod = drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'ind', invisible = 'ind', title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Variables factor map (MCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             if(!is.null(param$quali.sup)) {
               writeRmd("*", gettext("The factors in red are considered as active whereas those in green are illustrative",domain="R-FactoInvestigate"), ".*", file = file, sep = "")
             }
             writeRmd(what.drawn, file = file, sep = "")
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
  }
