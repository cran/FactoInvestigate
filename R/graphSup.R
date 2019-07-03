graphSup <-
function(res, file = "", dim = 1:2, Mselec = "cos2", Mcoef = 1, figure.title = "Figure", graph = TRUE, cex = 0.7, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(Mselec) & !is.character(Mselec)) {return(warning("the argument 'Mselec' should be a numeric or character vector"))}
    if(!is.numeric(Mcoef)) {return(warning("the argument 'Mcoef' must be numeric"))}
    if(Mcoef < 0) {return(warning("the argument 'Mcoef' must be positive"))}
    
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
             selec.res = selection(res, dim = dim, margin = 3, selec = Mselec, coef = Mcoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph) {
               plot.PCA(res, select = drawn, axes = dim[1]:dim[2], choix = 'ind', invisible = c('ind', 'ind.sup'), title = gettext("Qualitative factor map (PCA)",domain="R-FactoInvestigate"), cex = cex)
             }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.PCA(res, select = drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'ind', invisible = c('ind', 'ind.sup'), title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Qualitative factor map (PCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             writeRmd(what.drawn, file = file, sep = "")
           },
           
           CA = {
             selec.res = selection(res, dim = dim, margin = 3, selec = Mselec, coef = Mcoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph) {
               plot.CA(res, selectCol = drawn, axes = dim[1]:dim[2], choix = 'quanti.sup', title = gettext("Quantitative factor map (CA)",domain="R-FactoInvestigate"), cex = cex)
             }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.CA(res, selectCol = drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'quanti.sup', title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Quantitative factor map (CA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             writeRmd(what.drawn, file = file, sep = "")
           },
           
           CaGalt = {},
           
           MCA = {
             selec.res = selection(res, dim = dim, margin = 3, selec = Mselec, coef = Mcoef)
             drawn = selec.res[[1]]
             what.drawn = selec.res[[2]]
             
             if(graph) {
               plot.MCA(res, select = drawn, axes = dim[1]:dim[2], choix = 'quanti.sup', title = gettext("Quantitative factor map (MCA)",domain="R-FactoInvestigate"), cex = cex)
             }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.MCA(res, select = drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'quanti.sup', title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Quantitative factor map (MCA)",domain="R-FactoInvestigate"), "**", file = file, sep = "")
             writeRmd(what.drawn, file = file, sep = "")
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
  }
