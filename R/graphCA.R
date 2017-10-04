graphCA <-
function(res, file = "", dim = 1:2, Rselec = "cos2", Cselec = "cos2", Rcoef = 1, Ccoef = 1, figure.title = "Figure", graph = TRUE, cex = 0.7, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(Rselec) & !is.character(Rselec)) {return(warning("the argument 'Rselec' should be a numeric or character vector"))}
    if(!is.numeric(Cselec) & !is.character(Cselec)) {return(warning("the argument 'Cselec' should be a numeric or character vector"))}
    
    if(!is.numeric(Rcoef)) {return(warning("the argument 'Rcoef' must be numeric"))}
    if(!is.numeric(Ccoef)) {return(warning("the argument 'Ccoef' must be numeric"))}
    
    if(Rcoef < 0) {return(warning("the argument 'Rcoef' must be positive"))}
    if(Ccoef < 0) {return(warning("the argument 'Ccoef' must be positive"))}
    
    if(!is.numeric(cex)) {return(warning("the argument 'cex' must be numeric"))}
    if(cex < 0) {return(warning("the argument 'cex' must be positive"))}
    
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) | length(dim) != 2) {return(warning("the argument 'dim' has to be a 2 dimensional numeric vector"))}
    if(any(dim < 0)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    switch(analyse,
           CA = {
             selec.res = selection(res, dim = dim, margin = 1, selec = Rselec, coef = Rcoef)
             r.drawn = selec.res[[1]]
             r.what.drawn = selec.res[[2]]
             
             selec.res = selection(res, dim = dim, margin = 2, selec = Cselec, coef = Ccoef)
             c.drawn = selec.res[[1]]
             c.what.drawn = selec.res[[2]]
             
             
             if(graph) {
               plot.CA(res, selectRow = r.drawn, selectCol = c.drawn, axes = dim[1]:dim[2], choix = 'CA', invisible = c('var', 'quali'), title = gettext("Overlayed factor map (CA)"), cex = cex)
             }
             writeRmd(file = file)
             writeRmd(start = TRUE, options = options, file = file, end = "")
             dump("r.drawn", file = file, append = TRUE)
             dump("c.drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.CA(res, selectRow = r.drawn, selectCol = c.drawn, axes = ", dim[1], ":", dim[2],
                      ", choix = 'CA', invisible = c('var', 'quali'), title = '', cex = cex)", stop = TRUE, sep = "", file = file, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Overlayed factor map (CA)"), "**", file = file, sep = "")
             if(!is.null(param$row.sup)) {
               writeRmd("*", gettext("The rows in light blue are considered as active whereas those in dark blue are illustrative"), ".*", file = file, sep = "")
             }
             if(!is.null(param$col.sup)) {
               writeRmd("*", gettext("The columns in light red are considered as active whereas those in dark red are illustrative"), ".*", file = file, sep = "")
             }
             writeRmd(r.what.drawn, file = file, sep = "")
             writeRmd(c.what.drawn, file = file, sep = "")
           },
           
           CaGalt = {
             
           })
    
  }
