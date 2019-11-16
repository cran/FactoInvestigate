factoGraph <-
function(res, file = "", dim = 1:2, hab = NULL, ellipse = TRUE, Iselec = "contrib", Vselec = "cos2", Rselec = "cos2", Cselec = "cos2", Mselec = "cos2", 
Icoef = 1, Vcoef = 1, Rcoef = 1, Ccoef = 1, Mcoef = 1, figure.title = "Figure", graph = TRUE, cex = 0.7, codeGraphInd = NULL, codeGraphVar = NULL ,codeGraphCA = NULL,options = NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    i = 1
    switch(analyse,
           PCA = {
			   graphInd(res, file = file, dim = dim, Iselec = Iselec, Icoef = Icoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, codeGraphInd = codeGraphInd, options = options)
             if((hab != "none") %dim0% TRUE & !is.null(param$quali.sup)) {
               writeRmd(file = file)
               i = i + 1
               graphHab(res, file = file, dim = dim, hab = hab, ellipse = ellipse, Iselec = Iselec, Icoef = Icoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, options = options)
             }
             i = i + 1
               graphVar(res, file = file, dim = dim, Vselec = Vselec, Vcoef = Vcoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, codeGraphVar = codeGraphVar, options = options)
               if(!is.null(param$quali.sup)) {
                 i = i + 1
                 graphSup(res, file = file, dim = dim, Mselec = Mselec, Mcoef = Mcoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, options = options)
               }
           },
           
           CA = {
             graphCA(res, file = file, dim = dim, Rselec = Rselec, Cselec = Cselec, Rcoef = Rcoef, Ccoef = Ccoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, codeGraphCA = codeGraphCA, options = options)
             if((hab != "none") %dim0% TRUE & !is.null(param$quali.sup)) {
               writeRmd(file = file)
               i = i + 1
               graphHab(res, file = file, dim = dim, hab = hab, ellipse = ellipse, Rselec = Rselec, Cselec = Cselec, Rcoef = Rcoef, Ccoef = Ccoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, options = options)
             }
             if(!is.null(param$quanti.sup)) {
               i = i + 1
               graphSup(res, file = file, dim = dim, Mselec = Mselec, Mcoef = Mcoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, options = options)
             }
           },
           
           CaGalt = {},
           
           MCA = {
			 graphInd(res, file = file, dim = dim, Iselec = Iselec, Icoef = Icoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, codeGraphInd = codeGraphInd, options = options)
             if((hab != "none") %dim0% TRUE & !is.null(param$quali.sup)) {
               writeRmd(file = file)
               i = i + 1
               graphHab(res, file = file, dim = dim, hab = hab, ellipse = ellipse, Iselec = Iselec, Icoef = Icoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, options = options)
             }
             i = i + 1
             graphVar(res, file = file, dim = dim, Vselec = Vselec, Vcoef = Vcoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, codeGraphVar = codeGraphVar, options = options)
             if(!is.null(param$quanti.sup)) {
               i = i + 1
               graphSup(res, file = file, dim = dim, Mselec = Mselec, Mcoef = Mcoef, figure.title = paste(figure.title, i, sep = "."), graph = graph, cex = cex, options = options)
             }
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
  }
