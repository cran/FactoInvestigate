classif <-
function(res, file = "", dim = 1:2, nclust = -1, selec = "contrib", coef = 1, mmax = 10, nmax = 10, figure.title = "Figure", graph = TRUE, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(selec) & !is.character(selec)) {return(warning("the argument 'selec' should be a numeric or character vector"))}
    if(!is.numeric(coef)) {return(warning("the argument 'Icoef' must be numeric"))}
    if(coef < 0) {return(warning("the argument 'Mcoef' must be positive"))}
    
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) | length(dim) != 2) {return(warning("the argument 'dim' has to be a 2 dimensionnal numeric vector"))}
    if(any(dim < 0)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    selec.res = selection(res, dim = dim, margin = 1, selec = selec, coef = coef)
    drawn = selec.res[[1]]
    what.drawn = selec.res[[2]]
    
    switch(analyse,
           PCA = {
             res.hcpc = HCPC(res, nb.clust = nclust, graph = FALSE)
             writeRmd("res.hcpc = HCPC(res, nb.clust = ", nclust, ", graph = FALSE)", file = file, start = TRUE, stop = TRUE, options = "r, echo = FALSE", sep = "")
             
             CD = res.hcpc$desc.var #catdes(res.hcpc$data.clust, which(names(res.hcpc$data.clust) == "clust"))
             
             if(graph) {
               plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
             }
             
             writeRmd(file = file)
             writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')", file = file, stop = TRUE, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Ascending Hierarchical Classification of the individuals"), ".**", file = file, sep = "")
             writeRmd("*", gettext("The classification made on individuals reveals"), " ", length(levels(res.hcpc$data.clust$clust)), " ",gettext("clusters"),".*", file = file, sep = "", end = "\n\n")
             
             for(i in levels(res.hcpc$data.clust$clust)) { # pour chaque cluster dans le groupe
               writeRmd(file = file)
               if(length(rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]) != 0) {
                 if(length(rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]) == 1) {
                   writeRmd(gettext("The **cluster"), " ", i, "** ", gettext("is made of individuals such as"), " *", rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn], 
                            "*. ", gettext("This group is characterized by"), file = file, sep = "", end = " :\n\n")
                 } else {
                   ind.clust = rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]
                   if(length(ind.clust) > mmax) {
                     ind.clust = paste(paste(ind.clust[1:(mmax - 1)], collapse = "*, *"), ind.clust[mmax], sep = gettext("* and *"))
                   } else {
                     ind.clust = paste(paste(ind.clust[- length(ind.clust)], collapse = "*, *"), ind.clust[length(ind.clust)], sep = gettext("* and *"))
                   }
                   writeRmd(gettext("The **cluster"), " ", i, "** ", gettext("is made of individuals such as"), " *", ind.clust, "*. ", gettext("This group is characterized by"), file = file, sep = "", end = " :\n\n")
                 }
               } else {
                 writeRmd(gettext("The **cluster"), " ", i, "** ", gettext("is made of individuals sharing"), end = " :\n\n", file = file, sep = "")
               }
               if(!is.null(CD$quanti[[i]])) {
                 if(any(CD$quanti[[i]][, 1] > 0)) {
                   if(nrow(data.frame(CD$quanti[[i]])[CD$quanti[[i]][, 1] > 0,]) == 1) {
                     writeRmd("- ", gettext("high values for the variable"), " *", rownames(data.frame(CD$quanti[[i]])[CD$quanti[[i]][, 1] > 0,]), end = "*.\n", file = file, sep = "")
                   } else {
                     high.var = rownames(CD$quanti[[i]][CD$quanti[[i]][, 1] > 0,])
                     if(length(high.var) > nmax) {
                       high.var = paste(paste(high.var[1:(nmax - 1)], collapse = "*, *"), high.var[nmax], sep = gettext("* and *"))
                       writeRmd("- ", gettext("high values for variables like"), " *", high.var, "* (", gettext("variables are sorted from the strongest"), end = ").\n", file = file, sep = "")
                     } else {
                       high.var = paste(paste(high.var[- length(high.var)], collapse = "*, *"), high.var[length(high.var)], sep = gettext("* and *"))
                       writeRmd("- ", gettext("high values for the variables"), " *", high.var, "* (", gettext("variables are sorted from the strongest"), end = ").\n", file = file, sep = "")
                     }
                   }
                 }
                 if(any(CD$quanti[[i]][, 1] < 0)) {
                   if(nrow(data.frame(CD$quanti[[i]])[CD$quanti[[i]][, 1] < 0,]) == 1) {
                     writeRmd("- ", gettext("low values for the variable"), " *", rownames(data.frame(CD$quanti[[i]])[CD$quanti[[i]][, 1] < 0,]), end = "*.\n", file = file, sep = "") 
                   } else {
                     low.var = rownames(CD$quanti[[i]][CD$quanti[[i]][, 1] < 0,])
                     low.var = low.var[length(low.var):1] 
                     if(length(low.var) > nmax) {
                       low.var = paste(paste(low.var[1:(nmax - 1)], collapse = "*, *"), low.var[nmax], sep = gettext("* and *"))
                       writeRmd("- ", gettext("low values for variables like"), " *", low.var, "* (", gettext("variables are sorted from the weakest"), end = ").\n", file = file, sep = "")
                     } else {
                       low.var = paste(paste(low.var[- length(low.var)], collapse = "*, *"), low.var[length(low.var)], sep = gettext("* and *"))
                       writeRmd("- ", gettext("low values for the variables"), " *", low.var, "* (", gettext("variables are sorted from the weakest"), end = ").\n", file = file, sep = "")
                     }
                   }
                 }
               } else {
                 writeRmd("- ", gettext("variables whose values do not differ significantly from the mean"), end = ".\n", file = file)
               }
             }
           },
           
           CA = {
             res.hcpc = HCPC(res, nb.clust = nclust, graph = FALSE)
             writeRmd("res.hcpc = HCPC(res, nb.clust = ", nclust, ", graph = FALSE)", file = file, start = TRUE, stop = TRUE, options = "r, echo = FALSE", sep = "")
             
             CD = res.hcpc$desc.var #catdes(res.hcpc$data.clust, which(names(res.hcpc$data.clust) == "clust"))
             
             if(graph) {
               plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
             }
             
             writeRmd(file = file)
             writeRmd(file = file, start = TRUE, options = "r, echo = FALSE, fig.align = 'center', fig.height = 3.5, fig.width = 5.5", end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')", file = file, stop = TRUE, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Ascending Hierachical Classification of the rows"), ".**", file = file, sep = "")
             writeRmd("*", gettext("The classification made on rows reveals"), " ", length(levels(res.hcpc$data.clust$clust)), " clusters.*", file = file, sep = "", end = "\n\n")
             
             for(i in levels(res.hcpc$data.clust$clust)) { # pour chaque cluster dans le groupe
               writeRmd(file = file)
               if(length(rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]) != 0) {
                 if(length(rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]) == 1) {
                   writeRmd(gettext("The cluster"), " ", i, " ",gettext("is made of rows such as"), " *", rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn], 
                            "*. ", gettext("This group is characterized by"), file = file, sep = "", end = " :\n\n")
                 } else {
                   ind.clust = rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]
                   if(length(ind.clust) > mmax) {
                     ind.clust = paste(paste(ind.clust[1:(mmax - 1)], collapse = "*, *"), ind.clust[mmax], sep = gettext("* and *"))
                   } else {
                     ind.clust = paste(paste(ind.clust[- length(ind.clust)], collapse = "*, *"), ind.clust[length(ind.clust)], sep = gettext("* and *"))
                   }
                   writeRmd(gettext("The cluster"), " ", i, " ", gettext("is made of rows such as"), " *", ind.clust, "*. ", gettext("This group is characterized by"), file = file, sep = "", end = " :\n\n")
                 }
               } else {
                 writeRmd(gettext("The cluster"), i, gettext("is made of rows sharing"), end = " :\n\n", file = file)
               }
               if(!is.null(CD[[i]])) {
                 if(any(CD[[i]][, "v.test"] > 0)) {
                   if(nrow(data.frame(CD[[i]])[CD[[i]][, "v.test"] > 0,]) == 1) {
                     writeRmd("- ", gettext("high frequency for the factor"), " *", rownames(data.frame(CD[[i]])[CD[[i]][, 1] > 0,]), end = "*.\n", file = file, sep = "")
                   } else {
                     high.var = rownames(CD[[i]][CD[[i]][, "v.test"] > 0,])
                     if(length(high.var) > nmax) {
                       high.var = paste(paste(high.var[1:(nmax - 1)], collapse = "*, *"), high.var[nmax], sep = gettext("* and *"))
                       writeRmd("- ", gettext("high frequency for factors like"), " *", high.var, "* (", gettext("factors are sorted from the most common"), end = ").\n", file = file, sep = "")
                     } else {
                       high.var = paste(paste(high.var[- length(high.var)], collapse = "*, *"), high.var[length(high.var)], sep = gettext("* and *"))
                       writeRmd("- ", gettext("high frequency for the factors"), " *", high.var, "* (", gettext("factors are sorted from the most common"), end = ").\n", file = file, sep = "")
                     }
                   }
                 }
                 if(any(CD[[i]][, "v.test"] < 0)) {
                   if(nrow(data.frame(CD[[i]])[CD[[i]][, "v.test"] < 0,]) == 1) {
                     writeRmd("- ", gettext("low frequency for the factor"), " *", rownames(data.frame(CD[[i]])[CD[[i]][, 1] < 0,]), end = "*.\n", file = file, sep = "") 
                   } else {
                     low.var = rownames(CD[[i]][CD[[i]][, "v.test"] < 0,])
                     low.var = low.var[length(low.var):1] 
                     if(length(low.var) > nmax) {
                       low.var = paste(paste(low.var[1:(nmax - 1)], collapse = "*, *"), low.var[nmax], sep = gettext("* and *"))
                       writeRmd("- ", gettext("low frequency for factors like"), " *", low.var, "* (", gettext("factors are sorted from the rarest"), end = ").\n", file = file, sep = "")
                     } else {
                       low.var = paste(paste(low.var[- length(low.var)], collapse = "*, *"), low.var[length(low.var)], sep = gettext("* and *"))
                       writeRmd("- ", gettext("low frequency for the factors"), " *", low.var, "* (", gettext("factors are sorted from the rarest"), end = ").\n", file = file, sep = "")
                     }
                   }
                 }
               } else {
                 writeRmd("- ", gettext("factors whose frequency does not differ significantly from the mean"), end = ".\n", file = file)
               }
             }
           },
           
           CaGalt = {},
           
           MCA = {
             res.hcpc = HCPC(res, nb.clust = nclust, graph = FALSE)
             writeRmd("res.hcpc = HCPC(res, nb.clust = ", nclust, ", graph = FALSE)", file = file, start = TRUE, stop = TRUE, options = "r, echo = FALSE", sep = "")
             
             CD = res.hcpc$desc.var #catdes(res.hcpc$data.clust, which(names(res.hcpc$data.clust) == "clust"))
             
             if(graph) {
               plot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')
             }
             
             writeRmd(file = file)
             writeRmd(file = file, start = TRUE, options = options, end = "")
             dump("drawn", file = file, append = TRUE)
             writeRmd("par(mar = c(4.1, 4.1, 1.1, 2.1))\nplot.HCPC(res.hcpc, choice = 'map', draw.tree = FALSE, select = drawn, title = '')", file = file, stop = TRUE, end = "\n\n")
             
             writeRmd("**", figure.title, " - ", gettext("Ascending Hierarchical Classification of the individuals"), ".**", file = file, sep = "")
             writeRmd("*", gettext("The classification made on individuals reveals"), " ", length(levels(res.hcpc$data.clust$clust)), " clusters.*", file = file, sep = "", end = "\n\n")
             
             for(i in levels(res.hcpc$data.clust$clust)) { # pour chaque cluster dans le groupe
               writeRmd(file = file)
               if(length(rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]) != 0) {
                 if(length(rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]) == 1) {
                   writeRmd(gettext("The 1st cluster is made of individuals such as"), " *", rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]
                            , "*. ", gettext("This group is characterized by"), file = file, sep = "", end = " :\n\n")
                 } else {
                   ind.clust = rownames(res.hcpc$data.clust)[res.hcpc$data.clust$clust == i & rownames(res.hcpc$data.clust) %in% drawn]
                   if(length(ind.clust) > mmax) {
                     ind.clust = paste(paste(ind.clust[1:(mmax - 1)], collapse = "*, *"), ind.clust[mmax], sep = gettext("* and *"))
                   } else {
                     ind.clust = paste(paste(ind.clust[- length(ind.clust)], collapse = "*, *"), ind.clust[length(ind.clust)], sep = gettext("* and *"))
                   }
                   writeRmd(gettext("The cluster"), " ", i, " ", gettext("is made of individuals such as"), "*. ", gettext("This group is characterized by"), ind.clust, "*.", file = file, sep = "", end = " :\n\n")
                 }
               } else {
                 writeRmd(gettext("The cluster"), i, gettext("is made of individuals sharing"), end = " :\n\n", file = file)
               }
               if(!is.null(CD$category[[i]])) {
                 if(any(CD$category[[i]][, "v.test"] > 0)) {
                   if(nrow(data.frame(CD$category[[i]])[CD$category[[i]][, "v.test"] > 0,]) == 1) {
                     writeRmd("- ", gettext("high frequency for the factor"), " *", rownames(data.frame(CD$category[[i]])[CD$category[[i]][, 1] > 0,]), end = "*.\n", file = file, sep = "")
                   } else {
                     high.var = rownames(CD$category[[i]][CD$category[[i]][, "v.test"] > 0,])
                     if(length(high.var) > nmax) {
                       high.var = paste(paste(high.var[1:(nmax - 1)], collapse = "*, *"), high.var[nmax], sep = gettext("* and *"))
                       writeRmd("- ", gettext("high frequency for factors like"), " *", high.var, "* (", gettext("factors are sorted from the most common"), end = ").\n", file = file, sep = "")
                     } else {
                       high.var = paste(paste(high.var[- length(high.var)], collapse = "*, *"), high.var[length(high.var)], sep = gettext("* and *"))
                       writeRmd("- ", gettext("high frequency for the factors"), " *", high.var, "* (", gettext("factors are sorted from the most common"), end = ").\n", file = file, sep = "")
                     }
                   }
                 }
                 if(any(CD$category[[i]][, "v.test"] < 0)) {
                   if(nrow(data.frame(CD$category[[i]])[CD$category[[i]][, "v.test"] < 0,]) == 1) {
                     writeRmd("- ", gettext("low frequency for the factor"), " *", rownames(data.frame(CD$category[[i]])[CD$category[[i]][, 1] < 0,]), end = "*.\n", file = file, sep = "") 
                   } else {
                     low.var = rownames(CD$category[[i]][CD$category[[i]][, "v.test"] < 0,])
                     low.var = low.var[length(low.var):1] 
                     if(length(low.var) > nmax) {
                       low.var = paste(paste(low.var[1:(nmax - 1)], collapse = "*, *"), low.var[nmax], sep = gettext("* and *"))
                       writeRmd("- ", gettext("low frequency for factors like"), " *", low.var, "* (", gettext("factors are sorted from the rarest"), end = ").\n", file = file, sep = "")
                     } else {
                       low.var = paste(paste(low.var[- length(low.var)], collapse = "*, *"), low.var[length(low.var)], sep = gettext("* and *"))
                       writeRmd("- ", gettext("low frequency for the factors"), " *", low.var, "* (", gettext("factors are sorted from the rarest"), end = ").\n", file = file, sep = "")
                     }
                   }
                 }
               } else {
                 writeRmd("- ", gettext("factors whose frequency does not differ significantly from the mean"), end = ".\n", file = file)
               }
             }
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
    res.hcpc
  }
