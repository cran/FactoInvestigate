getParam <-
function(res) {
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC", "HCPCshiny"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD' or 'HCPC'"))}
    
    switch(analyse,
           PCA = {
             data = res$call$X
             ind = nrow(res$ind$coord)
             var = nrow(res$var$coord)
             
             quanti.sup = which(names(data) %in% names(res$call$quanti.sup)) %dim0% NULL
             quali.sup = res$call$quali.sup$numero
             ind.sup = res$call$ind.sup
             row.w = res$call$row.w
             row.w.init = res$call$row.w.init
             col.w = res$call$col.w
             scale = res$call$scale.unit
             ncp.mod = res$call$ncp
             
             
             if(!is.null(quali.sup)) {
               modalites = lapply(data[, quali.sup, drop = FALSE], levels)
             } else {
               modalites = NULL
             }
             
             list(data = data, ind = ind, var = var, quanti.sup = quanti.sup, quali.sup = quali.sup, ind.sup = ind.sup, 
                  row.w = row.w, row.w.init = row.w.init, col.w = col.w, scale = scale, ncp.mod = ncp.mod, modalites = modalites)
           },
           
           CA = {
             data = res$call$Xtot
             ind = res$call$N
             row = nrow(res$row$coord)
             col = nrow(res$col$coord)
             
             row.sup = which(rownames(data) %in% rownames(res$row.sup$coord)) %dim0% NULL
             col.sup = which(names(data) %in% rownames(res$col.sup$coord)) %dim0% NULL
             quanti.sup = res$call$quanti.sup
             quali.sup = res$call$quali.sup
             row.w = res$call$row.w
             ncp.mod = res$call$ncp
             
             if(!is.null(quali.sup)) {
               modalites = lapply(data[, quali.sup, drop = FALSE], levels)
             } else {
               modalites = NULL
             }
             
             list(data = data, ind = ind, row = row, col = col, row.sup = row.sup, col.sup = col.sup, 
                  quanti.sup = quanti.sup, quali.sup = quali.sup, row.w = row.w, ncp.mod = ncp.mod, modalites = modalites)
           },
           
           CaGalt = {},
           
           MCA = {
             data = res$call$X
             ind = nrow(res$ind$coord)
             var = nrow(res$var$coord)
             
             quanti.sup = res$call$quanti.sup
             quali.sup = res$call$quali.sup
             ind.sup = res$call$ind.sup
             row.w = res$call$row.w
             ncp.mod = res$call$ncp
             
             if(!is.null(quanti.sup)) {
               modalites = lapply(data[, - quanti.sup, drop = FALSE], levels)
             } else {
               modalites = lapply(data, levels)
             }
             
             list(data = data, ind = ind, var = var, quanti.sup = quanti.sup, quali.sup = quali.sup, ind.sup = ind.sup, 
                  row.w = row.w, ncp.mod = ncp.mod, modalites = modalites)
           },
           
           MFA = {
             data = res$call$X
             group = res$call$group 
			 type <- res$call$type
			 ind.sup <- res$call$ind.sup
			 row.w <- res$call$row.w 
			 num.group.sup <- res$call$num.group.sup
             ncp.mod <- res$call$ncp
             
             list(data = data, group=group, type = type, ind.sup = ind.sup, 
                  row.w = row.w, num.group.sup = num.group.sup, ncp.mod = ncp.mod)
		   },
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {
             list(data=res$call$t$res$X, nb.clust = res$call$t$nb.clust)
           })
    
  }
