selection <-
function(res, dim = 1:2, margin = 1, selec = "cos2", coef = 1) {
    if(!margin %in% 1:3) 
    {return(warning("the argument 'margin' should be an integer of value equal to 1, 2 or 3 (for individuals, active variables or illustrative variables selection)"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) | length(dim) != 2) {return(warning("the argument 'dim' has to be a 2 dimensionnal numeric vector"))}
    if(any(dim < 0)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    switch(analyse,
           PCA = {
             drawn = NULL
             what.drawn = NULL
             if(is.numeric(selec)) {
               drawn = selec
               if(margin == 1) {
                 what.drawn = paste("*", gettext("The labeled individuals are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               } else if(margin == 2) {
                 what.drawn = paste("*", gettext("The labeled variables are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               } else {
                 what.drawn = paste("*", gettext("The labeled factors are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               }
             } else {
               if(is.character(selec) & length(grep("contrib", selec)) != 0) {
                 k = as.numeric(unlist(strsplit(selec, " "))[2])
                 if(margin == 1) {
                   contrib = sort(apply(res$ind$contrib[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 } else if(margin == 2) {
                   contrib = sort(apply(res$var$contrib[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 } else {
                   return(warning("Illustrative variables does not contribute to the plane construction, thus they cannot be selected after that criterion"))
                 }
                 
                 if(is.na(k)) {
                   t = exp(-length(contrib) / 55) * length(contrib) * coef
                   
                   if(round(t, 0) == 0){
                     drawn = integer(0)
                   } else {
                     if(margin == 1) {
                       k = min(round(t, 0), max(round(t / 2, 0), length(contrib[contrib > mean(res$ind$contrib) ^ 2])))
                     } else {
                       k = min(round(t, 0), max(round(t / 2, 0), length(contrib[contrib > mean(res$var$contrib) ^ 2])))
                     }
                     drawn = names(contrib)[1:k]
                   }
                 } else {
                   if(k > 0 & k < 1) 
                   {drawn = names(contrib)[1:(length(contrib) * (1 - k))]} 
                   else 
                   {drawn = names(contrib)[1:k]}
                 }
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled individuals are those with the higher contribution to the plane construction",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled variables are those with the higher contribution to the plane construction",domain="R-FactoInvestigate"), ".*", sep = "")
                 }
               } else if(is.character(selec) & length(grep("cos2", selec)) != 0) {
                 k = as.numeric(unlist(strsplit(selec, " "))[2])
                 if(margin == 1) {
                   cos2 = sort(apply(res$ind$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 } else if(margin == 2) {
                   cos2 = sort(apply(res$var$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 } else {
                   cos2 = sort(apply(res$quali.sup$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 }
                 
                 if(margin == 1) {
                   if(!is.null(param$ind.sup)) {
                     cos2.sup = apply(res$ind.sup$cos2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                 } else if(margin == 2) {
                   if(!is.null(param$quanti.sup)) {
                     cos2.sup = apply(res$quanti.sup$cos2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                 }
                 
                 if(is.na(k)) {
                   t = exp(-length(cos2) / 55) * length(cos2) * coef
                   
                   if(round(t, 0) == 0){
                     drawn = integer(0)
                   } else {
                     k = min(round(t, 0), max(round(t / 2, 0), length(cos2[cos2 > 0.5])))
                     drawn = names(cos2)[1:k]
                   }
                 } else {
                   if(k > 0 & k < 1) 
                   {drawn = names(cos2)[cos2 >= k]} 
                   else 
                   {drawn = names(cos2)[1:k]}
                 }
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled individuals are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else if(margin == 2) {
                   what.drawn = paste("*", gettext("The labeled variables are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled factors are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 }
               } else if(is.character(selec)) {
                 drawn = selec
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled individuals are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 } else if(margin == 2) {
                   what.drawn = paste("*", gettext("The labeled variables are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled factors are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 }
               }
             }
           },
           
           CA = {
             drawn = NULL
             what.drawn = NULL
             if(is.numeric(selec)) {
               drawn = selec
               if(margin == 1) {
                 what.drawn = paste("*", gettext("The labeled rows are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               } else if(margin == 2) {
                 what.drawn = paste("*", gettext("The labeled columns are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               } else {
                 what.drawn = paste("*", gettext("The labeled variables are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               }
             } else {
               if(is.character(selec) & length(grep("contrib", selec)) != 0) {
                 k = as.numeric(unlist(strsplit(selec, " "))[2])
                 if(margin == 1) {
                   contrib = sort(apply(res$row$contrib[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 } else if(margin == 2) {
                   contrib = sort(apply(res$col$contrib[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 } else {
                   return(warning("Illustrative variables does not contribute to the plane construction, thus they cannot be selected after that criterion"))
                 }
                 
                 if(is.na(k)) {
                   t = exp(-length(contrib) / 55) * length(contrib) * coef
                   
                   if(round(t, 0) == 0){
                     drawn = integer(0)
                   } else {
                     if(margin == 1) {
                       k = min(round(t, 0), max(round(t / 2, 0), length(contrib[contrib > mean(res$row$contrib) ^ 2])))
                     } else {
                       k = min(round(t, 0), max(round(t / 2, 0), length(contrib[contrib > mean(res$col$contrib) ^ 2])))
                     }
                     drawn = names(contrib)[1:k]
                   }
                 } else {
                   if(k > 0 & k < 1) 
                   {drawn = names(contrib)[1:(length(contrib) * (1 - k))]} 
                   else 
                   {drawn = names(contrib)[1:k]}
                 }
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled rows are those with the higher contribution to the plane construction",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled columns are those with the higher contribution to the plane construction",domain="R-FactoInvestigate"), ".*", sep = "")
                 }
               } else if(is.character(selec) & length(grep("cos2", selec)) != 0) {
                 k = as.numeric(unlist(strsplit(selec, " "))[2])
                 if(margin == 1) {
                   cos2 = sort(apply(res$row$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 } else if(margin == 2) {
                   cos2 = sort(apply(res$col$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 } else {
                   cos2 = sort(apply(res$quanti.sup$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 }
                 
                 if(margin == 1) {
                   if(!is.null(param$row.sup)) {
                     cos2.sup = apply(res$row.sup$cos2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                 } else if(margin == 2) {
                   if(!is.null(param$col.sup)) {
                     cos2.sup = apply(res$col.sup$cos2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                   if(!is.null(param$quali.sup)) {
                     cos2.sup = apply(res$quali.sup$eta2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                 }
                 
                 if(is.na(k)) {
                   t = exp(-length(cos2) / 55) * length(cos2) * coef
                   
                   if(round(t, 0) == 0){
                     drawn = integer(0)
                   } else {
                     k = min(round(t, 0), max(round(t / 2, 0), length(cos2[cos2 > 0.5])))
                     drawn = names(cos2)[1:k]
                   }
                 } else {
                   if(k > 0 & k < 1) 
                   {drawn = names(cos2)[cos2 >= k]} 
                   else 
                   {drawn = names(cos2)[1:k]}
                 }
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled rows are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else if(margin == 2) {
                   what.drawn = paste("*", gettext("The labeled columns are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled variables are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 }
               } else if(is.character(selec)) {
                 drawn = selec
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled rows are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 } else if(margin == 2) {
                   what.drawn = paste("*", gettext("The labeled columns are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled variables are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 }
               }
             }
           },
           
           CaGalt = {},
           
           MCA = {
             drawn = NULL
             what.drawn = NULL
             if(is.numeric(selec)) {
               drawn = selec
               if(margin == 1) {
                 what.drawn = paste("*", gettext("The labeled individuals are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               } else {
                 what.drawn = paste("*", gettext("The labeled variables are those numbered",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
               }
             } else {
               if(is.character(selec) & length(grep("contrib", selec)) != 0) {
                 k = as.numeric(unlist(strsplit(selec, " "))[2])
                 if(margin == 1) {
                   contrib = sort(apply(res$ind$contrib[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 } else if(margin == 2) {
                   contrib = sort(apply(res$var$contrib[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 } else {
                   return(warning("Illustrative variables does not contribute to the plane construction, thus they cannot be selected after that criterion"))
                 }
                 
                 if(is.na(k)) {
                   t = exp(-length(contrib) / 55) * length(contrib) * coef
                   
                   if(round(t, 0) == 0){
                     drawn = integer(0)
                   } else {
                     if(margin == 1) {
                       k = min(round(t, 0), max(round(t / 2, 0), length(contrib[contrib > mean(res$ind$contrib) ^ 2])))
                     } else {
                       k = min(round(t, 0), max(round(t / 2, 0), length(contrib[contrib > mean(res$var$contrib) ^ 2])))
                     }
                     drawn = names(contrib)[1:k]
                   }
                 } else {
                   if(k > 0 & k < 1) 
                   {drawn = names(contrib)[1:(length(contrib) * (1 - k))]} 
                   else 
                   {drawn = names(contrib)[1:k]}
                 }
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled individuals are those with the higher contribution to the plane construction",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled variables are those with the higher contribution to the plane construction",domain="R-FactoInvestigate"), ".*", sep = "")
                 }
               } else if(is.character(selec) & length(grep("cos2", selec)) != 0) {
                 k = as.numeric(unlist(strsplit(selec, " "))[2])
                 if(margin == 1) {
                   cos2 = sort(apply(res$ind$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 } else if(margin == 2) {
                   cos2 = sort(apply(res$var$cos2[, dim, drop = FALSE], 1, sum), decreasing = TRUE)
                 } else {
                   cos2 = sort(apply(res$quanti.sup$coord[, dim, drop = FALSE], 1, function(x) {x[1] ^ 2 + x[2] ^ 2}), decreasing = TRUE)
                 }
                 
                 if(margin == 1) {
                   if(!is.null(param$ind.sup)) {
                     cos2.sup = apply(res$ind.sup$cos2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                 } else if(margin == 2) {
                   if(!is.null(param$quali.sup)) {
                     cos2.sup = apply(res$quali.sup$cos2[, dim, drop = FALSE], 1, sum)
                     cos2 = sort(c(cos2, cos2.sup), decreasing = TRUE)
                   }
                 }
                 
                 
                 if(is.na(k)) {
                   t = exp(-length(cos2) / 55) * length(cos2) * coef
                   
                   if(round(t, 0) == 0){
                     drawn = integer(0)
                   } else {
                     k = min(round(t, 0), max(round(t / 2, 0), length(cos2[cos2 > 0.5])))
                     drawn = names(cos2)[1:k]
                   }
                 } else {
                   if(k > 0 & k < 1) 
                   {drawn = names(cos2)[cos2 >= k]} 
                   else 
                   {drawn = names(cos2)[1:k]}
                 }
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled individuals are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled variables are those the best shown on the plane",domain="R-FactoInvestigate"), ".*", sep = "")
                 } 
               } else if(is.character(selec)) {
                 drawn = selec
                 if(margin == 1) {
                   what.drawn = paste("*", gettext("The labeled individuals are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 } else {
                   what.drawn = paste("*", gettext("The labeled variables are those named",domain="R-FactoInvestigate"), " ", paste(selec, collapse = ", "), ".*", sep = "")
                 }
               }
             }
           },
           
           MFA = {
             drawn <- NULL
             what.drawn <- NULL
           },
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
    list(drawn = drawn, what.drawn = what.drawn) 
  }
