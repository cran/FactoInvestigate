eigenRef <-
function(res, dim = NULL, q = 0.95, time = "10000L", parallel = TRUE) {
    if(any(dim < 0)) {return(warning("the parameter 'dim' should only take positive values"))}
    
    dim = unique(dim)
    if(!is.numeric(dim) & !is.null(dim)) {return(warning("the argument 'dim' has to be a numeric vector"))}
    if(any(dim < 1)) {return(warning("the 'dim' vector elements must all be positive"))}
    
    if(!is.character(time)) {return(warning("the argument 'time' has to be a character chain"))}
    if(length(grep("[sL]", time)) == 0) {return(warning("the argument 'time' must specifie the desired unity : add 's' for second or 'L' for the number of repetitions"))}
    
    if(length(grep("L", time)) != 0) {
      time = as.numeric(strsplit(time, "L")[[1]])
      rep = TRUE
    } else {
      time = as.numeric(strsplit(time, "s")[[1]])
      rep = FALSE
    }
    
    if(!is.numeric(q)) {return(warning("the argument 'q' has to be a numeric vector"))}
    if(q < 0 | q > 1) {return(warning("the 'q' value has to be set between 0 and 1"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    if(parallel) {
      nb.cores = detectCores() - 1
      if(nb.cores == 0) {nb.cores = 1}
      clust <- makeCluster(nb.cores)
      clusterEvalQ(clust, require(FactoMineR))
    }
    
    switch(analyse,
           PCA = {
             ind = param$ind
             var = param$var
             row.w = param$row.w
             col.w = param$col.w
             scale = param$scale
             
             if(parallel) {
               clusterExport(clust, c("ind", "var", "row.w", "col.w", "scale"), envir = environment())
               
               if(rep) {
                 Q = t(parSapply(clust, 1:time, function(x, ind, var, row.w, col.w, scale) {
                   X = matrix(rnorm(ind * var), ind, var)
                   mp = crossprod(row.w / sum(row.w), as.matrix(X))
                   sp = sqrt(crossprod(row.w / sum(row.w), as.matrix(X) ^ 2))
                   sp[sp <= 1e-16] <- 1
                   X <- t(t(X) - as.vector(mp))
                   if(scale) {
                     X <- t(t(X) / as.vector(sp))
                   }
                   svd.triplet(X, row.w = row.w, col.w = col.w)$vs^2
                 }, ind = ind, var = var, row.w = row.w, col.w = col.w, scale = scale))
                 
                 n = time
               } else {
                 t1 = Sys.time()
                 Q = t(parSapply(clust, 1:100, function(x, ind, var, row.w, col.w, scale) {
                   X = matrix(rnorm(ind * var), ind, var)
                   mp = crossprod(row.w / sum(row.w), as.matrix(X))
                   sp = sqrt(crossprod(row.w / sum(row.w), as.matrix(X) ^ 2))
                   sp[sp <= 1e-16] <- 1
                   X <- t(t(X) - as.vector(mp))
                   if(scale) {
                     X <- t(t(X) / as.vector(sp))
                   }
                   svd.triplet(X, row.w = row.w, col.w = col.w)$vs^2
                 }, ind = ind, var = var, row.w = row.w, col.w = col.w, scale = scale))
                 duree = as.numeric(difftime(Sys.time(), t1)) / 100
                 
                 n = round(time / duree) - 100
                 if(n < 1) {n = 1}
                 
                 Q = rbind(Q, t(parSapply(clust, 1:n, function(x, ind, var, row.w, col.w, scale) {
                   X = matrix(rnorm(ind * var), ind, var)
                   mp = crossprod(row.w / sum(row.w), as.matrix(X))
                   sp = sqrt(crossprod(row.w / sum(row.w), as.matrix(X) ^ 2))
                   sp[sp <= 1e-16] <- 1
                   X <- t(t(X) - as.vector(mp))
                   if(scale) {
                     X <- t(t(X) / as.vector(sp))
                   }
                   svd.triplet(X, row.w = row.w, col.w = col.w)$vs^2
                 }, ind = ind, var = var, row.w = row.w, col.w = col.w, scale = scale)))
                 
                 n = n + 100
               }
               stopCluster(clust)
             } else {
               if(rep) {
                 Q = t(sapply(1:time, function(x, ind, var, row.w, col.w, scale) {
                   X = matrix(rnorm(ind * var), ind, var)
                   mp = crossprod(row.w / sum(row.w), as.matrix(X))
                   sp = sqrt(crossprod(row.w / sum(row.w), as.matrix(X) ^ 2))
                   sp[sp <= 1e-16] <- 1
                   X <- t(t(X) - as.vector(mp))
                   if(scale) {
                     X <- t(t(X) / as.vector(sp))
                   }
                   svd.triplet(X, row.w = row.w, col.w = col.w)$vs^2
                 }, ind = ind, var = var, row.w = row.w, col.w = col.w, scale = scale))
                 
                 n = time
               } else {
                 t1 = Sys.time()
                 Q = t(sapply(1:100, function(x, ind, var, row.w, col.w, scale) {
                   X = matrix(rnorm(ind * var), ind, var)
                   mp = crossprod(row.w / sum(row.w), as.matrix(X))
                   sp = sqrt(crossprod(row.w / sum(row.w), as.matrix(X) ^ 2))
                   sp[sp <= 1e-16] <- 1
                   X <- t(t(X) - as.vector(mp))
                   if(scale) {
                     X <- t(t(X) / as.vector(sp))
                   }
                   svd.triplet(X, row.w = row.w, col.w = col.w)$vs^2
                 }, ind = ind, var = var, row.w = row.w, col.w = col.w, scale = scale))
                 duree = as.numeric(difftime(Sys.time(), t1)) / 100
                 
                 n = round(time / duree) - 100
                 if(n < 1) {n = 1}
                 
                 Q = rbind(Q, t(sapply(1:n, function(x, ind, var, row.w, col.w, scale) {
                   X = matrix(rnorm(ind * var), ind, var)
                   mp = crossprod(row.w / sum(row.w), as.matrix(X))
                   sp = sqrt(crossprod(row.w / sum(row.w), as.matrix(X) ^ 2))
                   sp[sp <= 1e-16] <- 1
                   X <- t(t(X) - as.vector(mp))
                   if(scale) {
                     X <- t(t(X) / as.vector(sp))
                   }
                   svd.triplet(X, row.w = row.w, col.w = col.w)$vs^2
                 }, ind = ind, var = var, row.w = row.w, col.w = col.w, scale = scale)))
                 
                 n = n + 100
               }
             }
             
           },
           
           CA = {
             row = param$row
             col = param$col
             ind = param$ind
             row.w = param$row.w
             
             if(parallel) {
               clusterExport(clust, c("ind", "col", "row", "row.w"), envir = environment())
               
               if(rep) {
                 Q = t(parSapply(clust, 1:time, function(x, ind, col, row, row.w) {
                   X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   while(any(apply(X, 1, sum) == 0) | any(apply(X, 2, sum) == 0)) {
                     X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, col = col, row = row, row.w = row.w))
                 
                 n = time
               } else {
                 t1 = Sys.time()
                 Q = t(parSapply(clust, 1:100, function(x, ind, col, row, row.w) {
                   X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   while(any(apply(X, 1, sum) == 0) | any(apply(X, 2, sum) == 0)) {
                     X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, col = col, row = row, row.w = row.w))
                 duree = as.numeric(difftime(Sys.time(), t1)) / 100
                 
                 n = round(time / duree) - 100
                 if(n < 1) {n = 1}
                 
                 Q = rbind(Q, t(parSapply(clust, 1:n, function(x, ind, col, row, row.w) {
                   X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   while(any(apply(X, 1, sum) == 0) | any(apply(X, 2, sum) == 0)) {
                     X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, col = col, row = row, row.w = row.w)))
                 
                 n = n + 100
               }
               stopCluster(clust)
             } else {
               if(rep) {
                 Q = t(sapply(1:time, function(x, ind, col, row, row.w) {
                   X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   while(any(apply(X, 1, sum) == 0) | any(apply(X, 2, sum) == 0)) {
                     X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, col = col, row = row, row.w = row.w))
                 
                 n = time
               } else {
                 t1 = Sys.time()
                 Q = t(sapply(1:100, function(x, ind, col, row, row.w) {
                   X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   while(any(apply(X, 1, sum) == 0) | any(apply(X, 2, sum) == 0)) {
                     X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, col = col, row = row, row.w = row.w))
                 duree = as.numeric(difftime(Sys.time(), t1)) / 100
                 
                 n = round(time / duree) - 100
                 if(n < 1) {n = 1}
                 
                 Q = rbind(Q, t(sapply(1:n, function(x, ind, col, row, row.w) {
                   X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   while(any(apply(X, 1, sum) == 0) | any(apply(X, 2, sum) == 0)) {
                     X <- table(factor(sample(row, ind, replace = TRUE), 1:row), factor(sample(col, ind, replace = TRUE), 1:col))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, col = col, row = row, row.w = row.w)))
                 
                 n = n + 100
               }
             }
             
           },
           
           CaGalt = {},
           
           MCA = {
             data = param$data
             ind = param$ind
             row.w = param$row.w
             
             factors = sapply(rownames(res$var$eta2), function(x) {length(levels(data[, x]))})
             
             if(parallel) {
               clusterExport(clust, c("ind", "factors", "row.w"), envir = environment())
               
               if(rep) {
                 Q = t(parSapply(clust, 1:time, function(x, ind, factors, row.w) {
                   X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   while(ncol(X) != sum(factors)) {
                     X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, factors = factors, row.w = row.w))
                 
                 n = time
               } else {
                 t1 = Sys.time()
                 Q = t(parSapply(clust, 1:100, function(x, ind, factors, row.w) {
                   X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   while(ncol(X) != sum(factors)) {
                     X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, factors = factors, row.w = row.w))
                 duree = as.numeric(difftime(Sys.time(), t1)) / 100
                 
                 n = round(time / duree) - 100
                 if(n < 1) {n = 1}
                 
                 Q = rbind(Q, t(parSapply(clust, 1:n, function(x, ind, factors, row.w) {
                   X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   while(ncol(X) != sum(factors)) {
                     X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, factors = factors, row.w = row.w)))
                 
                 n = n + 100
               }
               stopCluster(clust)
             } else {
               if(rep) {
                 Q = t(sapply(1:time, function(x, ind, factors, row.w) {
                   X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   while(ncol(X) != sum(factors)) {
                     X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, factors = factors, row.w = row.w))
                 
                 n = time
               } else {
                 t1 = Sys.time()
                 Q = t(sapply(1:100, function(x, ind, factors, row.w) {
                   X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   while(ncol(X) != sum(factors)) {
                     X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, factors = factors, row.w = row.w))
                 duree = as.numeric(difftime(Sys.time(), t1)) / 100
                 
                 n = round(time / duree) - 100
                 if(n < 1) {n = 1}
                 
                 Q = rbind(Q, t(sapply(1:n, function(x, ind, factors, row.w) {
                   X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   while(ncol(X) != sum(factors)) {
                     X <- tab.disjonctif(sapply(factors, function(x, ind) {as.factor(sample(x, ind, replace = TRUE))}, ind = ind))
                   }
                   X <- X * (row.w / sum(X))
                   svd.triplet(t(t(X / rowSums(X)) / colSums(X)) - 1, row.w = rowSums(X), col.w = colSums(X))$vs^2
                 }, ind = ind, factors = factors, row.w = row.w)))
                 
                 n = n + 100
               }
             }
             
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
    if(is.null(dim)) {
      inertia = apply(apply(Q, 1, function(x) {cumsum(x) / sum(x)}), 1, function(y, q) quantile(y, q, names = FALSE), q = q)
      names(inertia) = paste("dim", 1:ncol(Q), sep = ".")
    } else {
      inertia = apply(apply(Q, 1, function(x) {cumsum(x) / sum(x)}), 1, function(y, q) quantile(y, q, names = FALSE), q = q)[dim]
      names(inertia) = paste("dim", dim, sep = ".")
    }
    
    list(datasets = n, 
         quantile = q,
         inertia = inertia)
  }
