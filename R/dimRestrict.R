dimRestrict <-
function(res, file = "", rand = NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    obs = res$eig[,2]
    
    if(is.null(rand)) {
      ref = eigenRef(res, time = "10s")
      rand = c(ref$inertia[1], diff(ref$inertia)) * 100
    }
    
    dim = min(length(obs), length(rand))
    
    ncp = which(obs[1:dim] < rand[1:dim])[[1]] - 1
    
    return(ncp)
  }
