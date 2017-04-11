writeRmd <-
function(..., file = "", append = TRUE, sep = " ", end = "\n", dump = FALSE, start = FALSE, stop = FALSE, options = NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.logical(append)) {return(warning("the argument 'append' must be logical"))}
    if(!is.logical(dump)) {return(warning("the argument 'dump' must be logical"))}
    if(!is.logical(start)) {return(warning("the argument 'start' must be logical"))}
    if(!is.logical(stop)) {return(warning("the argument 'stop' must be logical"))}
    
    if(!is.character(options) & !is.null(options)) {return(warning("the argument 'options' must be a character chain"))}
    
    cat(file = file, append = append)
    if(start) {
      cat("```", file = file, append = TRUE)
      if(!is.null(options)) {cat("{", options, "}", sep = "", file = file, append = TRUE)}
      cat("\n", file = file, append = TRUE)
    }
    
    if(dump) {
      dump(..., file = file, append = TRUE)
    } else {
      cat(..., file = file, append = TRUE, sep = sep)
    }
    if(stop) {cat("\n```", file = file, append = TRUE)}
    
    cat(end, file = file, append = TRUE)
  }
