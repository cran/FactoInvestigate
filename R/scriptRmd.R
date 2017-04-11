scriptRmd <-
function(file, output = "code.R") {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to read"))}
    
    texte = readLines(file, n = -1) # On lit le fichier Rmd
    script = c()
    script.idx = grep("```", texte)
    for(i in 1:(length(script.idx)/2)) {
      for(j in (script.idx[2 * i - 1] + 1):(script.idx[2 * i] - 1))
        script[length(script) + 1] = texte[j]
    }
    
    script = script[-grep("mar = c", script)]
    script = script[-grep("'Workspace.RData'", script)]
    
    cat(script, file = output, sep = "\n", append = FALSE)
  }
