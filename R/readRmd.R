readRmd <-
function(file, document = "html_document") {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to read"))}
    
    if(any(!document %in% c("word_document", "pdf_document", "html_document"))) 
    {return(warning("the parameter 'document' should only take 'word_document', 'pdf_document' or 'html_document' as value"))}
    document = unique(document)
    
    exe = render(file, document, quiet = TRUE)
    sapply(exe, shell.exec)
  }
