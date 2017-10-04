readRmd <- function(file, document = "html_document") {

# Function 1
open_file <- function(path=NULL){

    OS = .Platform$OS.type

    if(OS == "unix"){
      if(Sys.info()["sysname"] == "Linux") OS = "linux"
      else OS = "mac"
    }

    switch(OS,
           windows = shell.exec(path),
           mac = system(paste0("open ", path)),
           linux = {
             if(interactive()) system(paste0("xdg-open ", path))
             else cat("File path: ", path, "\n")
             }
           )
}

# Function 2
# preview_site <- function (path) {
  # utils::browseURL(path)
# }


    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to read"))}
    
    if(any(!document %in% c("word_document", "pdf_document", "html_document"))) 
    {return(warning("the parameter 'document' should only take 'word_document', 'pdf_document' or 'html_document' as value"))}
    document = unique(document)
    
    exe = render(file, document, quiet = TRUE)
    sapply(exe, open_file)
  }
