inertiaDistrib <-
function(res, file = "", ncp = NULL, q = 0.95, time = "10000L", parallel = TRUE, figure.title = "Figure", graph = TRUE, options=NULL) {
    if(!is.character(file)) {return(warning("the parameter 'file' has to be a character chain giving the name of the .Rmd file to write in"))}
    
    if(!is.numeric(ncp) & !is.null(ncp)) {return(warning("the argument 'ncp' must be numeric"))}
    if(!is.null(ncp)) {if(ncp < 0) {return(warning("the argument 'ncp' must be positive"))}}
    
    if(!is.logical(graph)) {return(warning("the argument 'graph' must be logical"))}
    
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    ref = eigenRef(res, dim = NULL, q = q, time = time, parallel = parallel)
    
    Refaxe = ref$inertia[1] * 100
    Refplan = ref$inertia[2] * 100
    
    Qaxe = res$eig[1,3]
    Qplan = res$eig[2,3]
    
    q = ref$quantile
    try = ref$datasets
    
    writeRmd(gettext("The inertia of the first dimensions shows if there are strong relationships between variables and suggests the number of dimensions that should be studied",domain="R-FactoInvestigate"), end = ".\n\n", file = file)
    
    switch(analyse,
           PCA = {
             writeRmd(gettext("The first two dimensions of",domain="R-FactoInvestigate"), " ", gettext("analyse",domain="R-FactoInvestigate"), " ", gettext("express",domain="R-FactoInvestigate"), " **", round(Qplan, 2), "%** ", gettext("of the total dataset inertia",domain="R-FactoInvestigate"), " ; ",
                      gettext("that means that",domain="R-FactoInvestigate"), " ", round(Qplan, 2), "% ", gettext("of the individuals (or variables) cloud total variability is explained by the plane",domain="R-FactoInvestigate"), end = ".\n", sep = "", file = file)
             
             ref.text = paste("(", gettext("the reference value is the",domain="R-FactoInvestigate"), " ", q, gettext("-quantile of the inertia percentages distribution obtained by simulating",domain="R-FactoInvestigate"), " ", 
                              try, " ", gettext("data tables of equivalent size on the basis of a normal distribution",domain="R-FactoInvestigate"), ").", sep = "")
           },
           
           CA = {
             writeRmd(gettext("The first two dimensions of",domain="R-FactoInvestigate"), " ", gettext("analyse",domain="R-FactoInvestigate"), " ", gettext("express",domain="R-FactoInvestigate"), " **", round(Qplan, 2), "%** ", gettext("of the total dataset inertia",domain="R-FactoInvestigate"), " ; ",
                      gettext("that means that",domain="R-FactoInvestigate"), " ", round(Qplan, 2), "% ", gettext("of the rows (or columns) cloud total variability is explained by the plane",domain="R-FactoInvestigate"), end = ".\n", sep = "", file = file)
             
             ref.text = paste("(", gettext("the reference value is the",domain="R-FactoInvestigate"), " ", q, gettext("-quantile of the inertia percentages distribution obtained by simulating",domain="R-FactoInvestigate"), " ", 
                              try, " ", gettext("data tables of equivalent size on the basis of a uniform distribution",domain="R-FactoInvestigate"), ").", sep = "")
           },
           
           CaGalt = {},
           
           MCA = {
             writeRmd(gettext("The first two dimensions of",domain="R-FactoInvestigate"), " ", gettext("analyse",domain="R-FactoInvestigate"), " ", gettext("express",domain="R-FactoInvestigate"), " **", round(Qplan, 2), "%** ", gettext("of the total dataset inertia",domain="R-FactoInvestigate"), " ; ",
                      gettext("that means that",domain="R-FactoInvestigate"), " ", round(Qplan, 2), "% ", gettext("of the individuals (or variables) cloud total variability is explained by the plane",domain="R-FactoInvestigate"), end = ".\n", sep = "", file = file)
             
             ref.text = paste("(", gettext("the reference value is the",domain="R-FactoInvestigate"), " ", q, gettext("-quantile of the inertia percentages distribution obtained by simulating",domain="R-FactoInvestigate"), " ", 
                              try, " ", gettext("data tables of equivalent size on the basis of a uniform distribution",domain="R-FactoInvestigate"), ").", sep = "")
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
    
    
    if(Qplan > Refplan) {
      if(Qplan > Refplan + 20) {
        ref.comp = paste(gettext("This value is strongly greater than the reference value that equals",domain="R-FactoInvestigate"), " **", round(Refplan, 2), 
                         "%**, ", gettext("the variability explained by this plane is thus highly significant",domain="R-FactoInvestigate"), sep = "")
      } else {
        ref.comp = paste(gettext("This value is greater than the reference value that equals",domain="R-FactoInvestigate"), " **", round(Refplan, 2), 
                         "%**, ", gettext("the variability explained by this plane is thus significant",domain="R-FactoInvestigate"), sep = "")
      }
      
      if(Qplan >= 95) {
        writeRmd(gettext("This percentage is particularly high and thus the first plane perfectly represents the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
        writeRmd(ref.comp, file = file)
        writeRmd(ref.text, file = file, end = "\n\n")
        writeRmd(gettext("From these observations, it is absolutely not necessary to interpret the next dimensions",domain="R-FactoInvestigate"), end = ".\n", file = file)
      } else {
        if(Qplan >= 90) {
          writeRmd(gettext("This percentage is very high and thus the first plane represents very well the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
          writeRmd(ref.comp, file = file)
          writeRmd(ref.text, file = file, end = "\n\n")
          writeRmd(gettext("From these observations, it is not necessary to interpret the next dimensions",domain="R-FactoInvestigate"), end = ".\n", file = file)
        } else {
          if(Qplan >= 70) {
            writeRmd(gettext("This percentage is high and thus the first plane represents an important part of the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
            writeRmd(ref.comp, file = file)
            writeRmd(ref.text, file = file, end = "\n\n")
            writeRmd(gettext("From these observations, it is probably not useful to interpret the next dimensions",domain="R-FactoInvestigate"), end = ".\n", file = file)
          } else {
            if(Qplan >= 50) {
              writeRmd(gettext("This percentage is relatively high and thus the first plane well represents the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
              writeRmd(ref.comp, file = file)
              writeRmd(ref.text, file = file, end = "\n\n")
              writeRmd(gettext("From these observations, it should be better to also interpret the dimensions greater or equal to the third one",domain="R-FactoInvestigate"), end = ".\n", file = file)
            } else {
              if(Qplan >= 30) {
                writeRmd(gettext("This is an intermediate percentage and the first plane represents a part of the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
                writeRmd(ref.comp, file = file)
                writeRmd(ref.text, file = file, end = "\n\n")
                writeRmd(gettext("From these observations, it may be interesting to consider the next dimensions which also express a high percentage of the total inertia",domain="R-FactoInvestigate"), end = ".\n", file = file)
              } else {
                if(Qplan >= 20) {
                  writeRmd(gettext("This is a small percentage and the first plane just represents a part of the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
                  writeRmd(ref.comp, file = file)
                  writeRmd(ref.text, file = file, end = "\n\n")
                  writeRmd(gettext("From these observations, it is interesting to consider the next dimensions which also express a high percentage of the total inertia",domain="R-FactoInvestigate"), end = ".\n", file = file)
                } else {
                  writeRmd(gettext("This is a very small percentage and the first plane represents a small part of the data variability",domain="R-FactoInvestigate"), end = ".\n", file = file)
                  writeRmd(ref.comp, file = file)
                  writeRmd(ref.text, file = file, end = "\n\n")
                  writeRmd(gettext("From these observations, it is important to also interpret the dimensions greater or equal to the third one",domain="R-FactoInvestigate"), end = ".\n", file = file)
                }
              }
            }
          }
        }
      }
    } else {
      writeRmd(gettext("The inertia observed on the first plane is smaller than the reference value that equals",domain="R-FactoInvestigate"), 
               " **", round(Refplan, 2), "%**, ", gettext("therefore low in comparison",domain="R-FactoInvestigate"), sep = "", file = file)
      writeRmd(ref.text, file = file)
      
      if(Qaxe > Refaxe)  {
        writeRmd(gettext("However, the inertia related to the first dimension is greater than the reference value",domain="R-FactoInvestigate"), 
                 " **", round(Refaxe, 2), "%**.", sep = "", file = file)
        writeRmd(gettext("Even if the inertia projected on the first plane is not significant, these explained by the first dimension is significant",domain="R-FactoInvestigate"), end = ".\n", file = file)
      } else {
        writeRmd(gettext("Moreover, the inertia projected on the first dimension is smaller than the reference value",domain="R-FactoInvestigate"), 
                 " **", round(Refaxe, 2), "%**.", sep = "", file = file)
        writeRmd(gettext("The variability expressed by the",domain="R-FactoInvestigate"), gettext("analyse",domain="R-FactoInvestigate"), gettext("is thus **not** significant",domain="R-FactoInvestigate"), end = ".\n", file = file)
        
      }
    }
    
    if(graph) {
      barplot(res$eig[,2], names.arg = 1:nrow(res$eig), main = paste(gettext("Decomposition of the total inertia on the components of the ",domain="R-FactoInvestigate"), gettext("analyse",domain="R-FactoInvestigate")))
    }
    
    writeRmd(file = file)
    writeRmd("par(mar = c(2.6, 4.1, 1.1, 2.1))\nbarplot(res$eig[,2], names.arg = 1:nrow(res$eig))", file = file, end = "\n\n",
             start = TRUE, stop = TRUE, options = options)
    
    writeRmd("**", figure.title, " - ", gettext("Decomposition of the total inertia on the components of the ",domain="R-FactoInvestigate"), gettext("analyse",domain="R-FactoInvestigate"), "**", sep = "", file = file)
    
    
    if(res$eig[1, 2] > 80) {
      writeRmd("*", gettext("The first factor is largely dominant: it expresses itself",domain="R-FactoInvestigate") , " ", round(res$eig[1, 2], 2), 
               "% ", gettext("of the data variability",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
      writeRmd("*", gettext("Note that in such a case, the variability related to the other components might be meaningless, despite of a high percentage",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
    } else {
      if(res$eig[1, 2] > 50 & res$eig[1, 2] > 2 * res$eig[2, 2]) {
        writeRmd("*", gettext("The first factor is major: it expresses itself",domain="R-FactoInvestigate"), " ", round(res$eig[1, 2], 2), 
                 "% ", gettext("of the data variability",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
        writeRmd("*", gettext("Note that in such a case, the variability related to the other components might be meaningless, despite of a high percentage",domain="R-FactoInvestigate"), ".*", sep = "", file = file)
      }
    }
    writeRmd(file = file)
    
    aleat = c(ref$inertia[1], diff(ref$inertia)) * 100
    estim.ncp = dimRestrict(res, aleat, file = file)
    
    if(is.null(ncp)) {
      if(estim.ncp == 0) {
        writeRmd(gettext("An estimation of the right number of axis to interpret suggests to not interpret the analysis at all",domain="R-FactoInvestigate"), file = file, end = ".\n")
        writeRmd(gettext("Indeed, the amount of inertia of the first axis is not higher than that obtained by the",domain="R-FactoInvestigate"), " ", q, 
                 gettext("-quantile of random distributions",domain="R-FactoInvestigate"), " (", round(res$eig[1, 3], 2), "% ", gettext("against",domain="R-FactoInvestigate"), " ", round(ref$inertia[1] * 100, 2), file = file, end = "%).\n", sep = "")
        writeRmd(gettext("This observation suggests that no axis is carrying a real information",domain="R-FactoInvestigate"), file = file, end = ".\n")
      } else if(estim.ncp == 1) {
        writeRmd(gettext("An estimation of the right number of axis to interpret suggests to restrict the analysis to the description of the first",domain="R-FactoInvestigate"), 
                 estim.ncp, gettext("axis",domain="R-FactoInvestigate"), file = file, end = ".\n")
        writeRmd(gettext("These axis present an amount of inertia greater than those obtained by the",domain="R-FactoInvestigate"), " ", q, 
                 gettext("-quantile of random distributions",domain="R-FactoInvestigate"), " (", round(res$eig[1, 3], 2), "% ", gettext("against",domain="R-FactoInvestigate"), " ", round(ref$inertia[1] * 100, 2), file = file, end = "%).\n", sep = "")
        writeRmd(gettext("This observation suggests that only this axis is carrying a real information",domain="R-FactoInvestigate"), file = file, end = ".\n")
      } else {
        writeRmd(gettext("An estimation of the right number of axis to interpret suggests to restrict the analysis to the description of the first",domain="R-FactoInvestigate"), 
                 estim.ncp, gettext("axis",domain="R-FactoInvestigate"), file = file, end = ".\n")
        writeRmd(gettext("These axis present an amount of inertia greater than those obtained by the",domain="R-FactoInvestigate"), " ", q, 
                 gettext("-quantile of random distributions",domain="R-FactoInvestigate"), " (", round(res$eig[estim.ncp, 3], 2), "% ", gettext("against",domain="R-FactoInvestigate"), " ", round(ref$inertia[estim.ncp] * 100, 2), file = file, end = "%).\n", sep = "")
        writeRmd(gettext("This observation suggests that only these axis are carrying a real information",domain="R-FactoInvestigate"), file = file, end = ".\n")
      }
      writeRmd(gettext("As a consequence, the description will stand to these axis",domain="R-FactoInvestigate"), file = file, end = ".\n")
      ncp = estim.ncp
    } else {
      writeRmd(gettext("We can observe that the first",domain="R-FactoInvestigate"), " ", estim.ncp, " ", gettext("axis present an amount of inertia greater than those obtained by the",domain="R-FactoInvestigate"),
               " ", q, gettext("-quantile of random distributions",domain="R-FactoInvestigate"), " (", round(res$eig[estim.ncp, 3], 2), "% ", gettext("against",domain="R-FactoInvestigate"), " ", round(ref$inertia[estim.ncp] * 100, 2), file = file, end = "%).\n", sep = "")
      writeRmd(gettext("Thus, a wise decision would be to restrict the description to these only axis",domain="R-FactoInvestigate"), file = file, end = ".\n")
      writeRmd(gettext("However, we choosed to describe the first",domain="R-FactoInvestigate"), ncp, gettext("axis",domain="R-FactoInvestigate"), file = file, end = ".\n")
    }
    
    return(ncp)
  }
