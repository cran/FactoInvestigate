dimActive <-
function(res) {
    analyse = whichFacto(res)
    if(!analyse %in% c("PCA", "CA", "CaGalt", "MCA", "MFA", "DMFA", "FAMD", "GPA", "HCPC"))
    {return(warning("the parameter 'res' has to be an object of class 'PCA', 'CA', 'CaGalt', 'MCA', 'MFA', 'DMFA', 'FAMD', 'GPA' or 'HCPC'"))}
    param = getParam(res)
    
    switch(analyse,
           PCA = {
             c(param$ind, param$var)
           },
           
           CA = {
             c(param$row, param$col)
           },
           
           CaGalt = {},
           
           MCA = {
             c(param$ind, length(param$modalites))
           },
           
           MFA = {},
           
           HMFA = {},
           
           DMFA = {},
           
           FAMD = {},
           
           GPA = {},
           
           HCPC = {})
  }
