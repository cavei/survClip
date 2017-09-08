survivalcox <- function(coxObj, formula){
  originalCoxObj=coxObj
  coxObj <- na.omit(coxObj)
  coxRes <- survival::coxph(as.formula(formula), na.omit(coxObj))
  coxSummary <- summary(coxRes)
  zlist <- coxSummary$coefficients[,"Pr(>|z|)"]
  names(zlist) <- row.names(coxSummary$coefficients)
  pvalue <- coxSummary$logtest["pvalue"]
  return(list(pvalue=pvalue, zlist=zlist, coxObj=originalCoxObj))
}

computeDays <- function(timeTable) {
  if (NCOL(timeTable) != 2)
    stop("Data time table not matched.")

  days <- as.numeric(as.Date(timeTable[,2], format="%d/%m/%Y")-as.Date(timeTable[,1], format="%d/%m/%Y"))

  if (any(days < 0, na.rm = T)){
    stop(paste(paste(row.names(timeTable)[which(days<0)], collpase=" "), "have negative time.", sep=" "))
  }
  return(days)
}
