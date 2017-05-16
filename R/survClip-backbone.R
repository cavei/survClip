survivalcox <- function(coxObj, formula){
  coxRes <- survival::coxph(as.formula(formula), na.omit(coxObj))
  coxSummary <- summary(coxRes)
  zlist <- coxSummary$coefficients[,"Pr(>|z|)"]
  names(zlist) <- row.names(coxSummary$coefficients)
  pvalue <- coxSummary$logtest["pvalue"]
  return(list(pvalue=pvalue, zlist=zlist))
  # return(pvalue)
}

fullsurvivalcox <- function(coxObj, formula){
  coxObj <- na.omit(coxObj)
  coxRes <- survival::coxph(as.formula(formula), coxObj)
  # sink("/dev/null")
  # suppressWarnings(stepCoxRes <- stats::step(coxRes))
  # if (!is.null(summary(stepCoxRes)$coefficients))
  #   coxRes <- stepCoxRes
  # sink(NULL)
  coxSummary <- summary(coxRes)
  zlist <- coxSummary$coefficients[,"Pr(>|z|)"]
  names(zlist) <- row.names(coxSummary$coefficients)
  pvalue <- coxSummary$logtest["pvalue"]
  return(list(pvalue=pvalue, zlist=zlist))
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

medianSurv <- function(exp, min){
  m <- median(exp)
  class <- rep("high", length(exp))
  names(class) <- names(exp)
  class[exp < m] <- "low"
  hcount <- sum(class=="high", na.rm=T)
  lcount <- sum(class=="low", na.rm=T)
  if ((hcount <= min) || (lcount <=min)) {
    q <- 0
  } else {
    q <- 1
  }
  return(list(class=class, qual=q))
}
