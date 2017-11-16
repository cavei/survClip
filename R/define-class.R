survCliques <- setClass("survCliques", package = "survClip",
                    slots = c(alphas  = "numeric",
                              zlist   = "list", 
                              cliques = "list",
                              coxObjs = "list",
                              cliquesLoadings = "list",
                              cliquesExpr = "list"),
                    contains = "list"
)

setMethod("show",
          signature = "survCliques",
          definition = function(object) {
            sthis <- seq_len(min(length(object@alphas), 3))
            sthis <- order(object@alphas)[sthis]
            
            sigCliquesIdx = which(object@alphas <= 0.05)
            
            for (i in sthis) {
              cat(paste0("Cliques ",i, ": pvalue ", object@alphas[i], "\n"))
              cat("Clique is composed by the followings:\n")
              cat(paste(object@cliques[[i]], collapse=", "))
              cat("\n-+-\n")
            }
            
            if (length(sthis) < length(sigCliquesIdx)) {
              cat(paste0("There are other ", length(sigCliquesIdx)-length(sthis), " cliques with pvalue <= 0.05"))
            }
              
            invisible(NULL)
          })


survPath <- setClass("survPath", package = "survClip",
                     slots = c(pvalue = "numeric",
                               zlist = "numeric",
                               coxObj = "data.frame",
                               loadings = "matrix",
                               method  = "character"),
                     contains = "list"
)

setMethod("show",
          signature = "survPath",
          definition = function(object) {
            cat(paste0("Pathway processed with ", object@method, " method\n   pvalue: ", object@pvalue, "\n"))
            invisible(NULL)
          })

