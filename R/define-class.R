clipper <- setClass("clipper", package = "clipper",
                    slots = c(table = "data.frame"),
                    contains = "data.frame"
                    )


setMethod("show",
          signature = "clipper",
          definition = function(object) {
            cat("clipper results summary\n")
            cat("tatal number of subPath detected: ", NROW(object@table), sep="")
            invisible(NULL)
          })
