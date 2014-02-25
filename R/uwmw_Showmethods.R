# Show methods

setMethod('show','uwmwRes',
          function(object){
            ng <- length(object)
            if(ng==0){
              cat("Empty uwmwRes object.")
              return(invisible(NULL))
            }
            typestring <- switch(object@type,"H"="housekeeping","O"="overall")
            cat("unified Wilcoxon-Mann-Whitney test \n")
            cat("with",typestring,"normalization\n")
            cat("number of features:",ng,"\n")
            if(type(object)=='H')
              cat("number of housekeeping features:", 
                  length(housekeeping(object)),"\n\n")
            cat(groups.message(groupinfo(object)),"\n\n")
            
            invisible(NULL)
            
          })


# For uwmwEstimate
setMethod('show','uwmwEstimate',
          function(object){
            ng <- length(object)
            if(ng==0){
              cat("Empty uwmwEstimate object.")
              return(invisible(NULL))
            }
            typestring <- switch(type(object),"H"="housekeeping","O"="overall")
            cat("unified Wilcoxon-Mann-Whitney test \n")
            cat("with",typestring," normalization\n")
            cat("number of genes:",ng,"\n")
            if(type(object)=='H')
              cat("number of housekeeping features:", length(housekeeping),"\n")
            cat("Confidence limits represent a",object@confint*100,"% confidence interval.\n")
            if(!is.na(object@refest))
              cat("Reference estimate:",
                  object@refest,"(",object@refll,"-",object@reful,")\n\n")
            
            cat(groups.message(groupinfo(object)),"\n\n")
            
            print(as.matrix(object))
            
            invisible(NULL)
            
          })