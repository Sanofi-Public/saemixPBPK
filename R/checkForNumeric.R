
####################################################################################
####			               checkForNumeric   			####
####################################################################################

#' Function to verify which column in a dataframe are numeric
#' 
#' 
#' @param myD a dataframe to be checked.
#' @param verbose Logical.If \code{TRUE}, the result provided by the function is a text. If \code{verbose=TRUE}, the result will be a number indicating the number of column containing NAs. In this last case, the outcome of the function may be used by other functions.
#' @return Print result on screen (\code{verbose=TRUE}) or return the number of column
#' @author Donato Teutonico.
#' 
#' @examples
#' 
#' 
#' checkForNumeric(d)     # By default verbose=TRUE
#' checkForNumeric(d,verbose=FALSE)
#
#'
#' @export checkForNumeric 


checkForNumeric <-
function(myD,verbose=TRUE){

                              nCol <- c(1:ncol(myD))
                              c <- 0
                              for (i in nCol){
                                         myClass <- class(myD[,nCol[i]])
                                         if(myClass!="numeric" & myClass!="integer" ){ c<-c+1
                                         if (verbose) print(paste("Values in column ",i," (",names(myD)[i],")"," are not numeric",sep=""))}
                                            }
                              if(c==0&verbose) return("All the values in the dataset are numeric")
                              if (!verbose) return(c)
}

