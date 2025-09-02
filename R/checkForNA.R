
####################################################################################
####			               checkForNA                                        			####
####################################################################################

#' Function to verify the presence of NAs in dataframes
#' 
#' 
#' @param d a dataframe to be checked.
#' @param verbose Logical.If \code{TRUE}, the result provided by the function is a text. If \code{verbose=TRUE}, the result will be a number indicating the number of column containing NAs. In this last case, the outcome of the function may be used by other functions.
#' @return Print result on screen (\code{verbose=TRUE}) or return the number of column
#' @author Donato Teutonico.
#' 
#' @examples
#' 
#' 
#' checkForNA(d)     # By default verbose=TRUE
#' checkForNA(d,verbose=FALSE)
#
#'
#' @export checkForNA




checkForNA <-
function(d,verbose=TRUE){
                         c <-0
                         if(any(is.na(d))){
                                    w1 <- which(sapply(d,function(x) any(is.na(x))))
                                    c<-c+length(w1)
                                    if (verbose) {
                                                 print("The dataset contains NA values")
                                                 print(w1)
                                                 } else return(c)
                                           } 
                         else if (verbose) print("The dataset does not contain NA values")
                              else return(c)
}

