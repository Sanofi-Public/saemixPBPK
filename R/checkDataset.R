
####################################################################################
####			               checkDataset   			####
####################################################################################

#' Function to check the dataset and detect possible issues for its use with saemixPBPK
#' 
#' 
#' @param data a dataframe to be checked.
#' @return Print result on screen the result of the check
#' @author Donato Teutonico.
#' 
#' @examples
#' 
#' 
#' checkDataset(data)   
#
#'
#' @export checkDataset 


checkDataset <-
function(data){
  checkForNumeric(data, verbose=TRUE)
  checkForNA(data, verbose=TRUE)
  checkID(data)
}
            

checkID <- function(data){
          header <- colnames(data)
          header <- tolower(header)
          id_col <- grepl("^id$", header)
          unsorted <- is.unsorted(data[,id_col]) 
          if(unsorted) print("IDs are unsorted, they should be ordered monotonically")
          else print("IDs are sorted monotonically")
        
          myIDs <- sort(as.numeric(as.character(unique(data[,id_col]))))
          if (length(myIDs)==0) stop("The IDs column contains undefine characters")
          result <- unique(diff(myIDs))
          if(length(result)>=2) print("Ids are not consecutives, use the function renumberID")
          else print("IDs are consecutives")
}
