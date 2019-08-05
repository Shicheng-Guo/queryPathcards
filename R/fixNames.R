#'Fix Pathway Names
#'
#'Delete any unnecessary symbols within pathway names.
#'
#'@param names_list list of pathway names.
#'
#'@return modified list of pthway names.
#'
#'@export

fixNames <- function(names_list) {
  names_list <- gsub(" ", "_", names_list)
  names_list <- gsub("/", "", names_list)
  names_list <- gsub("[+]", "", names_list)
  names_list <- gsub("[:]", "", names_list)
  names_list <- gsub("[,]", "", names_list)
  names_list <- gsub("[.]", "", names_list)
  return(names_list)
}
