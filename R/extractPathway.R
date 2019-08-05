#'Extract Signaling Pathways
#'
#'Queries the PathCards website to find signaling pathways for each gene.
#'
#'@param gene any gene.
#'
#'@return list of pathways.
#'
#'@export

extractPathway <- function(gene) {
  url_template <- "https://pathcards.genecards.org/Search/Results?query="
  url_1 <- paste(url_template, gene, sep = "")

  # read html and get rid of empty columns.  Found nodes 'td' and 'a' using a chrome extension called
  # SelectorGadget.
  initial_data_html <- url_1 %>% read_html() %>% html_nodes("td") %>% html_nodes("a")

  # checking if gene didn't have a page for it
  if (length(html_text(initial_data_html)) == 0) {
    extracted_pathways <- "NA"
    return(extracted_pathways)

  } else {

    extracted_pathways <- initial_data_html %>% html_text()
    return(extracted_pathways)
  }
}
