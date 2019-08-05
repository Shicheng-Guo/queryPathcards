#'Query Pathcards Website
#'
#'Finds signaling pathways for each gene and then finds gene list for each signaling pathway.
#'
#'@param gene_list list of genes.
#'
#'@return dataframe with original genes and corresponding signaling pathways and a list of signaling pathways and corresponding gene lists.
#'
#'@export

queryPathcards <- function(gene_list) {
  # checking if input was only ligands
  if (dim(gene_list)[2] == 1) {
    gene_list[1] <- toupper(gene_list[1])
    pathway_matrix_main <- data.frame(Ligand = 0, Shared_Pathway = 0)
    for (i in 1:dim(gene_list)[1]) {
      gene_name <- gene_list[i]
      extracted_pathways <- extractPathway(gene_name)
      pathway_matrix_temp <- matrix("", nrow = length(extracted_pathways), ncol = 2)
      colnames(pathway_matrix_temp) <- c("Ligand", "Shared_Pathway")
      if (length(extracted_pathways) == 0) {
        pathway_matrix_temp[1, 1] <- gene_name
        pathway_matrix_temp[1, 2] <- "NA"
      } else {
        for (j in 1:length(extracted_pathways)) {
          pathway_matrix_temp[j, 1] <- gene_name
          pathway_matrix_temp[j, 2] <- extracted_pathways[j]
        }
      }
      pathway_matrix_main <- rbind(pathway_matrix_main, pathway_matrix_temp)
    }
  } else {
    # the url only works if the genes are all CAPS
    gene_list[, 1] <- toupper(gene_list[, 1])
    gene_list[, 2] <- toupper(gene_list[, 2])
    pathway_matrix_main <- data.frame(Ligand = 0, Receptor = 0, Shared_Pathway = 0)
    for (i in 1:dim(gene_list)[1]) {
      ligand_name <- gene_list[i, 1]
      receptor_name <- gene_list[i, 2]
      extracted_pathways_1 <- extractPathway(ligand_name)
      extracted_pathways_2 <- extractPathway(receptor_name)
      shared_pathways <- intersect(extracted_pathways_1, extracted_pathways_2)
      if (length(shared_pathways) == 0) {
        pathway_matrix_temp <- matrix("", nrow = 1, ncol = 3)
        colnames(pathway_matrix_temp) <- c("Ligand", "Receptor", "Shared_Pathway")
        pathway_matrix_temp[1, 1] <- ligand_name
        pathway_matrix_temp[1, 2] <- receptor_name
        pathway_matrix_temp[1, 3] <- "NA"
      } else {
        pathway_matrix_temp <- matrix("", nrow = length(shared_pathways), ncol = 3)
        colnames(pathway_matrix_temp) <- c("Ligand", "Receptor", "Shared_Pathway")
        for (j in 1:length(shared_pathways)) {
          pathway_matrix_temp[j, 1] <- ligand_name
          pathway_matrix_temp[j, 2] <- receptor_name
          pathway_matrix_temp[j, 3] <- shared_pathways[j]
        }
      }
      print(as.character(paste(ligand_name, receptor_name)))
      pathway_matrix_main <- rbind(pathway_matrix_main, pathway_matrix_temp)
    }
  }
  url_template_pathway <- "https://pathcards.genecards.org/card/"
  pathways <- unique(as.character(pathway_matrix_main[-1, ]$Shared_Pathway))

  if (identical(which(pathways %in% "NA"), integer(0)) == FALSE) {
    pathways <- pathways[-which(pathways %in% "NA")]
  }

  pathway_list <- vector("list", length(pathways))
  names(pathway_list) <- pathways
  pathways <- tolower(pathways)

  # removing unnecessary symbols from pathway names
  for (i in 1:length(pathways)) {
    name_simplified <- fixNames(pathways[i])
    url_3 <- paste(url_template_pathway, name_simplified, sep = "")
    gene_data <- url_3 %>% read_html() %>% html_nodes(".geneslist") %>% html_text() %>% str_remove_all("\r\n") %>%
      str_remove_all("via the multiplicity of each gene in the constituent pathways.") %>% str_remove_all("\\*Darkness represents the genes rank within the SuperPath,") %>%
      trimws() %>% as.data.frame()
    if (nrow(gene_data) >= 2) {
      gene_data <- gene_data[-1, ]
    } else {
      gene_data <- gene_data[1, ]
    }
    pathway_list[[i]] <- gene_data %>% as.character() %>% strsplit("\\s+")
    print(as.character(pathways[i]))
  }
  return(list(ligand_receptor_pathways = pathway_matrix_main[-1, ], pathway_genes = pathway_list))
}
