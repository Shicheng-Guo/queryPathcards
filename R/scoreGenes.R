#'Score Genes
#'
#'Scores each pathway with a specificed seurat object
#'
#'@param result_list list created by pathcardsQuery function
#'@param seurat_object seurat object specified by user
#'@param group_by_names specific categories user specifies the scores to be grouped by
#'
#'@return list of pathways with their entire score and pathways with their summarized scores
#'
#'@export

scoreGenes <- function(result_list, seurat_object, group_by_names) {

  pathway_names <- fixNames(names(result_list[[2]]))
  result_list[[1]]$Shared_Pathway <- fixNames(result_list[[1]]$Shared_Pathway)

  # score genes
  for (i in 1:length(pathway_names)) {
    scoring_genes <- unlist(result_list$pathway_genes[[i]])
    scoring_genes_reduced <- list(scoring_genes[which(scoring_genes %in% rownames(GetAssayData(object = seurat_object,
                                                                                               slot = "data", assay = "RNA")))])
    seurat_object <- AddModuleScore(seurat_object, features = scoring_genes_reduced, name = pathway_names[i],
                                    nbin = 24)
    print(paste("completed scoring of", pathway_names[i]))
  }

  # make them mergable taking out pathways and their scores along with other factors chosen by the user
  pathway_scores <- dplyr::select(seurat_object@meta.data, group_by_names, colnames(seurat_object@meta.data)[((ncol(seurat_object@meta.data) -
                                                                                                                 length(pathway_names)) + 1):ncol(seurat_object@meta.data)])
  # summarise them by the factors chosen by user
  summary_scores <- pathway_scores %>% group_by_at(group_by_names) %>% summarise_all(funs(mean))
  # Scoring adds number '1' to the end of each pathway, so I removed them to make it mergable
  colnames(summary_scores)[(ncol(summary_scores) - length(pathway_names) + 1):ncol(summary_scores)] <- gsub(".{1}$",
                                                                                                            "", colnames(summary_scores)[(ncol(summary_scores) - length(pathway_names) + 1):ncol(summary_scores)])
  summary_scores <- gather(data = summary_scores, key = Shared_Pathway, value = Score_Average, pathway_names)

  # merge
  summarized_result <- merge(result_list[[1]], summary_scores)
  return(list(full_data = pathway_scores, summarized_data = summarized_result))
}
