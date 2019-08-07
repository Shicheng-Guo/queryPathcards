#'Receptor Ligand Annotation From Seurat Object
#'
#'Uses the seurat object to identify certain ligand-receptor pairs who have expression levels higher than the ones set by the user.
#'
#'@param seurat.object seurat object set by the user.
#'@param group.by variables to group by set by the user.
#'@param receptor.ligand.mat the list of receptors and ligands.
#'@param mean.expression.threshold mean expression threshold set by the user.
#'@param percent.expression.threshold percent expression threshold set by the user.
#'
#'@return list of ligand-receptor pairs with expression levels.
#'
#'@export

receptor.ligand.annotation.from.object <- function(seurat.object, group.by, receptor.ligand.mat, mean.expression.theshold=0.1, percent.expression.threshold=.2)
{

  unique_receptors.tmp<-unique(as.character(receptor.ligand.mat$Receptor))
  unique_ligands.tmp<-unique(as.character(receptor.ligand.mat$Ligand))

  # Receptors

  receptor.matrix.subset.tmp <- FetchData(seurat.object, c(group.by,unique_receptors.tmp[which(unique_receptors.tmp %in% rownames(seurat.object@assays$RNA@data))]))

  receptor.matrix.subset.melt.tmp <- melt(receptor.matrix.subset.tmp, by=group.by)

  receptor_summary.tmp <- receptor.matrix.subset.melt.tmp %>% group_by(get(group.by),variable) %>% summarise(mean_expression=mean(value), percent_expression=nonzero(value)/n())

  receptor_summary_subset.tmp <- subset(receptor_summary.tmp, mean_expression > mean.expression.theshold & percent_expression > percent.expression.threshold)

  colnames(receptor_summary_subset.tmp)<-c(group.by,"variable","mean_expression","percent_expression")

  # ligands

  ligand.matrix.subset.tmp <- FetchData(seurat.object, c(group.by,unique_ligands.tmp[which(unique_ligands.tmp %in% rownames(seurat.object@assays$RNA@data))]))

  ligand.matrix.subset.melt.tmp <- melt(ligand.matrix.subset.tmp, by=group.by)

  ligand_summary.tmp <- ligand.matrix.subset.melt.tmp %>% group_by(get(group.by),variable) %>% summarise(mean_expression=mean(value), percent_expression=nonzero(value)/n())

  ligand_summary_subset.tmp <- subset(ligand_summary.tmp, mean_expression > mean.expression.theshold & percent_expression > percent.expression.threshold)

  colnames(ligand_summary_subset.tmp) <- c(group.by,"variable","mean_expression","percent_expression")

  final.interaction.mat<-receptor.ligand.pair.annotation(receptor_summary_subset.tmp, ligand_summary_subset.tmp, receptor.ligand.mat, group.by = group.by)

  return(final.interaction.mat)
}
