#'Annotate Receptor Ligand Pair
#'
#'Subfunction for the function receptor.ligand.annotation.from.object
#'
#'@export

receptor.ligand.pair.annotation<-function(Mat1, Mat2, L.R.Mat, group.by)
{
  #Counter for Row of Final Data Frame
  k=1

  #Final Data Frame to Return
  interactions.mat<-as.data.frame(matrix(ncol=2, nrow=0))
  colnames(interactions.mat) <- c("Receptor","Ligand")

  #Loop through Receptor DF to Find Possible Receptor / Ligand Interactions based on thresholded receptors and full R/L interactions list
  for(i in 1:nrow(Mat1))  {

    Cell.Type.1.tmp<-as.character(Mat1[[group.by]][i])
    Receptor.tmp<-as.character(Mat1$variable[i])
    receptor.interactions.tmp<-L.R.Mat[which(L.R.Mat$Receptor == Receptor.tmp),]
    #print(receptor.interactions.tmp)
    #Loop through thesholded Ligand DF to see which of the current receptor and its possible interactions match

    for(j in 1:nrow(Mat2))  {

      Ligand.tmp<-as.character(Mat2$variable[j])
      Cell.Type.2.tmp<-as.character(Mat2[[group.by]][j])
      ligand.pair.tmp<-receptor.interactions.tmp[which(receptor.interactions.tmp$Ligand == Ligand.tmp),]
      #print(ligand.pair.tmp)

      if(nrow(ligand.pair.tmp)!=0)  {

        ligand.pair.confirm<-ligand.pair.tmp
        interactions.mat[k,1]<-paste(Cell.Type.1.tmp, ligand.pair.confirm$Receptor, sep="_")
        interactions.mat[k,2]<-paste(Cell.Type.2.tmp, ligand.pair.confirm$Ligand, sep="_")
        #print(c(paste(Cell.Type.1.tmp, ligand.pair.confirm$Receptor, sep="_"),paste(Cell.Type.2.tmp, ligand.pair.confirm$Ligand, sep="_")))
        k=k+1

      } else {

        #print("Skip")

      }
    }
  }
  return(interactions.mat)
}

#' Nonzero
#'
#' Another auxillary function for receptor.ligand.annotation.from.object

nonzero <- function(x){sum(x > 0)}
