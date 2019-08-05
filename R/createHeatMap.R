#'Create Heat Map
#'
#'Creates Heat Map using ggplot2 with the scores.
#'
#'@param score_df data frame of scores of each signaling pathway.
#'@param y.var y variable specified by the user.
#'@param x.var x variable specified by the user.
#'
#'@return Heat Map of the scores.
#'
#'@export

createHeatMap <- function(score_df, y.var, x.var) {
  return(ggplot(score_df, aes_string(y = y.var, x = x.var, fill = "Score_Average")) + geom_tile() + theme_classic() +
           scale_fill_gradient2(low = "steelblue", mid = "gray90", high = "darkred") + theme(axis.text.x = element_text(angle = 90,
                                                                                                                        hjust = 1)))
}
