enrich <- function(clusters, 
                   pval = 0.01, 
                   lfc = 1, 
                   annoDB = 'xCell', 
                   plot = TRUE){
  
  require(xlsx)
  require(ggplot2)
  
  cat(paste("Annotation database selected:",annoDB,'\n'))
  
  cat('Loading cell type signatures annotation file \n')
  # After cloning the git repo it should be saved in the working directory
  load(file = './enrich/Data/adj_mat.RData')

  # Initializing list object, adding metadata and results from Fisher's test
  enrich_object <- list()
  enrich_object$annotation_adj_matrix <- adj.mat
  
  # Reading in cluster data
  clusters <- read.csv(file = clusters, header = TRUE, row.names = 1)
  
  # Selecting significantly differentially expressed cluster markers
  clusters <- clusters[which(clusters[,1] <= pval & abs(clusters[,2] >= lfc)),]

  rownames(clusters) <- toupper(rownames(clusters))
  enrich_object$clusters <- clusters
  
  cat('Computing Fishers Exact Test to evaluate statistical significance of overlap \n')
  no_of_clusters <- unique(clusters[,3])

  matrix_fisher <- data.frame()
  
  for(i in no_of_clusters){
    
    v <- sapply(colnames(adj.mat), FUN = function(x){
      
      present_cluster <- as.character(rownames(clusters)[which(clusters[,3] == i)])
      present_anno <- names(which(adj.mat[,x] == 1))
      intersection <- intersect(present_anno,present_cluster)
      cont_tbl <- matrix(c(abs(length(rownames(adj.mat)) - 
                                 length(unique(union(present_anno,present_cluster)))),
                           length(present_anno) - length(intersection),
                           length(present_cluster) - length(intersection),
                           length(intersection)), ncol = 2)
      rownames(cont_tbl) <- c('not_in_anno','in_anno')
      colnames(cont_tbl) <- c('not_in_cluster','in_cluster')

      # Performing Fisher's test
      res <- fisher.test(x = cont_tbl, alternative = 'greater')

      # Saving p.value, odds ratio, common genes between cluster markers 
      # and celltype signature, cluster ID, no. of genes in annotation and
      # no. of markers in cluster
      vec <- c(res$p.value, 
               res$estimate, 
               length(present_anno), 
               length(present_cluster), 
               ifelse(test = length(intersection) > 1, 
                      paste(intersection,collapse = '|'), intersection), 
               i)
      return(vec)
      
    })
    
    matrix_fisher <- rbind(matrix_fisher,t(v))
    
  }
  colnames(matrix_fisher) <- c('p.value', 'estimate.of.odds.ratio', 'present.annotation', 
                               'present.cluster', 'genes.in.intersection', 'cluster.ID')
  
  significant_fisher <- data.frame()
  
    for(i in unique(matrix_fisher$cluster.ID)){

      matrix_fisher_sub <- matrix_fisher[which(matrix_fisher$cluster.ID == i),]
      matrix_fisher_sub$p.value <- sapply(matrix_fisher_sub$p.value, as.vector)
      matrix_fisher_sub$p.value <- as.vector(sapply(matrix_fisher_sub$p.value, as.numeric))
      # Saving top 3 enriched pathways
      temp <- matrix_fisher_sub[order(matrix_fisher_sub$p.value, decreasing = FALSE),][c(1:3),]
      significant_fisher <- rbind.data.frame(significant_fisher, temp)
    }
  significant_fisher <- cbind.data.frame(significant_fisher, 
                                         rank = rep(c(1,2,3),nrow(significant_fisher)/3))
  
  enrich_object$matrix_fisher <- matrix_fisher
  enrich_object$significant_fisher <- significant_fisher

  if(plot == TRUE){
    
    g <- ggplot(significant_fisher, 
            aes(x = cluster.ID, y = -log10(p.value), fill = as.factor(rank))) + 
    geom_bar(stat = 'identity', colour = "black", position = position_dodge(width = 0.8),
             width = 0.5) + 
    coord_flip() + 
    geom_text(aes(label = rownames(significant_fisher)), size = 2, 
              position = position_dodge(width = 0.8), hjust = -0.25) + 
    ggtitle("Most significantly enriched cluster cell-type from Fisher's exact test") + 
    theme(plot.title = element_text(hjust = 0.5))
  
  print(g)
  }
  return(enrich_object)

}