# The below script updates the adjacency matrix
# file - path of the file from which the cell-type annotations are to be read in

make_adjacency_matrix <- function(file){
  
  require(xlsx)
  
  # setwd('C:/Users/user/')
  
  load(file = './enrich/data/adj_mat.RData')
  
  # updating the adjacency matrix with brain anatomy related annotations
  brain_anno_anatomical <- as.matrix(read.csv("E:/cell-type enrichment analysis/nn.4171-S13.csv", 
                                    header = TRUE, row.names = 1))
  
  dim(adj.mat) # 5079X489
  
  for(i in 1:ncol(brain_anno_anatomical)){
    
    topGenes <- head(sort(brain_anno_anatomical[,i], decreasing = TRUE), 100)
    genestoadd <- setdiff(names(topGenes),rownames(adj.mat))
    temp <- matrix(data = 0, nrow = length(genestoadd), ncol = ncol(adj.mat))
    rownames(temp) <- genestoadd
    adj.mat <- rbind(adj.mat,temp)
    temp <- matrix(data = 0, nrow = nrow(adj.mat),ncol = 1)
    colnames(temp) <- colnames(brain_anno_anatomical)[i]
    adj.mat <- cbind(adj.mat,temp)
    
    adj.mat[which(rownames(adj.mat) %in% names(topGenes)),ncol(adj.mat)] = 1
    
  }
  
  dim(adj.mat) # 5927X510
  
  write.csv(adj.mat,file = 'E:/cell-type enrichment analysis/adjacency matrix with xCell and brain annotations.csv')
  save(adj.mat,file = 'E:/cell-type enrichment analysis/adj_matv2.RData')
  
  # constructing the first adjacency matrix with xCell
  # file contains the path for the cell-type annotations for excel
  annoData <- read.xlsx(file, sheetIndex = 1)
  annoData = sapply(annoData, as.character)
  annoData[is.na(annoData)] <- NaN
  annoData <- as.data.frame(annoData)
  rownames(annoData) <- annoData[,1]
  # annoData <- annoData[,-1]
  annoData <- annoData[,-c(1,2)]
  
  cat('Constructing adjacency matrix \n')
  # Construct adjacency matrix -- genes x annotations/cell types

  signature <- c()
  
  for(i in 1:ncol(annoData)){
    signature = append(signature, as.character(annoData[,i]))
  }
  
  signature <- unique(signature)
  signature <- signature[-grep('NaN', signature)]
  adj.mat <- matrix(data = 0, nrow = length(signature), ncol = nrow(annoData))
  rownames(adj.mat) <- signature
  colnames(adj.mat) <- rownames(annoData)
  
  for(i in rownames(adj.mat)){
    temp <- which(annoData == i, arr.ind = TRUE)
    adj.mat[i, which(colnames(adj.mat) %in% rownames(temp))] = 1
  }
  
  save(adj.mat,file = 'E:/cell-type enrichment analysis/adj_mat.RData')

}