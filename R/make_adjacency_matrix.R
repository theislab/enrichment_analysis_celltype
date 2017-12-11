# The below script updates the adjacency matrix
# file - path of the file from which the cell-type annotations are to be read in. The first column
# should correspond to the cell-type (source) ID

make_adjacency_matrix <- function(file){
  
  require(xlsx)
  
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