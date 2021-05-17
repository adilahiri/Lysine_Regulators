get_nlargest <- function (data_mat,n){
  nlargest_index <- order(data_mat,decreasing= TRUE)[1:n]
  ind <- matrix(nrow=n,ncol=2)
  names<- NULL
  for (iter in 1:n){
    ind[iter,] <- which(data_mat==data_mat[nlargest_index[iter]],arr.ind = TRUE)
    names[iter]<- rownames(data_mat)[ind[iter]]
  }
  df <- as.data.frame(ind)
  df$values<- data_mat[nlargest_index]
  df$Names <-names
  return(df)
}