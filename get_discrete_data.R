get_discrete_data <- function (data,meth)
{
  data_discrete=matrix(nrow=nrow(data),ncol=ncol(data))
  col_num = ncol(data)
  for (iter in 1:col_num){
    temp <- arules::discretize(data[,iter],method= meth,labels = c(-1,0,1))
    data_discrete[,iter] <- as.numeric(levels(temp))[temp]
  }
  
  return(data_discrete)
}