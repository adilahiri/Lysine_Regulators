get_expected <- function (data){
  Names<- data$Names
  data$Names <- NULL
  
  colnames(data) <- c("inhibit","dormant","active")
  
  data<-apply(data,1,function(x) c(x[1]/sum(x),x[2]/sum(x),x[3]/sum(x)))
  
  data_trans<-as.data.frame(t(as.matrix(data)))
  data_trans$Names <- Names
  return(data_trans)
  
}