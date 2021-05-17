get_dirichlet <- function(data){
  prior_para <- c(1,1,1)  # set the prior to dirichlet (1,1,1) in the order of -1,0,1
  
  parameter_matrix <- matrix(ncol=3)
  node_param_list <- list()
  colnames(parameter_matrix)<- c("alpha_inhibit","alpha_dormant","alpha_active")
  
  node_param<-matrix(ncol=3)
  
  parent_matrix <-NULL
  Name_Vector <- NULL
  Nodes <-ncol(data)
  for (iter in 1: Nodes){
    if (iter >=1 & iter <=5) {# Nodes A,B,C,D,E
      parent_matrix <-NULL
    }
    else if(iter ==6){ # Node F
      parent_matrix<-cbind(data[,c("A","B","C","D","E")])
    }
    else if (iter == 7 | iter ==8){ # Nodes G and H
      parent_matrix<-cbind(data[,c("D","E","F")])
    }
    else if (iter ==9 | iter ==10){ # Nodes I and J
      parent_matrix <- cbind(data[,c("G","H")])
    }
    else if (iter == 11 | iter == 12) { # Nodes K and L 
      parent_matrix <- cbind(data[,c("I","J")])
    }
    else if (iter == 13 ) { # Node M
      parent_matrix <- cbind(data[,c("K","L")])
    } 
    else if (iter == 14){ # Node N
      parent_matrix <- cbind(data[,"M"])
    }
    node_param<-rbind( get_param(prior=prior_para, node=data[,iter],parent=parent_matrix))
    parameter_matrix<-rbind(parameter_matrix,node_param)
    Name_Vector <- c(Name_Vector,rep(colnames(data)[iter],dim(node_param)[1]))
  }
  parameter_matrix<-as.data.frame(parameter_matrix)
  parameter_matrix<-parameter_matrix[-1,]
  parameter_matrix$Names <- Name_Vector
  return(parameter_matrix)
  
}