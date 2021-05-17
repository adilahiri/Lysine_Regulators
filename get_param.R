get_param <- function(prior_param,node,parent){
  alpha_dirichlet_posterior <- matrix(ncol = 3)
  # Case when there are no parents
  if(is.null(parent)){
    counts_inhibit<- as.numeric(table(node)[1])
    counts_dormant<- as.numeric(table(node)[2])
    counts_active<- as.numeric(table(node)[3])
    
    alpha_inhibit <- prior_param[1] + counts_inhibit
    alpha_dormant <- prior_param[2] + counts_dormant
    alpha_active  <- prior_param[3] + counts_active
    alpha_dirichlet_posterior <- c(alpha_inhibit,alpha_dormant,alpha_active)
    
    
  }
  # Case when there is one parent
  else if (dim(parent)[2]==1) {
    parent_evidence_table <- expand.grid(c(-1,0,1))
    for (iter in 1:nrow(parent_evidence_table)){
      alpha_inhibit=prior_param[1]+length(which(parent==parent_evidence_table[iter,]& node ==-1))
      alpha_dormant =prior_param[2] +length(which(parent==parent_evidence_table[iter,] & node ==0))
      alpha_active = prior_param[3] +length(which(parent==parent_evidence_table[iter,]& node ==1))    
      alpha_dirichlet_posterior <- rbind(alpha_dirichlet_posterior,c(alpha_inhibit,alpha_dormant,alpha_active))
      
    }
    alpha_dirichlet_posterior<-alpha_dirichlet_posterior[-1,]
  }
  # Case when there are two parents
  else if (dim(parent)[2]==2) {
    parent_evidence_table <- expand.grid(c(-1,0,1),c(-1,0,1))
    for (iter in 1:nrow(parent_evidence_table)){
      alpha_inhibit=prior_param[1]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                          parent[,2]==parent_evidence_table[iter,2]& 
                                          node ==-1))
      alpha_dormant=prior_param[2]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                          parent[,2]==parent_evidence_table[iter,2]& 
                                          node ==0))
      
      alpha_active=prior_param[3] + length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                           parent[,2]==parent_evidence_table[iter,2] &
                                          node ==0))
      alpha_dirichlet_posterior <- rbind(alpha_dirichlet_posterior,c(alpha_inhibit,alpha_dormant,alpha_active))
      
    }
    alpha_dirichlet_posterior<-alpha_dirichlet_posterior[-1,]
    
  }
  # Case when there are three parents
  else if (dim(parent)[2]==3) {
    parent_evidence_table <- expand.grid(c(-1,0,1),c(-1,0,1),c(-1,0,1))
    for (iter in 1:nrow(parent_evidence_table)){
      alpha_inhibit=prior_param[1]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                              parent[,2]==parent_evidence_table[iter,2]& 
                                              parent[,3]==parent_evidence_table[iter,3]&
                                              node ==-1))
      alpha_dormant=prior_param[2]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                              parent[,2]==parent_evidence_table[iter,2]&
                                              parent[,3]==parent_evidence_table[iter,3]&
                                              node ==0))
      
      alpha_active=prior_param[3]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                              parent[,2]==parent_evidence_table[iter,2]&
                                              parent[,3]==parent_evidence_table[iter,3]&
                                              node ==0))
      alpha_dirichlet_posterior <- rbind(alpha_dirichlet_posterior,c(alpha_inhibit,alpha_dormant,alpha_active))
      
    }
    alpha_dirichlet_posterior<-alpha_dirichlet_posterior[-1,]
    
  }
  # Case when there are five parents
  else if (dim(parent)[2]==5) {
    parent_evidence_table <- expand.grid(c(-1,0,1),c(-1,0,1),c(-1,0,1),
                                         c(-1,0,1),c(-1,0,1))
    for (iter in 1:nrow(parent_evidence_table)){
      alpha_inhibit=prior_param[1]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                              parent[,2]==parent_evidence_table[iter,2]& 
                                              parent[,3]==parent_evidence_table[iter,3]&
                                              parent[,4]==parent_evidence_table[iter,4]&
                                              parent[,5]==parent_evidence_table[iter,5]&
                                              node ==-1))
      alpha_dormant=prior_param[2]+length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                              parent[,2]==parent_evidence_table[iter,2]&
                                              parent[,3]==parent_evidence_table[iter,3]&
                                              parent[,4]==parent_evidence_table[iter,4]&
                                              parent[,5]==parent_evidence_table[iter,5]&
                                              node ==0))
      
      alpha_active=prior_param[3] + length(which(parent[,1]==parent_evidence_table[iter,1]& 
                                              parent[,2]==parent_evidence_table[iter,2]&
                                              parent[,3]==parent_evidence_table[iter,3]&
                                              parent[,4]==parent_evidence_table[iter,4]&
                                              parent[,5]==parent_evidence_table[iter,5]&
                                              node ==0))
      alpha_dirichlet_posterior <- rbind(alpha_dirichlet_posterior,c(alpha_inhibit,alpha_dormant,alpha_active))
      
    }
    alpha_dirichlet_posterior<-alpha_dirichlet_posterior[-1,]
    
  }
  
  return(alpha_dirichlet_posterior)
}