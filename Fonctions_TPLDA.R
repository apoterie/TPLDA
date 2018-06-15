#Remark: works for the case when y has only 0 and 1 categorie.

#*******************************************
# LIBRARY USED
#*******************************************
library(penalizedLDA)
library(cvTools)
library(e1071)

#******************************************* 
#  entropy() and gini()
# ==> Entropy and Gini Index
#******************************************* 

entropy <- function(p) {
  if (any(p == 1)) return(0)# 
  -sum(p*log(p,2)) 
}

gini <- function(p) {
  sum(p * (1 - p))
}

# INPUT PARAMETERS
###     p = a vector including the class probabilities

# OUTPUT PARAMETERS
###   the impurity value  


#**************************************************************************************
# Tree_PLDA()
# ==> Elaboration of a TPLDA tree
#**************************************************************************************


Tree_PLDA <-function(data,group,crit = 1,case_min = 3,kfold = 3,penalty="No",grp.importance=TRUE) {
  names(data)<-"Y"
  nb_group<-length(unique(group[!is.na(group)]))
  group_size<-sapply(unique(group[!is.na(group)]),function(x)length(which(group==x)))
  if(grp.importance==TRUE){
    importance<-list()
    cimportance<-list()
    agreement<-list()
  }
  tree <- NULL
  parent<-c()
  depth<-c()
  var<-c()
  ldas<-list()
  n<-c()
  n_case<-c()
  n_noncase<-c()
  yval<-c()
  yhat<-c()
  prob<-c()
  pop<-list()
  splits<-list()
  action<-c()
  lambdas <- list()
  lambda<-c()
  taille_groupe<-c()
  n_zero<-c()
  nums <- list() 
  surrogates<-list()
  improvment<-c()
  parent[1] <- NA
  depth[1] <- 0
  n[1] <- dim(data)[1]
  pop[[1]] <- rownames(data)
  i <- 1 
  while (i <= length(n)) {
    node <- data[intersect(unlist(pop[[i]]), rownames(data)), ]
    yval[i] <- ifelse((length(which(node$Y == "1")) / n[i]) < 0.5, "0", "1")
    prob[i] <- round((length(which(node$Y == "1")) / n[i]), 4)
    n_case[i] <- length(which(node$Y == "1"))
    n_noncase[i] <- length(which(node$Y == "0"))
    
    #a stopping criterion fulfills ?
    if ((n_case[i] >= case_min) & (n_noncase[i] >= case_min)) { 
      lambdas[[i]] <- sapply(unique(group[!is.na(group)]),function(x)unlist(cv_ptlda(node, igroup = which(group == x),nfold = kfold,crit = crit,label = yval[i])$shrunk))
      splits[[i]] <- split_ptlda(node = node,group = group,label = yval[i],lambda=unlist(lambdas[[i]]),penalty=penalty,nb_group=nb_group,group_size=group_size)
      improvment[i]<-max(unlist(splits[[i]][[crit]]))
      if (improvment[i] > 0){
        action[i] <- 1
        if(length(which(unlist(splits[[i]][[crit]])==improvment[i]))>1){
          surrogate<-rep(NA,nb_group)
          ind_var<-sample(which(unlist(splits[[i]][[crit]])==improvment[i]),1,FALSE)
          surrogate[which(improvment[i]==unlist(splits[[i]][[crit]]))]<-1
          surrogate[ind_var]<-0
          surrogates[[i]]<-surrogate
        }else{
          ind_var<-which.max(unlist(splits[[i]][[crit]]))
        }
        var[i] <- unique(group[-1])[ind_var]
        pred<-as.numeric(as.character(unlist(splits[[i]]$pred[,ind_var])))
        if(grp.importance==TRUE){
          importance[[i]]<-unlist(splits[[i]][[crit]])
          agreement[[i]]<-sapply(1:nb_group,function(x)round(classAgreement(table(unlist(splits[[i]]$pred[,x]),unlist(splits[[i]]$pred[,ind_var])),match.names=TRUE)$diag,4))
          cimportance[[i]]<-unlist(importance[[i]])*unlist(agreement[[i]])
        }
        if (dim(table(pred)) == 2){
          parent <- c(parent, i, i)
          depth <- c(depth, depth[i] + 1, depth[i] + 1)
          n <- c(n, length(which(pred == "1")), length(which(pred == "0")))
          yhat[length(n) - 1] <- "1"
          yhat[length(n)] <- "0"
          pop[[length(n) - 1]] <- rownames(node[which(pred == "1"), ])
          pop[[length(n)]] <- rownames(node[which(pred == "0"), ])
          ldas[[i]] <- splits[[i]]$ldas[[ind_var]]
          nums[[i]] <- unlist(splits[[i]]$indices[[ind_var]])
          n_zero[i] <-group_size[ind_var]-unlist(splits[[i]]$n_effec)[ind_var]
          taille_groupe[i]<-group_size[ind_var]
          lambda[i]<-ifelse(cv==TRUE,unlist(lambdas[[i]])[ind_var],0)
        }else{
          action[i] <- -3 
          n_zero[i] <-group_size[ind_var]-unlist(splits[[i]]$n_effec)[ind_var]
          taille_groupe[i]<-group_size[ind_var]
          lambda[i]<-ifelse(cv==TRUE,unlist(lambdas[[i]])[ind_var],0)
        }
      }else{
        action[i] <- -2 
        var[i] <- NA
        n_zero[i] <-NA
        taille_groupe[i]<-NA
        lambda[i]<-NA
        improvment[i]<-NA
      }
    }else{
      action[i] <- -1 
      improvment[i]<-NA
      var[i] <- NA
      n_zero[i] <-NA
      taille_groupe[i]<-NA
      lambda[i]<-NA
    }
    i<-i+1
    pred<-NULL
    ind_var<-NULL
  }
  if(grp.importance==FALSE){
    importance<-NULL
    cimportance<-NULL
    agreement<-NULL
  }
  leave <- ifelse(action < 0, "*", "")
  node<-1:length(depth)
  tree <-as.data.frame(cbind(action,var,lambda=round(lambda,3),taille_groupe,n_zero,depth, parent,n, n_case,n_noncase,yhat,yval,prob,leave,node))
  return(list(tree=tree, ldas=ldas, nums=nums, splits=splits,lambdas=lambdas,surrogates=surrogates,improvment=improvment,importance=importance,cimportance=cimportance,agreement=agreement))
}

# INPUT PARAMETERS
###     data: the training sample. Note that it is a data frame with n lines (for observations) and p columns (for variables). The first column is the output vector named "Y" and with the lable "0" and "1". 
###                            The p-1 others variables are continuous.
###     group: a vector of size p indicating the group of each variable. The first component of the vector is NA, it is for the output variable Y.
###     crit: the used impurity function. Write 1 for Gini, 2 for entropy and 3 for the misclassification rate. The impurity function used by default is Gini index.
###     case_min: the minimum number of observations that must exist in a node in order for a split to be attempted.
###     kfold: the number of folds in the cross-validation (the cross-validation is used to select the tuning paramter lambda in the Penalized linear discriminant analysis)
###     penalty: the penalty function used in the splitting criterion. Write "No" for no penalty (the default value), "Size" for the penalty function pen(d_j)=1/d_j (where d_j is the size of the j-th group of inputs),
###                                                               "Root.size" for the penalty function pen(d_j)=1/sqrt(d_j) and "Log" for the penalty function pen(d_j)=1/max{log(d_j),1}. 
###     group.importance: a boolean. Should importance of predictors be assessed?


# OUTPUT VALUES, a list including
###    tree : data frame with one row for each node in the tree. The row.names of tree contain the (unique) node numbers. It is also contained in the column named "node" in tree.
###           Columns of tree include "action" (1 if split, -1 if no split and a stopping criterion fulfilled,-2 and -3 no split because splitting does not improve the node impurity),
###           "var"  (group number of the splitting groups), "lambda" (the penalty used in plda for the splitting group), "taille_groupe" (size of the splitting group), "n_zero" 
###           (number of zero in the discriminant vector), "depth" (node depth in the tree), "parent" (the ancestors), "n" (node size), "n_case" (number of observations in class "Y=1"
###           in the node), "n_noncase" (number of observations in class "Y=0" in the node), "yhat" (class label in the node assigned by the plda), "yval" (class label in the tree, it 
###           is the majority vote), "prob" (empirical probability that an observation in the node belongs to the class "Y=1"), "leave" (the leaves/terminal nodes are indicated by a star "*").
###    lda: a list including all the penalized lda used to built the tree
###    nums : a list including the indices of the inputs variables used in each plda 
###    splits : a list including the results of the function split_plda() 
###    lambdas : a list containing the lambda vector for each plda
###    surrogates : if there are several "best" groups when coosing the splitting group. the vector gives tie groups.
###    improvment : the vector giving the penalized decrease in node impurity resulting from splitting the node. 
###    importance : a list giving, for each group and each node, the penalized decrease in node impurity resulting from splitting the node by the given group. 
###    agreement  : a list giving for each group and each node, the agreement probability between the split based on the group and the best split 
###    cimportance: a list giving for each group and each node, the product importance*agreement. It is used to measure the importance of each group.


#**************************************************************************************  
# cv_ptlda()  
# ==> Perform cross-validation for penalized linear discriminant analysis
#**************************************************************************************  

cv_ptlda <-function(node,igroup = igroup, nfold = 3,lambdas = c(1e-4,1e-3,0.01,0.1,1,10),crit=1,label) {
  w0 <- which(node$Y == "0")
  w1 <- which(node$Y == "1")
  samples0 <- cvFolds(length(w0),K = nfold, R = 1,type = "random")
  samples1 <- cvFolds(length(w1), K = nfold,R = 1, type = "random")
  cverror <- matrix(rep(NA, length(lambdas) * 3), ncol = 3, nrow = length(lambdas))
  tables<-list()
  for (i in 1:length(lambdas)) {
    lda<-NULL
    pred <- rep(NA, length(node$Y))
    for (k in 1:nfold) {
      train <-node[c(w0[samples0$subset[samples0$which != k]], w1[samples1$subset[samples1$which !=k]]), ]
      test <-node[c(w0[samples0$subset[samples0$which == k]], w1[samples1$subset[samples1$which == k]]), ]
      num <- igroup
      for (l in igroup) {
        if ((var(train[, l]) == 0) |(var(train[train$Y == 1, l]) == 0) | (var(train[train$Y == 0, l]) == 0) | (mean(train[train$Y == 1, l]) == mean(train[train$Y == 0, l]))) { num <- num[-which(num == l)]}
      }
      if (length(num) > 0) {
        lda <- PenalizedLDA(as.matrix(train[, num]),as.numeric(as.character(train$Y)) + 1,type = "standard",lambda = lambdas[i], K = 1)
        pred[c(w0[samples0$subset[samples0$which == k]], w1[samples1$subset[samples1$which == k]])] <- as.numeric(predict(lda, as.matrix(test[, num]))$ypred) - 1
      } 
    }
    if(2%in%pred){pred<-ifelse(pred==2,"1","0")}
    tab <- table(pred, node$Y)
    tables[[i]]<-tab
    label<-as.factor(label)
    levels(label)<-levels(as.factor(node$Y))
    if (dim(tab)[1] == 2) {
      p <- sum(tab[, 2]) / (sum(tab[, 2]) + sum(tab[, 1]))
      prop_n0 <- sum(tab[1, ]) / (sum(tab[, 2]) + sum(tab[, 1]))
      prop_n1 <- sum(tab[2, ]) / (sum(tab[, 2]) + sum(tab[, 1]))
      cverror[i, 1] <- gini(prop.table(table(node$Y))) - (prop_n0 * gini(prop.table(tab[1, ])) + prop_n1 *gini(prop.table(tab[2, ])))
      cverror[i, 2] <- entropy(prop.table(table(node$Y))) - (prop_n0 * entropy(prop.table(tab[1, ])) + prop_n1 *entropy(prop.table(tab[2, ])))
      cverror[i, 3] <-  ((length(which(label != node$Y)) -  tab[1, 2] -  tab[2, 1]) / length(node$Y))
    }
    else{
      cverror[i, 1] <- 0
      cverror[i, 2] <- 0
      cverror[i, 3] <- 0
    }
  }
  if (length(cverror[is.na(cverror[, crit]), crit]) == length(cverror[, crit])) {
    shrunk <- 0
  }else{ 
    shrunk <- mean(lambdas[which(max(cverror[!is.na(cverror[, crit]), crit]) == cverror[!is.na(cverror[, crit]), crit])])
  }
  if(max(cverror[!is.na(cverror[, crit]), crit])==0){
    shrunk<-0
  }
  return(list(shrunk=shrunk, cverror=cverror))
}


# INPUT PARAMETERS
###     node : the node. Note that it is a data frame with n lines (for observations) and p columns (for variables). The first column is the output vector named "Y" and with the lable "0" and "1". 
###                            The p-1 others variables are continuous.
###     igroup : a vector containing the indices of the variables used to split the node
###     nfold : the number of folds in the cross-validation
###     lambdas : a vector of lambda values to be considered.grille de valeurs possible pour lambda. 
###               By default lambda=c(1e-4, 1e-3, 1e-2, .1, 0.5, 1, 5, 10) (= values used in the Penalizedlda package)
###     crit : the used impurity function. Write 1 for Gini, 2 for entropy and 3 for the misclassification rate. The impurity function used by default is Gini index.                                                             "Root.size" for the penalty function pen(d_j)=1/sqrt(d_j) and "Log" for the penalty function pen(d_j)=1/max{log(d_j),1}. 
###     label : the class value of the node (the majority class)


# OUTPUT VALUES, a list including
###    shrunk : the optimal value according to the chosen criterion
###    cverror : a matrix containing for each lambda (in row) the value of each criterion (Gini index, entropy, Misclassification rate)


#**************************************************************************************  
# split_ptlda() 
# ==> Select a split for each group
#************************************************************************************** 


split_ptlda <- function(node, group, label, lambda,penalty="No",nb_group,group_size) {
  Gain_Gini<-rep(0,nb_group)
  Gain_Ent<-rep(0,nb_group)
  Gain_Clas<-rep(0,nb_group)
  pred<-matrix(rep(NA,dim(node)[1]*nb_group),ncol=nb_group)
  ldas <- list()
  n_effec <- rep(0, nb_group)
  indicesList<-list()
  for (j in 1: nb_group) {
    indices <- which(group == unique(group[!is.na(group)])[j])
    for (l in which(group == unique(group[!is.na(group)])[j])) {
      if (var(node[, l]) == 0 |(mean(node[node$Y == 1, l]) == mean(node[node$Y == 0, l])) |(var(node[node$Y == 1, l]) == 0) |(var(node[node$Y == 0, l]) == 0)) {
        indices <- indices[-which(indices == l)]
      }
    }
    indicesList[[j]]<-indices
    if (length(indices) > 0) {
      lda <-PenalizedLDA( as.matrix(node[, indices]),as.numeric(as.character(node$Y)) + 1,type = "standard",lambda = lambda[j],K = 1)
      pred[, j] <- as.numeric(predict(lda, node[, indices])$ypred) - 1
      n_effec[j] <- length(which(lda$discrim != 0))
      ldas[[j]] <- lda
      lda<-NULL
    }
  }
  for (j in 1:nb_group) {
    if (sum(pred[, j] == 2, na.rm = T) > 0) {
      pred[, j] <- ifelse(pred[, j] == 2, "1", "0")
    }
    tab <- table(pred[, j], node$Y)
    label<-as.factor(label)
    levels(label)<-levels(as.factor(node$Y))
    if (dim(tab)[1] == 2) {
      p <- sum(tab[, 2]) / (sum(tab[, 2]) + sum(tab[, 1]))
      prop_n0 <- sum(tab[1, ]) / (sum(tab[, 2]) + sum(tab[, 1]))
      prop_n1 <- sum(tab[2, ]) / (sum(tab[, 2]) + sum(tab[, 1]))
      Gain_Gini[j] <- length(node$Y)*(gini(prop.table(table(node$Y))) - (prop_n0 * gini(prop.table(tab[1, ])) + prop_n1 *gini(prop.table(tab[2, ]))))
      Gain_Ent[j] <- length(node$Y)*(entropy(prop.table(table(node$Y))) - (prop_n0 * entropy(prop.table(tab[1, ])) + prop_n1 *entropy(prop.table(tab[2, ]))))
      Gain_Clas[j] <- length(node$Y)*(((length(which(label != node$Y)) -  tab[1, 2] -  tab[2, 1]) / length(node$Y)))
    }
    else{
      Gain_Gini[j] <- 0
      Gain_Ent[j] <- 0
      Gain_Clas[j] <- 0
    }
    if(penalty=="Size"){
      Gain_Gini[j]<-Gain_Gini[j]/group_size[j]
      Gain_Ent[j]<-Gain_Ent[j]/group_size[j]
      Gain_Clas[j]<-Gain_Clas[j]/group_size[j]
    }
    if(penalty=="Root.size"){
      Gain_Gini[j]<-Gain_Gini[j]/sqrt(group_size[j])
      Gain_Ent[j]<-Gain_Ent[j]/sqrt(group_size[j])
      Gain_Clas[j]<-Gain_Clas[j]/sqrt(group_size[j])
    }
    if(penalty=="Log"){
      if(group_size[j]>1){
        Gain_Gini[j]<-Gain_Gini[j]/(log(group_size[j]))
        Gain_Ent[j]<-Gain_Ent[j]/(log(group_size[j]))
        Gain_Clas[j]<-Gain_Clas[j]/(log(group_size[j])) 
      }
    }
  }
  
  return(list(Gain_Gini=Gain_Gini, Gain_Ent=Gain_Ent, Gain_Clas=Gain_Clas, ldas=ldas, n_effec=n_effec,pred=pred,indicesList=indicesList))
}


# INPUT PARAMETERS
###     node : the node. Note that it is a data frame with n lines (for observations) and p columns (for variables). The first column is the output vector named "Y" and with the lable "0" and "1". 
###                            The p-1 others variables are continuous.
###     group: a vector of size p indicating the group of each variable. The first component of the vector is NA, it is for the output variable Y.
###     label : the class value of the node (the majority class)
###     lambda : a vector with the values of lambda for each group
###     penalty : the penalty function used in the splitting criterion. Write "No" for no penalty (the default value), "Size" for the penalty function pen(d_j)=1/d_j (where d_j is the size of the j-th group of inputs),
###     nb_group : a real indicationg number of groups
###     group_size : a vector indicating the group sizes


# OUTPUT VALUES, a list including
###    Gain_Gini : a vector containing the penalized decrease in Gini index for each group
###    Gain_Entropie : a vector containing the penalized decrease in Entropy for each group
###    Gain_Clas :a vector containing the penalized decrease in misclassification rate for each group. (This criterion is not recommended).
###    ldas : liste des ldas estim'es pour chaque groupe
###    n_effec : a vector containing the number of non-zero coefficient in the discriminant vector for eacg group
###    pred : a matrix containing the predction of each obersavtion in the node for each group
###    indicesList : a list containing for each group the indices of the non-zero coefficient in the plda


# DETAILS          
### the function computes the penalized decrease in node impurity weighted by the node size 
###     i.e. Delta Q(t) = n_t*(Q(t) - (n_t0/n_t)*Q(t0) - (n_t1/n_t)*Q(t1)) 
### ==> It is the global penalized decrease in node impurity (within the constant n = the total number of observation in the training sample)


#**************************************************************************************  
# Prediction_penalizedLDA() 
# ==> Predictions from a TPLDA tree
#************************************************************************************** 


Prediction_penalizedLDA <- function(new, tree, ldas, nums) {
  P <- dim(new)[1]
  pred <- rep(NA, P)
  score<-rep(NA,P)
  noeuds <- rep(NA, P)
  res<-sapply(1:P,predFonctionPLDA,new=new,ldas=ldas,nums=nums,tree=tree)
  
  return(as.data.frame(cbind(
    hat.Y=as.numeric(as.character(res[1,])),
    Y=as.numeric(as.character(new$Y)),
    noeuds=as.numeric(as.character(res[2,])),
    score=as.numeric(as.character(res[3,])))))
}


# INPUT PARAMETERS
###     new : a data frame containing the values at which predictions are required. It must contains the variables that the data frame used to built the tree.
###     ldas : the object "ldas" outputted by the function Tree_PLDA()
###     nums : the object "nums" outputted by the function Tree_PLDA()
###     tree : the object "tree" outputted by the function Tree_PLDA()


# OUTPUT VALUES, a list including
###    hat.Y : the prediction vector
###    Y : the vector containing the true output value
###    noeuds : the vector containing the number of the terminal node 
###    score : the vector containing the score value for each observation


#**************************************************************************************  
# predFonctionPLDA() 
# ==> Function called by Prediction_penalizedLDA()
#************************************************************************************** 


predFonctionPLDA<-function(x, new,ldas, nums,tree){
  pred<-NA
  score<-NA
  noeuds<-NA
  i <- 1
  indice_node<-which(as.numeric(as.character(tree$node))==i)
  while (as.numeric(as.character(tree$action[indice_node])) == 1) {
    if (length(unlist(nums[[indice_node]])) > 0) {
      temp.pred <-as.numeric(predict(ldas[[i]], new[x, unlist(nums[[i]])])$ypred) - 1
    } 
    indice_node<-as.numeric(as.character(which(as.numeric(as.character(tree$parent))==as.numeric(as.character(i)) & as.numeric(as.character(tree$yhat))==temp.pred)))
    i<-as.numeric(as.character(tree$node[indice_node]))
  }
  if (as.numeric(as.character(tree$action[indice_node])) < 1) {#node that contains the observation 
    pred <- as.numeric(as.character(tree$yval[indice_node]))
    score<- as.numeric(as.character(tree$prob[indice_node]))
    noeuds <- i
    i <- 1
  }
  return(c(as.numeric(as.character(pred)),noeuds,as.numeric(as.character(score))))
}

# INPUT PARAMETERS
###     x : the index of an observation in new 
###     new : a data frame containing the values at which predictions are required. It must contains the variables that the data frame used to built the tree.
###     ldas : The object "ldas" outputted by the function Tree_PLDA()
###     nums : The object "nums" outputted by the function Tree_PLDA()
###     tree : The object "tree" outputted by the function Tree_PLDA()

# OUTPUT VALUES, a list including
###    hat.Y : the predicted class for the observation
###    noeuds : the terminal node which contains the observation 
###    score : the score for the observation


#**************************************************************************************  
# tree_seq_PLDA() 
# ==> built a sequence of subtrees with respect to a TPLDA tree and the depth notion
#************************************************************************************** 

tree_seq_PLDA<-function(treeTPLDA){
  d_max<-max(as.numeric(as.character(unlist(treeTPLDA$depth))))
  TPLDA_seq<-list()
  for (i in 0:d_max) {
    tree <- as.data.frame(treePLDA[which(as.numeric(as.character(treePLDA$depth)) <= i), ])
    tree$action <-ifelse(as.numeric(as.character(tree$depth)) == i, -1,as.numeric(as.character(tree$action)))
    tree$var <-ifelse(as.numeric(as.character(tree$depth)) == i,NA,as.numeric(as.character(tree$var)))
    tree$leave <-ifelse(as.numeric(as.character(tree$action)) == -1, "*", "")
    TPLDA_seq[[i+1]] <- tree
  }
  return(TPLDA_seq)
}

# INPUT PARAMETERS
###     tree : The object "tree" outputted by the function Tree_PLDA()

# OUTPUT VALUES, a list including
###    TPLDA_seq : a sequence of trees. Each tree is summed up by a data frame data frame with one row for each node in the tree.

#**************************************************************************************  
# impurete_plda() 
# ==> compute the impurity of a tree. Used to select the final tree 
#************************************************************************************** 


impurete_plda <- function(validation, tree_seq,treePLDA) {
  N_tree <- length(tree_seq)
  N <- dim(validation)[1]
  pred <- list()
  sum_noeuds <- list()
  impurete <- matrix(rep(NA, N_tree * 3), nrow = N_tree, ncol = 3)
  for (k in 1:N_tree) {
    predictions <-as.data.frame(Prediction_penalizedLDA(new=validation,tree=as.data.frame(tree_seq[[k]]),ldas=treePLDA$ldas,nums=treePLDA$nums)[,-4])
    names(predictions) <- c("hat.Y", "Y", "noeuds")
    n <- as.numeric(table(as.numeric(as.character(predictions$noeuds))))
    n1 <- tapply(as.numeric(as.character(predictions$Y)), as.numeric(as.character(predictions$noeuds)),sum)
    p1 <- n1 / n
    p0 <- 1 - p1
    predictions$error.pred <-as.numeric(as.character(apply(predictions[, c("hat.Y", "Y")], 1, function(x)ifelse(x[1] != x[2], "1", "0"))))
    misclass <-tapply(as.numeric(as.character(predictions$error.pred)), as.numeric(as.character(predictions$noeuds)), sum)
    summaryNoeuds <-as.data.frame(cbind(as.numeric(as.character(names(table(as.numeric(as.character(predictions$noeuds)))))), n, n1, p1, p0, misclass))
    names(summaryNoeuds) <-c("nom_noeuds","N","N[Y=1]","P[Y=1]","P[Y=0]","P[hat.Y!=Y]")
    impurete[k, 1] <- sum(apply(summaryNoeuds, 1, function(x)x[2] * gini(x[c(4, 5)]))) / N
    impurete[k, 2] <-sum(apply(summaryNoeuds, 1, function(x)x[2] * entropy(x[c(4, 5)]))) / N
    impurete[k, 3] <-sum(summaryNoeuds[, 6]) / (N)
    impurete<-as.data.frame(impurete)
    names(impurete)<-c("Gini","Information","Misclass")
    pred[[k]]<-predictions
    sum_noeuds[[k]]<-summaryNoeuds
    predictions<-NULL
    summaryNoeuds<-NULL
  }
  list(impurete = impurete,pred = pred,summary_noeuds = sum_noeuds)
}

# INPUT PARAMETERS
###     validation : a data frame containing the values at which predictions are required. It must contains the variables that the data frame used to built the tree.
###     tree-seq :  the list outputted by the function tree_seq_PLDA() 
###     treePLDA : the list outputted by the function Tree_PLDA()

# OUTPUT VALUES, a list including
###    impurity : a matrix containing the impurity values for each impurity criterion and each subtree in the sequence 
###    pred : a list containing the prediction vector for each subtree
###    summary_noeuds :  a list containing for each subtree a data frame which provide some information about the nodes in the subtrees : 
###                      "nom_noeuds" (=node name), "N" (node size), "N[Y=1]" (number of observations in the node with the label "Y=1"), 
###                      P[Y=1] (probability that an observation into the node belongs to the class "Y=1"), P[Y=0] (probability that an 
###                      observation into the node belongs to the class "Y=0"), P[hat.Y!=Y] (node misclassifcation rate)

