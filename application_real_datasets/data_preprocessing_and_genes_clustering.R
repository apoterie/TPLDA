# =========================================================================================
#
#         Application to microarray data
#
# =========================================================================================



rm(list=ls())
path<-""# path to datasets


library(stringr)
library(fastICA)
#=========================================================================================
#
# data loading
#  
# ========================================================================================
### Leukemia (Golub et al. (1999))
load(paste(path,"Leukemia_data_fusion.RData",sep=""))
for(i in 2:ncol(leuk)){leuk[,i]<-as.numeric(as.character(leuk[,i]))}

### Lymphoma (Shipp et al. (2002))
load(paste(path,"shipp.RData",sep=""))
Ylym<-ifelse(shipp$y=="FL","1","0")#"1"= FL  et "0"=DLBCL
Xlym<-shipp$x
lym<-as.data.frame(cbind(Ylym,Xlym))
names(lym)[1]<-"Y"
for(i in 2:ncol(lym)){lym[,i]<-as.numeric(as.character(lym[,i]))} 
table(Ylym);dim(Xlym)

### Colon (Alon et al. (1999))
load(paste(path,"alon.RData",sep=""))
YCol<-ifelse(alon$y=="t","1","0")#"1"= t  et "0"=n
XCol<-alon$x
col<-as.data.frame(cbind(YCol,XCol))
names(col)[1]<-"Y"
for(i in 2:ncol(col)){col[,i]<-as.numeric(as.character(col[,i]))}  


#=========================================================================================
# Data preprocessing
# (Sewak et al. (2009))
# ========================================================================================

## A completer suivant le jeux de données : leuk, pros ou lun
X<-#the matrix with the genes in column and the samples in lines
Y<-#the response variable


## Data range (a ceiling of 16000 and a floor of 20)
X<-apply(X,2,function(x){ifelse(x<20,20,ifelse(x>16000,16000,x))})


## Keep only the 25% of genes with the grestest variation across th sample
cvar <- function(x){(sd(x)/mean(x))*100}
CV<-apply(X,2,cvar);length(CV) # Results in %
summary(abs(CV))
CV_threshold<-quantile(abs(CV),prob=0.75);print(CV_threshold)
indices<-which(abs(CV)>=CV_threshold);length(indices)
Xc<-X[,indices]

## Index of the selected genes
num_genes<-paste("Gene_",indices,sep="")


#=========================================================================================
# Genes clustering (Lee & Batzoglou (2003))
# ========================================================================================
nb_factors<-# number of groups

ica<-fastICA(Xc,nb_factors)
image.plot(1:nb_factors,1:dim(Xc)[2],ica$A,xlab="Factors",ylab="Genes",zlim=range(ica$A))

threshold_ica<-# percentage of genes retained

## Elaboration of the groups: we retain only the threshold_ica% of genes with the largest loads (in absolute value)
matrix_group<-matrix(ncol=ncol(Xc),nrow=nb_factors)
for(i in 1:nb_factors){
  matrix_group[i,]<-ifelse(ica$A[i,]>=quantile(ica$A[i,],probs=1-(threshold_ica/2)),1,ifelse(ica$A[i,]<quantile(ica$A[i,],probs=(threshold_ica/2)),1,0))
}


## Elaboration of the dataset with the groups
X_clust<-Xc[,which(matrix_group[1,]==1)]
group<-c(NA,rep(1,sum(matrix_group[1,])))
names_X<-as.character(num_genes)[which(matrix_group[1,]==1)]
for(i in 2:nb_factors){
  X_clust<-cbind(X_clust,Xc[,which(matrix_group[i,]==1)])
  group<-c(group,rep(i,sum(matrix_group[i,])))
  names_X<-c(names_X,as.character(num_genes)[which(matrix_group[i,]==1)])
}
table(group)
dataset<-as.data.frame(cbind(Y,X_clust))
names(dataset)[1]<-"Y"

#dataset$Y<-ifelse(dataset$Y=="2","1","0")
names(dataset)[2:dim(dataset)[2]]<-paste("N",names_X,"_G",group[-1],sep="")
for(l in 2:ncol(dataset)){dataset[,l]<-as.numeric(as.character(dataset[,l]))}
dataset[,-1]<-scale(dataset[,-1],scale=T,center=T)


