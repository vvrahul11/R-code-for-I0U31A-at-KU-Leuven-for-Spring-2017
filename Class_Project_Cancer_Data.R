##  Script runs post being in the directory with the file
load("~/Desktop/ALL_DESKTOP_MATERIAL/ALL_KU_LEUVEN_COURSE_MATERIALS/STATISTICAL METHODS IN BIOINFORMATICS/ASSIGNMENT/VIJVER(1).Rdata")
bcData <- data
ncol(bcData[,-1]) ## how many genes/contigs exlucding first column since it is sample identification
nrow(bcData) ## how many samples
DM <- subset(bcData,meta=="DM") ## all the DM samples 
NODM <- subset(bcData,meta=="NODM") ## all the NODM samples
nrow(DM) ## count of DM samples which was 78
nrow(NODM) ## count of NODM samples which was 110
bc.Data.No.Sample.Type <- bcData[,-1]
bc.Data.Subset <- bc.Data.No.Sample.Type[,1:40]
cor.Between.Genes <- cor(bc.Data.No.Sample.Type, use = "pairwise.complete.obs")
for(i in 1:ncol(cor.Between.Genes)){
  for (j in 1: nrow(cor.Between.Genes)){
    if(cor.Between.Genes[i,j] <= 0.6){
      cor.Between.Genes[i,j] <- NA
    }
  }
 
}
## In genomes, genes regulate each other. Hence high levels of correlation of expressiom
## between genes is expected. In fact in the correlation matrix there were correlation levels 
## greater than 0.7 suggesting several genes in this data had correlated expression
## Thus we cannot include all the genes for our model because it will inflate the 
## coefficients for certain feature genes in our model giving inflated t test statistic which will not be
## trustable. As an example AL157502 has a correlation of 0.84 with NM_003617 which is a strong positive correlation.
## Similarly between Contig

set.seed(1)
bc.Data.Subset.Sample.Type <- bcData[,1:20]
Sample_Id <- bc.Data.Subset.Sample.Type[,1] 
train = sample(188,94) ## total 188 samples with 94 to train and 94 to predict
train
trainingData <- bc.Data.Subset.Sample.Type[train,]
testingData <- bc.Data.Subset.Sample.Type[-train,]
library(plyr)
Contig29982_RC <- bc.Data.Subset.Sample.Type[,c(1,2)]
Contig29982_RC.test.data <- bc.Data.Subset.Sample.Type[c(1,2)][-train,]
Contig29982_RC.training.data <- bc.Data.Subset.Sample.Type[c(1,2)][train,]
boxplot(Contig29982_RC~ meta , data = Contig29982_RC )
Contig29982_RC.model <- glm(meta ~ Contig29982_RC,data=Contig29982_RC.training.data,family=binomial)
summary(Contig29982_RC.model)

Contig42854.test.data <- bc.Data.Subset.Sample.Type[c(1,4)][-train,]
Contig42854.training.data <- bc.Data.Subset.Sample.Type[c(1,4)][train,]
Contig42854 <- bc.Data.Subset.Sample.Type[,c(1,4)]
boxplot(Contig42854~ meta , data = Contig42854,xlab="Patient category", ylab="Contig 42854 expression" ) ## higher contig 42854 in no dm
## boxplot suggests that COntig42854 expresison is higher in NoDM 
Contig42854.model <- glm(meta ~ Contig42854,data=Contig42854.training.data,family=binomial)
## glm model suggests that Contig42854 is a significant predicitor of DM vs NoDM outcome
summary(Contig42854.model)
Contig42854.probabilities = predict(Contig42854.model,Contig42854.test.data,type="response")
contrasts(Contig42854$meta)
glm.pred <- rep("DM", 188)
glm.pred[Contig42854.probabilities > 0.5 ] <- "NoDM"
table(glm.pred,Contig42854$meta)
## We see that that the model had very high values in the 2*2 table for cells NoDM/NODM and DM/DM 
## which should not be the case for an ideal classifier
## Hence single gene prediction does not work as independent expressions are poor determiners
## of cancer outcome
## RIDGE REGRESSION (alpha=0) STill to answer how many genes are used

library(glmnet)

coef(cv.ridge,s=cv.ridge$lambda.min)
##LASSO REGRESSION (alpha=1) STill to answer how many genes are used
## the datasets are wrong, use test and train
require(glmnet)
whole.bcData.MM <- model.matrix(meta~.,data=bcData)[,-1]
whole.bcData.label <- bcData$meta
lambda.Used = 10^seq(10,-2,length=100) 
set.seed(10)
training.Sample <- sample(1:nrow(whole.bcData.MM),nrow(whole.bcData.MM)/2)
## LASSO AND RIDGE WITH NO CROSS VALIDATION
testing.Sample <- (-training.Sample)
training.Sample.label <- bcData$meta[training.Sample]

lasso.no.cv <- glmnet(whole.bcData.MM[training.Sample,],bcData$meta[training.Sample],family=c("binomial"),alpha=1,lambda = lambda.Used)
ridge.no.cv <- glmnet(whole.bcData.MM[training.Sample,],bcData$meta[training.Sample],family=c("binomial"),alpha=0,lambda = lambda.Used)
plot(ridge.no.cv)
plot(lasso.no.cv)

lasso.no.cv.pred <- predict(lasso.no.cv,s=4,newx=whole.bcData.MM[testing.Sample,],type='response')
predict(lasso.no.cv,type="coefficients",s=lambda.Used)
## LASSO AND RIDGE WITH CROSS VALIDATION
lasso.cv.for.best.lambda.10.fold <- cv.glmnet(whole.bcData.MM[training.Sample,],bcData$meta[training.Sample],family=c("binomial"),alpha=1,lambda = lambda.Used,nfolds = 10)
lasso.cv.for.best.lambda.10.fold$nzero
ridge.cv.for.best.lambda.10.fold <- cv.glmnet(whole.bcData.MM[training.Sample,],bcData$meta[training.Sample],family=c("binomial"),alpha=0,lambda = lambda.Used,nfolds=10)
ridge.cv.for.best.lambda.10.fold$nzero
ridge.cv.for.best.lambda.10.fold$cvm
lasso.cv.for.best.lambda.5.fold <- cv.glmnet(whole.bcData.MM[training.Sample,],bcData$meta[training.Sample],family=c("binomial"),alpha=1,lambda = lambda.Used,nfolds = 5)
lasso.cv.for.best.lambda.5.fold$nzero
ridge.cv.for.best.lambda.5.fold <- cv.glmnet(whole.bcData.MM[training.Sample,],bcData$meta[training.Sample],family=c("binomial"),alpha=0,lambda = lambda.Used,nfolds=5)
ridge.cv.for.best.lambda.5.fold$nzero
ridge.cv.for.best.lambda.5.fold$cvm
lasso.with.10.fold.cv.lambda.prediction <- predict(lasso.no.cv,s=lasso.cv.for.best.lambda.10.fold$lambda.min,newx=whole.bcData.MM[testing.Sample,])
lasso.with.5.fold.cv.lambda.prediction <- predict(lasso.no.cv,s=lasso.cv.for.best.lambda.5.fold$lambda.min,newx=whole.bcData.MM[testing.Sample,])
ridge.with.10.fold.cv.lamda.prediction <- predict(ridge.no.cv,s=ridge.cv.for.best.lambda.10.fold$lambda.min,newx=whole.bcData.MM[testing.Sample,])
ridge.with.5.fold.cv.lamda.prediction <- predict(ridge.no.cv,s=ridge.cv.for.best.lambda.5.fold$lambda.min,newx=whole.bcData.MM[testing.Sample,])

## Trying PCA again
pcr.fit <- pcr(meta~.,data=bcData,scale=TRUE,validation="CV")
summary(pcr.fit)
validationplot(pcr.fit,val.type="RMSEP")
plot(pcr.fit,"loadings",comps=1:2,legendpos="topright",labels=1:2)
abline(h=0)
pcr.fit <- pcr(meta~.,data=bcData,subset=training.Sample,scale=TRUE,validation="CV")
validationplot(pcr.fit,val.type = "MSEP")
pcr.pred <- predict(pcr.fit,whole.bcData.MM[testing.Sample,],ncomp=10)
testingData <- bcData[testing.Sample,]
testing.Table.With.Pred <- cbind(pcr.pred,testingData)
testedPredictions <- ifelse(testing.Table.With.Pred$model_pred >= 1.5 , 2,1 )
testing.Table.With.Pred.Final <- cbind(testedPredictions,testing.Table.With.Pred)
table(testing.Table.With.Pred.Final$testedPredictions,testing.Table.With.Pred.Final$meta)
## final model predicted 21/39 correctly class 1 and 45/55 of class 2