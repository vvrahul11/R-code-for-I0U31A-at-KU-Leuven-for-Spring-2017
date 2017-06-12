## MUST INCLUDE THE SA HEART DATA FILE WITH THE ASSIGNMENT

heart.data <- read.table(file="SAheart.data.txt",header=TRUE,sep="\t",row.names=1,check.names = FALSE)

## checking for missing data
check.for.missing <- is.na(heart.data) ## no missing data
summary(heart.data) ## a basic summary of the data, "chd" column to be ignored as they are indicator variab;e
nrow(heart.data) ## total no of observations
heart.data$chd <- as.numeric(heart.data$chd)
heart.data.only.numeric <- heart.data[-c(5,10)]
corr.check <- cor(heart.data.only.numeric,use="pairwise.complete.obs")
for (i in 1:ncol(corr.check)){
  for(j in 1:nrow(corr.check)){
    if(corr.check[i,j] <= 0.7){
      corr.check[i,j] <- NA
    }
  }
} ##obesity and adiposity relatively highly correlated 
heart.data$fam.hist.new <- heart.data$famhist
heart.data$fam.hist.new <- ifelse(heart.data$famhist=="Present",1,0)
family.history.glm <- glm(heart.data$chd~heart.data$fam.hist.new,data=heart.data,family="binomial"(link="logit"))
summary(family.history.glm)
pred.heart.attack <- predict(family.history.glm,type="response")## probability returned
exp(coef(family.history.glm)) ## returns 3.2
## backward selection to select optimal generalised linear model
library(leaps)
##install.packages("bestglm")
library(bestglm)
heart.data.no.pred <- heart.data[,-5]
attach(heart.data.no.pred)
heart.data$chd <- as.factor(heart.data$chd)
heart.data$famhist <- NULL
set.seed(1)
train.data = sample(1:dim(heart.data)[1],dim(heart.data)[1]*0.5)
test.data <- (-train.data)
heart.training.data <- heart.data[train.data,]
library(boot)
library(caret)
## Fit model and to determine its accuracy do cross validation
attach(heart.data)
## now taking the model with significant t test for coefficient for backward selection method
full.model <- glm(formula=chd~sbp+tobacco+ldl+adiposity+typea+obesity+alcohol+age+fam.hist.new,family = binomial,data=heart.data) 
summary(full.model)
AIC(full.model)
### FITTING GLM WITH BACKWARD STEPWISE SELECTION  ###
backwards.2 <- step(full.model,trace=1,direction = "backward") ## choosing model by aic
model.fitted <- formula(backwards.2) ##chd ~ tobacco + ldl + typea + age + fam.hist.new
drop.results <- drop1(full.model,test="Chisq") ## chd~ tobacco+ldl+typea+age+fam.hist.new (these are suggested to be kept)
drop.results1 <- drop1(full.model,test="LRT") ## chd~ tobacco+ldl+typea+age+fam.hist.new (same result as above)
train.data.set <- heart.data[train.data,][c(2,3,5,8,10,9)]
model.to.fit <- glm(chd~ tobacco+ldl+typea+age+fam.hist.new,family=binomial(link="logit"),data=train.data.set)
test.data.set <- heart.data[test.data,][c(2,3,5,8,10,9)]
predictions.glm = predict(model.to.fit,newdata = test.data.set[,-6],type="response") ## predictions as per reduced model
test.data.set$model.predictions <- predictions.glm
test.data.set$new.label <- ifelse(test.data.set$model.predictions >= 0.5,1,0)
confusion.matrix.glm <- table(test.data.set$new.label,test.data.set$chd)
sum(diag(confusion.matrix.glm))/sum(confusion.matrix.glm) 
confusion.matrix.glm ## the confusion matrix for model fit after backward selection
## CROSS VALIDATION TO CHECK HOW FINAL FITTED MODEL PERFORMED
cross.validation.glm <- cv.glm(data=test.data.set,model.to.fit,K=5)$delta
## selected model had adjusted error rate of 18% or accuracy of 82%
## dropped out of model is sbp, adiposity,obesity,alcohol
### TESTING SOME OF THE VARIABLES THAT DROPPED OUT OF MODEL AS INDIVIDUAL PREDICTOR OF OUTCOME ###
sbp.model <- glm(formula=chd~sbp,family=binomial(link="logit"))
summary(sbp.model) ## highly significant
boxplot(sbp~chd,data=heart.data,names=c("no chd","chd"),xlab="Disease State",ylab="Sbp")
adiposity.model <- glm(formula=chd~adiposity,family=binomial(link="logit"))
summary(adiposity.model) ## highly significant
boxplot(adiposity~chd,data=heart.data,names=c("no chd","chd"),xlab="Disease State",ylab="Adiposity")
obesity.model <- glm(formula=chd~obesity,family=binomial(link="logit"))
summary(obesity.model) ## significant
alcohol.model <- glm(formula=chd~alcohol,family=binomial(link="logit"))
summary(alcohol.model) ## not significant
alcohol.model$model
### GAM MODEL THROUGH BACKWARD STEPWISE SELECTION ###
library(gam)
levels <- c("Heart.Disease","no.Heart.Disease")
heart.data$fam.hist.new <- as.numeric(heart.data$fam.hist.new)
full.additive.model <- gam(formula=chd~s(age,4)+s(obesity,4)+s(typea,4)+s(ldl,4)+s(sbp,4)+s(tobacco,4)+s(adiposity,4)+fam.hist.new,family = binomial(link="logit"),data=heart.data,trace=TRUE) 
par(mfrow=c(2,3))
plot.gam(full.additive.model,se=TRUE,col="purple")
summary(full.additive.model)
all.gam.models <- list()
heart.data.rearranged <- heart.data[c(1,2,3,4,5,6,7,8,10,9)]
## all possible models range of variables
library(MASS)
## step.gam for best gam selection by backward stepwise selection
full.additive.model.df3 <- gam(formula=chd~s(age,3)+s(obesity,3)+s(typea,3)+s(ldl,3)+s(sbp,3)+s(tobacco,3)+s(adiposity,3)+fam.hist.new,family = binomial(link="logit"),data=heart.data,trace=TRUE) 
gam.backward.df3 <- step.gam(full.additive.model.df3,scope=list("chd"=~1+s(age,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3)+s(typea,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3)+s(typea,3)+s(ldl,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3)+s(typea,3)+s(ldl,3)+s(sbp,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3)+s(typea,3)+s(ldl,3)+s(sbp,3)+s(tobacco,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3)+s(typea,3)+s(ldl,3)+s(sbp,3)+s(tobacco,3)+s(adiposity,3),
                                                                "chd"=~1+s(age,3)+s(obesity,3)+s(typea,3)+s(ldl,3)+s(sbp,3)+s(tobacco,3)+s(adiposity,3)+fam.hist.new),
                                                                direction = "backward" , trace=TRUE)
## Step:2 chd ~ s(age, 3) + s(obesity, 3) + s(typea, 3) + s(ldl, 3) + s(tobacco,      3) + fam.hist.new ; AIC= 489.2912 
test.data.df3.gam.sig.features <- heart.data.no.pred[test.data,][c("age","obesity","typea","ldl","tobacco","fam.hist.new")]
gam.model.prediction.df3 <- predict(gam.backward.df3,newdata = test.data.df3.gam.sig.features,type="response") 
test.data.df3.gam.sig.features$gamPrediction <- gam.model.prediction.df3
test.data.df3.gam.sig.features$chd <- heart.data$chd[test.data] ## 75% correctly predicted 
test.data.df3.gam.sig.features$predicted.outcome <- ifelse(test.data.df3.gam.sig.features$gamPrediction >= 0.5, 1,0)
table(test.data.df3.gam.sig.features$predicted.outcome,test.data.df3.gam.sig.features$chd) ## 75% acccuracy of prediction


