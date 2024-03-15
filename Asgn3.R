
## Install/Load libraries ####

library(glmnet)
library(pROC)
library(ggplot2)

## Exploratory data analysis ####
#Load data, initially converted from .xlsx to in .csv format using google sheets 
covid_data <- read.csv("Immunologic profiles of patients with COVID-19.xlsx - Immunologic profiles of patient.csv")

#Initial exploration of dataset
str(covid_data)

summary(covid_data)
#of note, IL.10 has a mean of 330.63, but a maximum of 15237.08 (2 orders of magnitude difference, and patient had severe covid)

#check if coefficients are standardized, they do not appear to be standardized by thankfully cv.glmnet() does it automatically 
apply(covid_data[, c(2, 5:32)],2,mean)
apply(covid_data[, c(2, 5:32)],2,sd)

## Model development ####

#set seed for reproducibility
set.seed(1717)

#Convert 'Severity' to a factor variable
#resp column is now column number 33, factor with 2 levels (1,2) wherre mild is 1 and severe is 2.  
covid_data$resp <- as.factor(covid_data$Severirty)

#split dataset into training and test sets 
#training index samples from 1 to the total number of rows in covid_data dataframe, taking only 75% number of rows, without replacement
train_index <- sample(1:nrow(covid_data), 0.75*nrow(covid_data), replace = FALSE)

#covid_data dataset is indexed using samples in train_index, allowing for random selection of 75% of dataset rows for training set 
train_data <- covid_data[train_index, ]

#covid_data is indexed for any rows that are not present in train_data, allowing for random selection of 25% of covid_data dataset rows for test set
test_data <- covid_data[-train_index, ]

#predictor variables passed onto x_train variable (age and blood serum features) 
x_train <- train_data[ ,c(2, 5:32)]

#response variable is 'Severity" as a factor variable
y_train <- train_data[ ,33]

#predictor variables passed onto x_test variable (age and blood serum features) 
x_test <- test_data[ ,c(2, 5:32)]

#response variable is 'Severity" as a factor variable
y_test <- test_data[ ,33]

#For model training, use seq(0, 1, .01) as the search grid for alpha
alphas <- seq(0,1,.01)
str(alphas)

#generate empty model to be populated with information
elastic_model <- matrix(NA,ncol=3,nrow=length(alphas))

#setting fold for replicabilty
elastic_net_foldid <- cv.glmnet(x = as.matrix(x_train),
                                y = y_train,
                                nfolds = 10,
                                family = "binomial",
                                type.measure = "deviance",
                                keep = TRUE
                                )

#foldid of initial elastic net regression set to foldid variable
foldid <- elastic_net_foldid$foldid

#view different folds
foldid

#see if the folds have been properly interspersed
table(foldid)


#for loop to search for ideal alpha value 
for (i in 1:length(alphas)){
  
  #for each alpha from 0 to 1, run elastic net regression 
  elastic_fit <- cv.glmnet(x = as.matrix(x_train), 
                           y = y_train,
                           alpha=alphas[i], 
                           foldid = foldid, 
                           family = "binomial", 
                           type.measure = "deviance")
  #for each elastic net regression determine the minimum cvm, the mean cross-validated error and assign to 'indx' 
  indx <- which.min(elastic_fit$cvm)
  #append to respective alpha row, minimum mean cross-validation error, and upper/lower curve (+/- standard error of cvm)  
  elastic_model[i,] <- c(elastic_fit$cvm[indx], elastic_fit$cvlo[indx],elastic_fit$cvup[indx])
}

#name elastic_model dataframe columns their respective identities
colnames(elastic_model) <- c("cvmin","cvlo","cvup")

#visualize each characteristic of a model given their alpha values 
par(mfrow=c(1,1))
plot(alphas,elastic_model[,1],pch=16,main="cvmin generated from elastic net regression using alpha values 0 to 1", ylab = "cvmin")
plot(alphas,elastic_model[,2],pch=16,main="cvlo generated from elastic net regression using alpha values 0 to 1", ylab = "cvmin - standard error of cvm")
plot(alphas,elastic_model[,3],pch=16,main="cvup generated from elastic net regression using alpha values 0 to 1", ylab = "cvmin + standard error of cvm")


### The best model? It appears to be alpha = 0.21
optimal_alpha <- alphas[which.min(elastic_model[,1])]

#Now we need to determine the optimal lambda value (lambda.min) using the same elastic net regression, while utilizing the previously determined 'best model' alpha value  
set.seed(1717)
cv <- cv.glmnet(x = as.matrix(x_train), 
                y = y_train,
                alpha = optimal_alpha,
                foldid = foldid,
                family = "binomial",
                type.measure = "deviance")

#return the value of the lowest lambda generated 
cv$lambda.min

#view plot, Binomial deviance vs log(lambda), the left dotted line is lambda.min, the right dotted line is lambda.1se
plot(cv)

## fitted values and test set predictions ####

#predict response values using model generated using elastic net regression, where newx is a matrix of training set predictor values, using the minimum lambda value found in model.   
prds.train <- predict(cv,newx = as.matrix(x_train), type = "response", s=cv$lambda.min)[,1]

#predict response values using model generated using elastc net regression, where newx is a matrix of test set predictor values, using the minimum lambda value found in model. 
prds.test <- predict(cv,newx = as.matrix(x_test), type = "response", s=cv$lambda.min)[,1]

## Prediction accuracy

#working with training set, builds a ROC curve
dev.train <- roc(y_train, prds.train)

#area under the curve 0.9143
dev.train

#Visualize the ROC curve for training dataset 
#each point on the curve is associated with threshold value
par(mfrow=c(1,1))
plot(dev.train)

#working with training set, builds a ROC curve
dev.test <- roc(y_test, prds.test)

#area under the curve is 0.8778, is not able not as accurately predict covid severity outcomes in unseen data compared to training set but still good 
dev.test

#view ROC curve for test dataset
#each point on the curve is associated with threshold value
par(mfrow=c(1,1))
plot(dev.test)


max(cv$cvm)
which.max(cv$cvm)

## Sensitivity, specificity, and classification cutoff

#check contingency table of response variable in training set 
table(y_train) 

#view table of predicted response, any value above 0.5 is considered 1, any value below 0.5 is considered 0
table(as.numeric(prds.train>.5))

#generating the confusion matrix
conf.mat1 <- table(y=y_train,yhat=as.numeric(prds.train>.5))

#view confusion matrix
conf.mat1

#use function to compute sensitivity and specifity 
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

#threshold is too skewed towards sensitivity
sn.sp(conf.mat1)

#manually determining an ideal threshold value
sn.sp(table(y=y_train,yhat=as.numeric(prds.train>.55)))

sn.sp(table(y=y_train,yhat=as.numeric(prds.train>.70)))

sn.sp(table(y=y_train,yhat=as.numeric(prds.train>.65)))

sn.sp(table(y=y_train,yhat=as.numeric(prds.train>.60)))

sn.sp(table(y=y_train,yhat=as.numeric(prds.train>.63)))

sn.sp(table(y=y_train,yhat=as.numeric(prds.train>.635)))

## Determining the threshold rigourously #### 
#using training set, generate a matrix of sensitivity and specificity 
snsp.train <- cbind(dev.train$sensitivities,dev.train$specificities)

#view head of matrix
head(snsp.train)

#uses min-max approach to try to maximize both sensitivity and specificity
indx <- which.max(apply(snsp.train,1,min))

#indicates which row contains maximized sensitivity and specificity
indx

#index row with maximized sensitivity and specificity 
snsp.train[indx,]

#this is the t value, if phat is greater than 0.6304922 patient has severe covid
cutoff <- dev.train$thresholds[indx]
cutoff

#the intersection of blue dotted lines, is the threshold that is closest to (0,1) the ideal that threshold goes to.
plot(dev.train, main = "Sensitivity vs Specificity plot using training set")
abline(h=snsp.train[indx,1],v=snsp.train[indx,2], col='blue', lty=2)

#using training set, generate a matrix of sensitivity and specificity
snsp.test <- cbind(dev.test$sensitivities, dev.test$specificities)

#uses min-max approach to try to maximize both sensitivity and specificity
indx2 <- which.max(apply(snsp.test,1,min))  ### min-max approach!

#indicates which row contains maximized sensitivity and specificity
indx2

#index row with maximized sensitivity and specificity 
snsp.test[indx2,]

#the intersection of blue dotted lines, is the threshold that is closest to (0,1) the ideal that threshold goes to
plot(dev.test, main = "Sensitivity vs Specificity plot using test set")
abline(h=snsp.test[indx2,1],v=snsp.test[indx2,2], col='blue', lty=2)

#this is the threshold, t value, if phat is greater than 0.6429771 patient has severe covid
cutoff2 <- dev.test$thresholds[indx2]
cutoff2

#Indicating which cytokines are most effective at predicting disease progression in patients ####

hist(prds.train)

#retrieve the coefficients of the best model (according to lambda.min)
coef.min <- coef(cv , s = cv$lambda.min)[,1]

#transfer to dataframe
coefficients <- data.frame(coef.min)

#remove intercept, retaining rownames
coefficients <- coefficients[-1, ,drop = FALSE]

#determine the standard deviation of the predictive features
sd <- apply(x_train,2, sd)

#transfer to dataframe
sd <- as.data.frame(sd)

#generate standardized coefficients by multiplying with respective standard deviation values
standardized <- coefficients * sd

#return only absolute values
standardized <- abs(standardized)

#remove any rows/covarites that have been reduced to 0.0 
standardized <- standardized[standardized != 0, , drop = FALSE]

#order rows according to standardized coefficient, in descending manner
standardized <- standardized[ order(-standardized$coef.min), , drop = FALSE]

## Determination of how AGE correlates with Covid 

#summary of AGE feature from whole dataset
summary(covid_data$AGE)

#scatter plot of how AGE correlates with covid severity 
p <- ggplot(data = covid_data, aes(x = AGE, y = resp, color = resp)) +
  geom_point(size = 2) + 
  scale_color_manual(values = c("blue", "red"), labels = c("Mild", "Severe"), name = "Severity") +
  labs(x = "X axis", y = "Y axis", title = "Scatter plot with legend")

# Show the plot
print(p)


covid_data$resp
class(covid_data$resp)
#view the possible skew of age feature
hist(covid_data$AGE)

correlation <- cor.test(as.numeric(covid_data$AGE), as.numeric(covid_data$resp), method = "pearson")

print(correlation)

boxplot(covid_data$AGE ~ covid_data$resp, xlab = "Severity", ylab = "Age",
        main = "Age Distribution by Severity")

