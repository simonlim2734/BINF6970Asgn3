## Problem 1 ####
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
par(mar = c(3,2,5,1))
plot(cv, main =  "Binomial deviance vs log(lambda), where nfolds = 10")

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

#Determining accuracy using confusion matrix 
accuracy <- sum(diag(conf.mat1))/sum(conf.mat1)

#check contingency table of response variable in training set 
table(y_test) 

#view table of predicted response, any value above 0.5 is considered 1, any value below 0.5 is considered 0
table(as.numeric(prds.test>.5))

#generating the confusion matrix
conf.mat1_test <- table(y=y_test,yhat=as.numeric(prds.test>.5))

#view confusion matrix
conf.mat1_test

#Determining accuracy using confusion matrix 
accuracy_test <- sum(diag(conf.mat1_test))/sum(conf.mat1_test)

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
plot(dev.train, main = "Sensitivity vs Specificity plot using training set, where nfolds = 10")
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
plot(dev.test, main = "Sensitivity vs Specificity plot using test set, where nfolds = 10")
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

#determine the correlation between age and covid severity using cor.test()
#Pearson's product moment correlation used, as the age is a continuous variable
#and severity resp is discrete (0 is mild and 1 is severe)

correlation <- cor.test(as.numeric(covid_data$AGE), (as.numeric(covid_data$resp) - 1))

print(correlation)

#boxplot, where the x axis is the severity of covid, and the y axis is the age of patients 
boxplot(covid_data$AGE ~ covid_data$resp, xlab = "Severity", ylab = "Age",
        main = "Age Distribution by Severity")

## Using nfolds = 20 ####

set.seed(1717)

#generate empty model to be populated with information
elastic_model_20 <- matrix(NA,ncol=3,nrow=length(alphas))

#setting fold for replicabilty
elastic_net_foldid_20 <- cv.glmnet(x = as.matrix(x_train),
                                y = y_train,
                                nfolds = 20,
                                family = "binomial",
                                type.measure = "deviance",
                                keep = TRUE
)

#foldid of initial elastic net regression set to foldid variable
foldid_20 <- elastic_net_foldid_20$foldid

#view different folds
foldid_20

#see if the folds have been properly interspersed
table(foldid_20)


#for loop to search for ideal alpha value 
for (i in 1:length(alphas)){
  
  #for each alpha from 0 to 1, run elastic net regression 
  elastic_fit_20 <- cv.glmnet(x = as.matrix(x_train), 
                           y = y_train,
                           alpha=alphas[i], 
                           foldid = foldid, 
                           family = "binomial", 
                           type.measure = "deviance")
  #for each elastic net regression determine the minimum cvm, the mean cross-validated error and assign to 'indx' 
  indx_20 <- which.min(elastic_fit$cvm)
  #append to respective alpha row, minimum mean cross-validation error, and upper/lower curve (+/- standard error of cvm)  
  elastic_model_20[i,] <- c(elastic_fit_20$cvm[indx_20], elastic_fit_20$cvlo[indx_20],elastic_fit_20$cvup[indx_20])
}

#name elastic_model dataframe columns their respective identities
colnames(elastic_model_20) <- c("cvmin","cvlo","cvup")

#visualize each characteristic of a model given their alpha values 
par(mfrow=c(1,1))
plot(alphas,elastic_model_20[,1],pch=16,main="cvmin generated from elastic net regression using alpha values 0 to 1 where nfolds = 20", ylab = "cvmin")
plot(alphas,elastic_model_20[,2],pch=16,main="cvlo generated from elastic net regression using alpha values 0 to 1 where nfolds = 20", ylab = "cvmin - standard error of cvm")
plot(alphas,elastic_model_20[,3],pch=16,main="cvup generated from elastic net regression using alpha values 0 to 1 where nfolds = 20", ylab = "cvmin + standard error of cvm")


### The best model? It appears to be alpha = 1 
optimal_alpha_20 <- alphas[which.min(elastic_model_20[,1])]

#Now we need to determine the optimal lambda value (lambda.min) using the same elastic net regression, while utilizing the previously determined 'best model' alpha value  
set.seed(1717)
cv_20 <- cv.glmnet(x = as.matrix(x_train), 
                y = y_train,
                alpha = optimal_alpha_20,
                foldid = foldid_20,
                family = "binomial",
                type.measure = "deviance")

#return the value of the lowest lambda generated 
cv_20$lambda.min

#view plot, Binomial deviance vs log(lambda), the left dotted line is lambda.min, the right dotted line is lambda.1se
par(mar = c(3,2,5,1))
plot(cv_20, main = "Binomial deviance vs log(lambda), where nfolds = 20")

## fitted values and test set predictions ####

#predict response values using model generated using elastic net regression, where newx is a matrix of training set predictor values, using the minimum lambda value found in model.   
prds.train_20 <- predict(cv_20,newx = as.matrix(x_train), type = "response", s=cv_20$lambda.min)[,1]

#predict response values using model generated using elastc net regression, where newx is a matrix of test set predictor values, using the minimum lambda value found in model. 
prds.test_20 <- predict(cv_20,newx = as.matrix(x_test), type = "response", s=cv_20$lambda.min)[,1]

## Prediction accuracy

#working with training set, builds a ROC curve
dev.train_20 <- roc(y_train, prds.train_20)

#area under the curve 0.9143
dev.train_20

#Visualize the ROC curve for training dataset 
#each point on the curve is associated with threshold value
par(mfrow=c(1,1))
plot(dev.train_20)

#working with training set, builds a ROC curve
dev.test_20 <- roc(y_test, prds.test_20)

#area under the curve is 0.8778, is not able not as accurately predict covid severity outcomes in unseen data compared to training set but still good 
dev.test_20

#view ROC curve for test dataset
#each point on the curve is associated with threshold value
par(mfrow=c(1,1))
plot(dev.test_20)


max(cv_20$cvm)
which.max(cv_20$cvm)

## Sensitivity, specificity, and classification cutoff

#check contingency table of response variable in training set 
table(y_train) 

#view table of predicted response, any value above 0.5 is considered 1, any value below 0.5 is considered 0
table(as.numeric(prds.train_20>.5))

#generating the confusion matrix
conf.mat1_20 <- table(y=y_train,yhat=as.numeric(prds.train_20>.5))

#view confusion matrix
conf.mat1_20

#Determining accuracy using confusion matrix 
accuracy_20 <- sum(diag(conf.mat1_20)/sum(conf.mat1_20))

#check contingency table of response variable in test set 
table(y_test) 

#view table of predicted response, any value above 0.5 is considered 1, any value below 0.5 is considered 0
table(as.numeric(prds.test_20>.5))

#generating the confusion matrix
conf.mat1_20_test <- table(y=y_test,yhat=as.numeric(prds.test_20>.5))

#view confusion matrix
conf.mat1_20_test

#Determining accuracy using confusion matrix 
accuracy_20_test <- sum(diag(conf.mat1_20_test)/sum(conf.mat1_20_test))

#use function to compute sensitivity and specifity 
sn.sp <- function(mat){
  sn <- mat[2,2]/sum(mat[2,])
  sp <- mat[1,1]/sum(mat[1,])
  return(unlist(list(sensitivity=sn, specificity=sp)))
}

#threshold is too skewed towards sensitivity
sn.sp(conf.mat1_20)

## Determining the threshold rigourously #### 
#using training set, generate a matrix of sensitivity and specificity 
snsp.train_20 <- cbind(dev.train_20$sensitivities,dev.train_20$specificities)

#view head of matrix
head(snsp.train_20)

#uses min-max approach to try to maximize both sensitivity and specificity
indx_20 <- which.max(apply(snsp.train_20,1,min))

#indicates which row contains maximized sensitivity and specificity
indx_20

#index row with maximized sensitivity and specificity 
snsp.train_20[indx_20,]

#this is the t value, if phat is greater than 0.599932 patient has severe covid
cutoff_20 <- dev.train_20$thresholds[indx_20]
cutoff_20

#the intersection of blue dotted lines, is the threshold that is closest to (0,1) the ideal that threshold goes to.
plot(dev.train_20, main = "Sensitivity vs Specificity plot using training set where nfolds = 20")
abline(h=snsp.train_20[indx_20,1],v=snsp.train_20[indx_20,2], col='blue', lty=2)

#using training set, generate a matrix of sensitivity and specificity
snsp.test_20 <- cbind(dev.test_20$sensitivities, dev.test_20$specificities)

#uses min-max approach to try to maximize both sensitivity and specificity
indx2_20 <- which.max(apply(snsp.test_20,1,min))  ### min-max approach!

#indicates which row contains maximized sensitivity and specificity
indx2_20

#index row with maximized sensitivity and specificity 
snsp.test_20[indx2_20,]

#the intersection of blue dotted lines, is the threshold that is closest to (0,1) the ideal that threshold goes to
plot(dev.test_20, main = "Sensitivity vs Specificity plot using test set where nfolds = 20")
abline(h=snsp.test_20[indx2_20,1],v=snsp.test_20[indx2_20,2], col='blue', lty=2)

#this is the threshold, t value, if phat is greater than 0.6429771 patient has severe covid
cutoff2_20 <- dev.test$thresholds[indx2]
cutoff2_20

#Indicating which cytokines are most effective at predicting disease progression in patients ####

hist(prds.train_20)

#retrieve the coefficients of the best model (according to lambda.min)
coef.min_20 <- coef(cv_20 , s = cv_20$lambda.min)[,1]

#transfer to dataframe
coefficients_20 <- data.frame(coef.min_20)

#remove intercept, retaining rownames
coefficients_20 <- coefficients_20[-1, ,drop = FALSE]

#determine the standard deviation of the predictive features
sd <- apply(x_train,2, sd)

#transfer to dataframe
sd <- as.data.frame(sd)

#generate standardized coefficients by multiplying with respective standard deviation values
standardized_20 <- coefficients_20 * sd

#return only absolute values
standardized_20 <- abs(standardized_20)

#remove any rows/covarites that have been reduced to 0.0 
standardized_20 <- standardized_20[standardized_20 != 0, , drop = FALSE]

#order rows according to standardized coefficient, in descending manner
standardized_20 <- standardized_20[ order(-standardized_20$coef.min), , drop = FALSE]

##### Problem 2 ####
##### Adding and manipulating data for PCA ####
load('geneexpression2.rda')

# Parse out the status and cell types from the row names
dat$status <- sub('^(...).*', '\\1', rownames(dat))
dat$status <- ifelse(dat$status == "HEA", "Healthy", "Melanoma")

dat$cell_type <- sub('.*_(.*)_(.*)', '\\1', rownames(dat))
dat$cell_type <- ifelse(dat$cell_type == "EFFE", "Effector", 
                        ifelse(dat$cell_type == "MEM", "Memory", "Naive"))

head(dat)

# Perform PCA on the dataset without the status and cell_type columns
gene_data <- dat[, !(colnames(dat) %in% c("status", "cell_type"))]
pca_result <- prcomp(gene_data, scale. = TRUE)

# Retrieve eigenvalues and eigenvectors
eigenvalues <- pca_result$sdev^2
eigenvectors <- pca_result$rotation

# Number of PCs to consider for the scree plot
num_pcs <- length(eigenvalues)

loading_scores_squared <- eigenvectors^2


##### Plotting Biplot ####

# Calculate the proportion of variance explained by principal components
prop_varex <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Extract cell types from row names
cell_types <- sub('(.*)_(.*)_(.*)', '\\2', rownames(dat))

# Create an empty plot space for the biplot
plot(pca_result$x[,1], pca_result$x[,2], type='n', 
     xlab = paste0("PC1 (", round(prop_varex[1]*100), "%)"), 
     ylab = paste0("PC2 (", round(prop_varex[2]*100), "%)"), 
     main = "PCA Biplot")

# Draw the text for each point
text(pca_result$x[,1], pca_result$x[,2], cell_types, cex=0.7)

# Calculate mean coordinates for each cell type
mean_coords <- aggregate(pca_result$x[, 1:2], list(cell_types), mean)

# Draw arrows from origin to mean coordinates for each cell type
for (i in 1:nrow(mean_coords)) {
  arrows(0, 0, mean_coords[i, 2], mean_coords[i, 3], length=0.1, col='blue')
  text(mean_coords[i, 2], mean_coords[i, 3], labels = mean_coords[i, 1], pos=2, cex=0.7, col='blue')
}

##### Adding cell status ####

# Define colors and shapes based on status and cell type
colors <- ifelse(dat$status == "Healthy", "blue", "green")
shapes <- ifelse(dat$cell_type == "Effector", 16, # Circle
                 ifelse(dat$cell_type == "Memory", 17, # Triangle
                        18)) # Square for Naive

# Create an empty plot space for the biplot
plot(pca_result$x[,1], pca_result$x[,2], type='n', 
     xlab = paste0("PC1 (", round(prop_varex[1]*100), "%)"), 
     ylab = paste0("PC2 (", round(prop_varex[2]*100), "%)"), 
     main = "PCA Biplot")

# Add points to the plot with the specified colors and shapes
points(pca_result$x[,1], pca_result$x[,2], pch=shapes, col=colors)

# Legend for the plot to interpret the colors and shapes
legend("topright" , inset=c(-0.3, -0.3), legend=c("Healthy", "Melanoma", "Effector", "Memory", "Naive"),
       col=c("blue", "green", "black", "black", "black"), pch=c(16, 16, 16, 17, 18),
       title="Legend")

# Calculate mean coordinates for each combination of cell type and status
mean_coords <- aggregate(pca_result$x[, 1:2], list(dat$status, dat$cell_type), mean)

# Draw arrows from origin to mean coordinates for each combination
for (i in 1:nrow(mean_coords)) {
  arrows(0, 0, mean_coords[i, 3], mean_coords[i, 4], length=0.1, col='red')
  text(mean_coords[i, 3], mean_coords[i, 4], labels=paste(mean_coords[i, 1], mean_coords[i, 2], sep=" - "), pos=2, cex=0.7, col='blue')
}

#### Explaining variance ####

# Calculate cumulative variance
cum_varex <- cumsum(prop_varex)

# Plot the scree plot
par(mfrow=c(1,2))  # Split the plotting area into a 1x2 grid

# First plot (on the left)
barplot(prop_varex, main="Scree Plot", xlab="Principal Component", ylab="Proportion of Variance Explained")

# Second plot (on the right)
plot(cum_varex, xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", type='b', 
     main="Cumulative Variance Plot", col="blue", pch=16)

par(mfrow=c(1,1))  # Reset to default single plotting area

