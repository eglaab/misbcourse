#
# R code for the machine learning analyses in the MISB course "Top-Down Systems Biology"
#
# Required input data in the working folder:
#
# - Pre-processed microarray data from previous day (available as Rdata-files from Moodle)
#

#
# Package installations
# (Please make sure you have a recent version of R > 3.4 - otherwise install the current software as follows:
# - Install the current version of R (3.6) for your operating system from https://ftp.gwdg.de/pub/misc/cran/
# - Install the current version of R-Studio (1.2) from: https://www.rstudio.com/products/rstudio/download/
#

# install R-packages for clustering
if(!require('cluster'))
{
	install.packages('cluster')
	require('cluster')
}

if(!require('mclust'))
{
	install.packages('mclust')
	require('mclust')
}

# install R-packages for classification
if(!require('randomForest'))
{
	install.packages('randomForest')
	require('randomForest')
}

if(!require('e1071'))
{
	install.packages('e1071')
	require('e1071')
}




# set the location of your working directory (note that there are differences between Windows & Mac concerning the use of back slash "\" vs. forward slash "/")

# format for Mac & Linux systems
setwd('/Users/set/your/current/working/directory/here')

# format for Windows
setwd('C:/set/your/current/working/directory/here')



# load datasets from previous days
# (if you haven't saved the data from day 1, you can get the pre-processed data as Rdata-files from Moodle)
#

load(file="zhang_preprocessed.Rdata") # zhangvsn, zhang_outcome_final
load(file="moran_preprocessed.Rdata") # moranvsn, moran_outcome_final



#
# Sample clustering analyses
#

# Function for variance filtering
# variance filter (filters the gene expression matrix X to only retain the genes with highest variance
# Parameter 1: X = input gene expression matrix
# Parameter 2: filtsize = number of genes with highest variance to retain after filtering
var_filter = function(X, filtsize=1000)
{
	filtsize <- min(nrow(X), as.numeric(filtsize))
	variances <- apply(X, 1, var)
	o <- order(variances, decreasing=TRUE)
	Xfilt <- (X[o,])[1:filtsize,]
	return(Xfilt)
}

# filter the expression matrices to only retain the top 2000 genes with the highest variance
zhang_filt = var_filter(zhangvsn, filtsize=2000)
moran_filt = var_filter(moranvsn, filtsize=2000)

# the dimensions are reduced to 2000 genes
dim(zhang_filt)
#[1] 2000   26
dim(moran_filt)
#[1] 2000   39




#
# k-Means clustering
#

# set a seed number for the random number generator to make the results reproducible
set.seed(1234)

# Zhang et al. dataset

# compute Euclidean distance matrix
distmat_zhang = dist(t(zhang_filt), method = "euclidean")

# two-group clustering (k = 2) with 100 random initializations to obtain robust results
kclust2_zhang <- kmeans(t(zhang_filt), 2, nstart = 100)

# show clustering results (two clusters indexed with no. 1 and 2):
kclust2_zhang$cluster

# PCA plot of the clustering results (requires cluster package loaded, see above)
clusplot(distmat_zhang, kclust2_zhang$cluster, diss = TRUE, main="Zhang et al. - PCA clustering plot (ellipses correspond to clusters)")

# three-group clustering (k = 3) with 100 random initializations to obtain robust results
kclust3_zhang <- kmeans(t(zhang_filt), 3, nstart = 100)

# show clustering results (two clusters indexed with no. 1 and 2):
kclust3_zhang$cluster

# PCA plot of the clustering results (requires cluster package loaded, see above)
clusplot(distmat_zhang, kclust3_zhang$cluster, diss = TRUE, main="Zhang et al. - PCA clustering plot (ellipses correspond to clusters)")


# Moran et al. dataset

# compute Euclidean distance matrix
distmat_moran = dist(t(moran_filt), method = "euclidean")

# two-group clustering (k = 2) with 100 random initializations to obtain robust results
kclust2_moran <- kmeans(t(moran_filt), 2, nstart = 100)

# show clustering results (two clusters indexed with no. 1 and 2):
kclust2_moran$cluster

# PCA plot of the clustering results (requires cluster package loaded, see above)
clusplot(distmat_moran, kclust2_moran$cluster, diss = TRUE, main="moran et al. - PCA clustering plot (ellipses correspond to clusters)")

# three-group clustering (k = 3) with 100 random initializations to obtain robust results
kclust3_moran <- kmeans(t(moran_filt), 3, nstart = 100)

# show clustering results (two clusters indexed with no. 1 and 2):
kclust3_moran$cluster

# PCA plot of the clustering results (requires cluster package loaded, see above)
clusplot(distmat_moran, kclust3_moran$cluster, diss = TRUE, main="moran et al. - PCA clustering plot (ellipses correspond to clusters)")



#
# Hierarchical clustering
#

# Compute average linkage hierarchical clustering (Zhang et al.)
hcldat_zhang = hclust(distmat_zhang, method="average")

# Show cluster dendrogram
plot(hcldat_zhang)

# turn hierarchical clustering into "flat" clustering by cutting the tree at different levels:
hcl2_zhang = cutree(hcldat_zhang,2)
hcl3_zhang = cutree(hcldat_zhang,3)

# Compute average linkage hierarchical clustering (Moran et al.)
hcldat_moran = hclust(distmat_moran, method="average")

# Show cluster dendrogram
plot(hcldat_moran)

# turn hierarchical clustering into "flat" clustering by cutting the tree at different levels:
hcl2_moran = cutree(hcldat_moran,2)
hcl3_moran = cutree(hcldat_moran,3)



# Internal cluster validity assessment - Average Silhouette width

# average silhouette width (Zhang et al., k-means, k = 2)
kclust2_zhang_score <- mean((silhouette(kclust2_zhang$cluster, distmat_zhang))[,3])
kclust2_zhang_score

# average silhouette width (Zhang et al., k-means, k = 3)
kclust3_zhang_score <- mean((silhouette(kclust3_zhang$cluster, distmat_zhang))[,3])
kclust3_zhang_score

# average silhouette width (Zhang et al., hclust, k = 2)
hcl2_zhang_score  <- mean((silhouette(hcl2_zhang, distmat_zhang))[,3])
hcl2_zhang_score

# average silhouette width (Zhang et al., hclust, k = 3)
hcl3_zhang_score  <- mean((silhouette(hcl3_zhang, distmat_zhang))[,3])
hcl3_zhang_score

# Plot of silhouette widths for Zhang et al. data
plot(silhouette(kclust3_zhang$cluster, distmat_zhang), main="Silhouette plot of best clustering result")


# average silhouette width (moran et al., k-means, k = 2)
kclust2_moran_score <- mean((silhouette(kclust2_moran$cluster, distmat_moran))[,3])
kclust2_moran_score

# average silhouette width (moran et al., k-means, k = 3)
kclust3_moran_score <- mean((silhouette(kclust3_moran$cluster, distmat_moran))[,3])
kclust3_moran_score

# average silhouette width (moran et al., hclust, k = 2)
hcl2_moran_score  <- mean((silhouette(hcl2_moran, distmat_moran))[,3])
hcl2_moran_score

# average silhouette width (moran et al., hclust, k = 3)
hcl3_moran_score  <- mean((silhouette(hcl3_moran, distmat_moran))[,3])
hcl3_moran_score

# Plot of silhouette widths for moran et al. data
plot(silhouette(kclust2_moran$cluster, distmat_moran), main="Silhouette plot of best clustering result")


# External cluster validity assessment - Adjusted rand index

# Zhang et al. - k-Means, k = 2
adjustedRandIndex(kclust2_zhang$cluster, zhang_outcome_final)

# Zhang et al. - k-Means, k = 3
adjustedRandIndex(kclust3_zhang$cluster, zhang_outcome_final)

# Zhang et al. - k-Means, k = 2
adjustedRandIndex(hcl2_zhang, zhang_outcome_final)


# Moran et al. - k-Means, k = 3
adjustedRandIndex(hcl3_moran, moran_outcome_final)

# Moran et al. - k-Means, k = 2
adjustedRandIndex(kclust2_moran$cluster, moran_outcome_final)

# Moran et al. - k-Means, k = 3
adjustedRandIndex(kclust3_moran$cluster, moran_outcome_final)

# Moran et al. - k-Means, k = 2
adjustedRandIndex(hcl2_moran, moran_outcome_final)

# Moran et al. - k-Means, k = 3
adjustedRandIndex(hcl3_moran, moran_outcome_final)



#
# Sample classification analyses
#

# Evaluation functions

# sensitivity
sensitivity <- function(tp, fn){
  return(tp/(tp+fn))
}

# specificity
specificity <- function(tn, fp){
  return(tn/(tn+fp))
}

# Matthew's correlation coefficient (=MCC)
corcoeff <- function(tp, tn, fp, fn){
	return(  ((tp*tn)-(fp*fn))/(sqrt((tn+fn)*(tn+fp)*(tp+fn)*(tp+fp))))
}

# Huberty's proportional chance criterion - p-value calculation for classification problems
ppc <- function(tp, fp, tn, fn)
{

		total <- tp+fp+fn+tn
		
		c_pro <- ((tp+fn)/total)*((tp+fp)/total) + ((tn+fp)/total) *((tn+fn)/total)
		
		acc <- (tp+tn)/total
		
		cat('\nc_pro:',c_pro,'acc:',acc,'\n')
		
		
		pval <- NULL
		if(c_pro > acc)
		{
		 pval <- 1
		} else if(c_pro != 1)
		{
			z <- (acc-c_pro)/sqrt(c_pro*(1-c_pro)/total)
		
			pval <- pnorm(-abs(z))
			cat('\np-value: ',pval,'\n')
			cat('\np-value (rounded): ',format(pval,digits=2),'\n')
	
		} else {
		  pval <- 1
		}

  return (pval)	
}


# make the outcome variables numeric

zhang_numout = ifelse(zhang_outcome_final=="disease state: Control",0,1)
moran_numout = ifelse(moran_outcome_final=="control",0,1)

set.seed(1234)

# Random Forest model for Zhang et al. data using 250 trees
rfmod_zhang = randomForest(t(zhangvsn), factor(zhang_numout), ntree=250, keep.forest=TRUE)

# show model evluation based on out-of-bag samples
rfmod_zhang

# compute performance statistics
sensitivity(rfmod_zhang$confusion[2,2], rfmod_zhang$confusion[2,1])
specificity(rfmod_zhang$confusion[1,1], rfmod_zhang$confusion[1,2])
corcoeff(rfmod_zhang$confusion[2,2], rfmod_zhang$confusion[1,1], rfmod_zhang$confusion[1,2], rfmod_zhang$confusion[2,1])
ppc(rfmod_zhang$confusion[2,2], rfmod_zhang$confusion[1,2], rfmod_zhang$confusion[1,1], rfmod_zhang$confusion[2,1])

# which variables were most informative for the prediction (multivariate feature selction - see day 1 lecture):
head(rfmod_zhang$importance[order(rfmod_zhang$importance, decreasing=T),])


# Random Forest model for Moran et al. data using 250 trees
rfmod_moran = randomForest(t(moranvsn), factor(moran_numout), ntree=250, keep.forest=TRUE)

# show model evluation based on out-of-bag samples
rfmod_moran

# compute performance statistics
sensitivity(rfmod_moran$confusion[2,2], rfmod_moran$confusion[2,1])
specificity(rfmod_moran$confusion[1,1], rfmod_moran$confusion[1,2])
corcoeff(rfmod_moran$confusion[2,2], rfmod_moran$confusion[1,1], rfmod_moran$confusion[1,2], rfmod_moran$confusion[2,1])
ppc(rfmod_moran$confusion[2,2], rfmod_moran$confusion[1,2], rfmod_moran$confusion[1,1], rfmod_moran$confusion[2,1])


# which variables were most informative for the prediction (multivariate feature selction - see day 1 lecture):
head(rfmod_moran$importance[order(rfmod_moran$importance, decreasing=T),])



# Support vector machine classification

# Run linear SVM - evaluate using leave-one-out cross-validation (Zhang et a.)
svmmod = svm(t(zhangvsn), factor(zhang_numout), kernel="linear", cross=ncol(zhangvsn))

# show the cross-validated accuracy (as percentage)
svmmod$tot.accuracy

# confusion matrix
conf_zhang = table(zhang_numout, ifelse(svmmod$accuracies==100,zhang_numout,1-zhang_numout))
conf_zhang

# compute performance statistics
sensitivity(conf_zhang[2,2], conf_zhang[2,1])
specificity(conf_zhang[1,1], conf_zhang[1,2])
corcoeff(conf_zhang[2,2], conf_zhang[1,1], conf_zhang[1,2], conf_zhang[2,1])
ppc(conf_zhang[2,2], conf_zhang[1,2], conf_zhang[1,1], conf_zhang[2,1])


# Run linear SVM - evaluate using leave-one-out cross-validation (Moran et a.)
svmmod = svm(t(moranvsn), factor(moran_numout), kernel="linear", cross=ncol(moranvsn))

# show the cross-validated accuracy (as percentage)
svmmod$tot.accuracy

# confusion matrix
conf_moran = table(moran_numout, ifelse(svmmod$accuracies==100,moran_numout,1-moran_numout))
conf_moran

sensitivity(conf_moran[2,2], conf_moran[2,1])
specificity(conf_moran[1,1], conf_moran[1,2])
corcoeff(conf_moran[2,2], conf_moran[1,1], conf_moran[1,2], conf_moran[2,1])
ppc(conf_moran[2,2], conf_moran[1,2], conf_moran[1,1], conf_moran[2,1])

#
# Other online resources for machine learning analysis of omics data:
# ArrayMining: www.arraymining.net
# PathVar: www.pathvar.embl.de
#


# For reproducibility: show and save information on all loaded R package versions
sessionInfo()
