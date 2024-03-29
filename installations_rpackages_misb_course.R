#
# Installation and loading of R-packages
#

#
# Package installations
# (Please make sure you have a recent version of R > 3.4 - otherwise install the current software as follows:
# - Install the current version of R (3.6) for your operating system from https://ftp.gwdg.de/pub/misc/cran/
# - Install the current version of R-Studio (1.2) from: https://www.rstudio.com/products/rstudio/download/
#

#
# Install this script using the following command:
# source('/PATH-TO-YOUR-WORKING-DIRECTORY/installations_rpackages_misb_course.R')
#


# install R-package for pathway analysis
if(!require('clusterProfiler'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager", suppressUpdates=TRUE)

	BiocManager::install("clusterProfiler")
	require('clusterProfiler')
}

if(!require('GSEABase'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("GSEABase", suppressUpdates=TRUE)
	require('GSEABase')
}

# install GEOquery to access datasets from the GEO database
if(!require('GEOquery'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install("GEOquery", suppressUpdates=TRUE, ask = FALSE)
	require('GEOquery')
}

# install Limma package for statistical analyses
if(!require('limma'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")

	BiocManager::install("limma", suppressUpdates=TRUE)
	require('limma')
}

# install R-packages for text-mining
if(!require('XML'))
{
	install.packages('XML')
	require('XML')
}
if(!require('httr'))
{
	install.packages('httr')
	require('httr')
}

if(!require('easyPubMed'))
{
	install.packages('easyPubMed')
	require('easyPubMed')
}

# load annotation package for gene ID conversion

# for old R version:
# source("http://bioconductor.org/biocLite.R")
# biocLite("hgu133a.db")

if(!require('hgu133a.db'))
{
	if (!requireNamespace("BiocManager", quietly = TRUE))
	    install.packages("BiocManager")
	
	BiocManager::install("hgu133a.db", suppressUpdates=TRUE)
	require('hgu133a.db')
}


# load R-packages for quality control
if(!require('arrayQualityMetrics'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("arrayQualityMetrics", suppressUpdates=TRUE)
  install.packages("gridSVG")
  # install.packages("https://cran.r-project.org/src/contrib/Archive/gridSVG/gridSVG_1.4-3.tar.gz", repos=NULL)
  require('arrayQualityMetrics')
}
if(!require('Biobase'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

  BiocManager::install("Biobase", suppressUpdates=TRUE)
  require('Biobase')
}

# load R-packages for power calculation
if(!require('impute'))
{
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 BiocManager::install("impute", suppressUpdates=TRUE)
 require('impute')
}

if(!require('samr'))
{
  install.packages("samr")
  require('samr')
}
options(error=NULL)

# load R-package for Variance stabilizing normalization
if(!require('vsn'))
{
	if (!requireNamespace("vsn", quietly = TRUE))
	    install.packages("BiocManager")
	BiocManager::install("vsn", suppressUpdates=TRUE)
	require('vsn')
}

# load R-package for meta-analysis
if(!require('metaMA'))
{
	install.packages('metaMA')
	require('metaMA')
}


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



#
# Testing installations
#

if(!require('clusterProfiler'))
{
	print('The clusterProfiler package is not successfully installed!')
}

if(!require('GSEABase'))
{
	print('The GSEABase package is not successfully installed!')
}

if(!require('limma'))
{
	print('The limma package is not successfully installed!')
}

if(!require('XML'))
{
	print('The XML package is not successfully installed!')
}

if(!require('httr'))
{
	print('The httr package is not successfully installed!')
}

if(!require('easyPubMed'))
{
	print('The easyPubMed package is not successfully installed!')
}

if(!require('hgu133a.db'))
{
	print('The hgu133a.db package is not successfully installed!')
}

if(!require('arrayQualityMetrics'))
{
	print('The arrayQualityMetrics package is not successfully installed!')
}

if(!require('Biobase'))
{
	print('The Biobase package is not successfully installed!')
}

if(!require('impute'))
{
	print('The impute package is not successfully installed!')
}

if(!require('samr'))
{
	print('The samr package is not successfully installed!')
}

if(!require('vsn'))
{
	print('The vsn package is not successfully installed!')
}

if(!require('metaMA'))
{
	print('The metaMA package is not successfully installed!')
}

if(!require('cluster'))
{
	print('The cluster package is not successfully installed!')
}

if(!require('mclust'))
{
	print('The mclust package is not successfully installed!')
}

if(!require('randomForest'))
{
	print('The randomForest package is not successfully installed!')
}

if(!require('limma'))
{
	print('The limma package is not successfully installed!')
}

if(!require('GEOquery'))
{
	print('The GEOquery package is not successfully installed!')
}

if(!require('e1071'))
{
	print('The e1071 package is not successfully installed!')
}


