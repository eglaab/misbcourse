#
# R code for the pathway, network and text-mining analyses in the MISB course "Top-Down Systems Biology"
#
# Required input data in the working folder:
#
# - Affymetrix annotations for the HG-U133A array, downloadable from Moodle: HG-U133A.na35.annot.csv.zip
# (alternatively, this and other annotations can only be obtained via free registration on the Affymetrix website: http://www.affymetrix.com/support/technical/byproduct.affx?product=hugene-1_0-st-v1)
#
# - Pathway annotations from MSigDb, downloadable from Moodle: msigdb_pathway_annoations.zip
# (alternatively, these and other annotations can be obtained from http://software.broadinstitute.org/gsea/msigdb/ after free registration)

#
# Package installations
# (Please make sure you have a recent version of R > 3.4 - otherwise install the current software as follows:
# - Install the current version of R (3.6) for your operating system from https://ftp.gwdg.de/pub/misc/cran/
# - Install the current version of R-Studio (1.2) from: https://www.rstudio.com/products/rstudio/download/
#


# install R-package for pathway analysis
if(!require('clusterProfiler'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("clusterProfiler")
  require('clusterProfiler')
}

if(!require('GSEABase'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("GSEABase")
  require('GSEABase')
}


# install Limma package for statistical analyses
if(!require('limma'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("limma")
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

# install R annotation package
if(!require('hgu133a.db'))
{
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("hgu133a.db")
  require('hgu133a.db')
}



# set the location of your working directory (note that there are differences between Windows & Mac concerning the use of back slash "\" vs. forward slash "/")

# format for Mac & Linux systems
setwd('/Users/set/your/current/working/directory/here')

# format for Windows
setwd('C:/set/your/current/working/directory/here')



# load datasets from day 1
# (if you haven't saved the data from day 1, you can get the pre-processed data as Rdata-files from Moodle)
#

load(file="moran_preprocessed.Rdata") # moranvsn, moran_outcome_final
load(file="zhang_preprocessed.Rdata") # zhangvsn, zhang_outcome_final

# check identity of rownames between the two datasets
all(rownames(zhangvsn) == rownames(moranvsn))
#[1] TRUE


#
# Convert Gene IDs to Gene Symbols in order to be able to apply pathway analyses
#


#
# unzip Affymetrix annotation file (downloaded from Moodle into the working directory, see above):
#
# In Windows:
# - unzip the file "HG-U133A.na36.annot.csv.zip" manually
#
# In Mac/Linux:
# - use the following line of code:
system('unzip HG-U133A.na36.annot.csv.zip')


# read annotations file (ignoring comments)
# get annotations directly from the GEO database (if still loaded from the previous day):
# annot = fData(gset)
annot = read.csv("HG-U133A.na36.annot.csv", comment.char="#")
head(annot[,c(1,14:15)])

# map probes to microarray rownames
mapids = match(rownames(zhangvsn), annot$Probe.Set.ID)

# check if all IDs were mapped successfully
any(is.na(mapids))
#[1] FALSE 
# ok, no missing IDs

# extract gene symbols corresponding to microarray Probe IDs (take always the first symbol mapped)
gene_symbols = sapply( as.character(annot$Gene.Symbol[mapids]) , function(x) strsplit(x, " /// ")[[1]][1])



#
# Convert expression matrix with Affymetrix IDs to Gene Symbol matrix (if multiple probes match to a gene, take the max. average value probe as representative for the gene) 
# 

# Function to convert probe-based expression matrix to gene-based expression matrix
# Parameters:
#   matdat = input matrix with probe rownames,
#   mat_conv = vector with gene symbols corresponding to probe rownames (NA for missing conversions)
probe2genemat <- function(matdat, mat_conv)
{
  
  if(nrow(matdat) != length(mat_conv))
  {
    stop("Matrix does not have the same number of rows as the gene name vector")
  }
  
  # take max expression vector (max = maximum of mean exp across samples), if gene occurs twice among probes
  unq <- unique(mat_conv)
  if(any(is.na(unq))){
    unq <- unq[-which(is.na(unq))]
  }
  mat <- matrix(0.0, nrow=length(unq), ncol=ncol(matdat))
  for(j in 1:nrow(mat))
  {
    ind <- which(unq[j]==mat_conv)
    
    # show conversion progress, every 1000 probes
    if(j %% 1000 == 0){ 
      print(j)
    }
    
    # 1-to-1 probe to gene symbol matching
    if(length(ind) == 1)
    {
      mat[j,] = as.numeric(as.matrix(matdat[ind,]))
    } else if(length(ind) > 1){
      
      # multiple probes match to one gene symbol
      curmat = matdat[ind,]
      
      # compute average expression per row -> select row with max. avg. expression
      avg = apply(curmat, 1, mean)
      mat[j,] = as.numeric(as.matrix(matdat[ind[which.max(avg)],]))
    }
  }
  rownames(mat) = unq
  
  return(mat)
}

# Run the conversion from probe matrix to gene matrix (Zhang data)
zhang_symb = probe2genemat(zhangvsn, gene_symbols)
# get the original column names
colnames(zhang_symb) = colnames(zhangvsn)

# show the dimensions of the new gene expression matrix
dim(zhang_symb)
# the matrix has less rows than the probe expression matrix, as expected

# Run the conversion from probe matrix to gene matrix (Moran data)
moran_symb = probe2genemat(moranvsn, gene_symbols)
colnames(moran_symb) = colnames(moranvsn)

dim(moran_symb)



#
# Differential expression analysis at the gene level (instead of probe level)
#

# Zhang et al. DEG Analysis

zhang_label = ifelse(zhang_outcome_final == "disease state: Control","control","parkinson")
design <- model.matrix(~ -1+factor(zhang_label))
colnames(design) <- unique(zhang_label)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(zhang_symb, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_zhang <- topTable(eb, n = nrow(zhang_symb)) 

head(ttable_zhang)

# no. of genes with adj. p-value below 0.05
length(which(ttable_zhang$adj.P.Val < 0.05))

zhang_degs = rownames(ttable_zhang)[which(ttable_zhang$adj.P.Val < 0.05)]



# Moran et al. DEG Analysis

moran_label = moran_outcome_final
design <- model.matrix(~ -1+factor(moran_label))
colnames(design) <- unique(moran_label)

# compute simple linear model fit to microarray data (not robust)
fit <- lmFit(moran_symb, design)

contrast.matrix = makeContrasts(parkinson-control, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

eb <- eBayes(fit2)

# extract the ranking table and show the top-ranked genes
ttable_moran <- topTable(eb, n = nrow(moran_symb)) 

head(ttable_moran)

# no. of genes with adj. p-value below 0.05
length(which(ttable_moran$adj.P.Val < 0.05))

moran_degs = rownames(ttable_moran)[which(ttable_moran$adj.P.Val < 0.05)]


#
# Pathway analyses
#


# Load pathway definitions from MSigDB:
# Decompress the file msigdb_pathway_annoations.zip obtained from Moodle (see above)
#

msigdb_go_pathways = read.gmt("c5.all.v6.2.symbols.gmt")
msigdb_kegg_pathways = read.gmt("c2.cp.kegg.v6.2.symbols.gmt")
msigdb_reactome_pathways = read.gmt("c2.cp.reactome.v6.2.symbols.gmt")
msigdb_biocarta_pathways = read.gmt("c2.cp.biocarta.v6.2.symbols.gmt")
msigdb_positional = read.gmt("c1.all.v6.2.symbols.gmt")

# Inspect top of the pathway annotation table for GO:
head(msigdb_go_pathways)


#
# Apply classical Fisher's Exact test (significance-of-overlap computation) - Zhang et al.
#

# in case of lack of memory:
# rm(zhangvsn)
# rm(moranvsn)
# gc()

fisher_go_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_go_pathways)
head(fisher_go_zhang[,1:6])

fisher_kegg_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_kegg_pathways)
head(fisher_kegg_zhang[,1:6])

fisher_biocarta_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_biocarta_pathways)
head(fisher_biocarta_zhang[,1:6])
# no gene can be mapped

fisher_reactome_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_reactome_pathways)
head(fisher_reactome_zhang[,1:6])

fisher_positional_zhang <- enricher(zhang_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_positional)
head(fisher_positional_zhang[,1:6])


#
# Apply classical Fisher's Exact test (significance-of-overlap computation) - Moran et al.
#

fisher_go_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_go_pathways)
head(fisher_go_moran[,1:6])

fisher_kegg_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_kegg_pathways)
head(fisher_kegg_moran[,1:6])

fisher_reactome_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_reactome_pathways)
head(fisher_reactome_moran[,1:6])

fisher_positional_moran <- enricher(moran_degs, universe = gene_symbols, pAdjustMethod = "BH", pvalueCutoff=1.0, qvalueCutoff = 0.2, TERM2GENE = msigdb_positional)
head(fisher_positional_moran[,1:6])




#
# Apply GSEA - Zhang et al.
#

ranked_genelst = ttable_zhang$B
names(ranked_genelst) = rownames(ttable_zhang)

gsea_go_zhang = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                     maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_go_pathways,
                     TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_go_zhang[,1:7])

gsea_kegg_zhang = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                       maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_kegg_pathways,
                       TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_kegg_zhang[,1:7])

gsea_reactome_zhang = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                           maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_reactome_pathways,
                           TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_reactome_zhang[,1:7])  

gsea_positional_zhang = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                             maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_positional,
                             TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_positional_zhang[,1:7])


#
# Apply GSEA - Moran et al.
#

ranked_genelst = ttable_moran$B
names(ranked_genelst) = rownames(ttable_moran)

gsea_go_moran = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                     maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_go_pathways,
                     TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_go_moran[,1:7])

gsea_kegg_moran = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                       maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_kegg_pathways,
                       TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_kegg_moran[,1:7])

gsea_reactome_moran = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                           maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_reactome_pathways,
                           TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_reactome_moran[,1:7])  

gsea_positional_moran = GSEA(ranked_genelst, exponent = 1, minGSSize = 10,
                             maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", TERM2GENE = msigdb_positional,
                             TERM2NAME = NA, verbose = TRUE, seed = 1, by = "fgsea")
head(gsea_positional_moran[,1:7])


#
# Note: Since we used multiple pathway databases and multiple tests, we would normally still need to adjust the p-values across all databases and tests!
# (here omitted to save time)
#


#
# Network-based pathway analysis: EnrichNet
#


#
# - Goto www.enrichnet.org
# - Copy up top 100 genes to the EnrichNet web-interface (using the code below to move gene names to the clipboard)
# - Set the Identifier format to "HGNC symbol" and try out a few analyses
#

# Zhang et al.

# Copy to clipboard 

# Mac version
clip <- pipe("pbcopy", "w")                       
write.table(zhang_degs, file=clip, sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE) 
close(clip)

# Windows version
write.table(zhang_degs, "clipboard", sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE) 


# Moran et al.

# Mac version
clip <- pipe("pbcopy", "w")                     
write.table(moran_degs[1:100], file=clip, sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE) 
close(clip)

# Windows version
write.table(moran_degs[1:100], "clipboard", sep = '\t', row.names = FALSE, col.names= FALSE, quote = FALSE) 


#
# Cytoscape / Jepetto: Download Cytoscape from https://cytoscape.org/download.html and install,
# then start Cytoscape, go to "Apps" - "Apps Manager", search for "JEPETTO", install Plugin "JEPETTO"
#


#
# Text-mining analysis: Pointwise mutual information scoring
#


# Counting term frequency in pubmed (requires internet access!)
pubmed_count <- function(query)
{
  res = get_pubmed_ids(query)
  
  count <- res$Count
  return(as.numeric(count))
}


# pointwise mutual information between two PubMed search terms
pmi  <- function(term1, term2)
{
  # get current size of Pubmed (for normalization)
  cur_size <- pubmed_count("1800:2100[dp]")
  
  count1 <- pubmed_count(term1)
  count2 <- pubmed_count(term2)
  count12 <- pubmed_count(paste(term1,'+AND+',term2,sep=""))
  #count12 <- pubmed_count(gsub(" ","+",paste(term1,'AND',term2,sep="+")))
  
  return( log( count12 * cur_size / (count1 * count2)) )
}


# Zhang et al.
gene_scores = numeric(length(zhang_degs))
names(gene_scores) = zhang_degs
for(i in 1:length(gene_scores))
{
  cat(paste("Current gene:",names(gene_scores)[i]))
  gene_scores[i] = pmi(names(gene_scores)[i], "Parkinson")
  cat(paste(" - PMI-score:",gene_scores[i],"\n"))
}

# Show sorted scores
gene_scores[order(gene_scores, decreasing=T)]


# Moran et al.
gene_scores = numeric(50)
names(gene_scores) = moran_degs[1:50]
for(i in 1:length(gene_scores))
{
  cat(paste("Current gene:",names(gene_scores)[i]))
  gene_scores[i] = pmi(names(gene_scores)[i], "Parkinson")
  cat(paste(" - PMI-score:",gene_scores[i],"\n"))
}

# Show sorted scores
gene_scores[order(gene_scores, decreasing=T)]

#
# Chilibot: Go to http://chilibot.net
# Search for relation between two genes / proteins / keywords (compare top-ranked genes vs. "Parkinson": 
# Try: NR4A2 & Parkinson
# Try: DRD2 & Parkinson
# Try: DDC & Parkinson
#


#
# DisGeNet: Go to http://www.disgenet.org
# Choose "Search" - select "genes"
# Copy up to 50 top-ranked genes (see code below)
# Search and go to "Summary of Gene-Disease Associations"
# Filter by condition term of interest (e.g. Parkinson's disease)
# Click on the link for the number of associated PubMed IDs (N. PMIDs)
# Click on the link for the number of associated SNPs (N. SNPs)
#


# Copy to clipboard 

# Mac version
clip <- pipe("pbcopy", "w")                       
write.table(t(zhang_degs), file=clip, sep = '::', row.names = FALSE, col.names= FALSE, quote = FALSE) 
close(clip)

# Windows version
write.table(t(zhang_degs), "clipboard", sep = '::', row.names = FALSE, col.names= FALSE, quote = FALSE) 


# Moran et al.

# Mac version
clip <- pipe("pbcopy", "w")                     
write.table(t(moran_degs[1:50]), file=clip, sep = '::', row.names = FALSE, col.names= FALSE, quote = FALSE) 
close(clip)

# Windows version
write.table(t(moran_degs[1:50]), "clipboard", sep = '::', row.names = FALSE, col.names= FALSE, quote = FALSE) 


#
# Other resources for enrichment analysis (require only gene list as input):
# DAVID website for pathway enrichment analysis: https://david.ncifcrf.gov/summary.jsp
# g:Profiler -- a web server for functional enrichment analysis and conversions of gene lists: https://biit.cs.ut.ee/gprofiler/gost 
#

# Copy to clipboard for DAVID and g:Profiler website:

# Zhang et al.
                      
# Mac version
clip <- pipe("pbcopy", "w")                       
write.table(t(zhang_degs), file=clip, sep = '\n', row.names = FALSE, col.names= FALSE, quote = FALSE) 
close(clip)

# Windows version
write.table(t(zhang_degs), "clipboard", sep = '\n', row.names = FALSE, col.names= FALSE, quote = FALSE) 

                      
# Moran et al.

# Mac version
clip <- pipe("pbcopy", "w")                     
write.table(t(moran_degs[1:50]), file=clip, sep = '\n', row.names = FALSE, col.names= FALSE, quote = FALSE) 
close(clip)

# Windows version
write.table(t(moran_degs[1:50]), "clipboard", sep = '\n', row.names = FALSE, col.names= FALSE, quote = FALSE) 



# optionally, save the session
save.image()
# reload the session data
#load(".RData")


# For reproducibility: show and save information on all loaded R package versions
sessionInfo()
