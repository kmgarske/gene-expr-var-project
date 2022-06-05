## get the proportion of gene regions made up of the genome annotations ##


#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb


######################
## set up workspace ##

library(tidyverse)

# read in gene ensembl IDs and symbols
hg19genes <- read.table( "./supportingFiles/hg19_ensembl_prot_linc_basicGeneInfo.txt", sep = "\t", header = T)

# get the list of genome annotation files 
myFiles <- list.files( pattern = "hg19geneRegions_bp_")


##################################################
## summarize genome annotations in gene regions ##

finalDF <- unique( data.frame( ensID = hg19genes$ensID))

for ( thisFile in myFiles) { 
  
  annName <- gsub( "hg19geneRegions_bp_", replacement = "", thisFile) %>% gsub( "_10kb.txt", "", .) # name of this annotation

  tmpDF <- read.table( thisFile, sep = "\t", header = F) # read in the data for this annotation
  colnames( tmpDF) <- c( "chr", "start", "end", "ensID", "ann", "bp")
  
  geneL <- unique( data.frame( cbind( ensID = tmpDF$ensID, geneL = tmpDF$end - tmpDF$start))) # gene lengths for calculating proportions
  
  # aggregate the total number of bp spanned by this annotation for each gene and calculate the proportion of the gene region in that annotation
  annDF <- aggregate( tmpDF$bp, by = list( ensID = tmpDF$ensID), FUN = sum)
  annDF <- merge( annDF, geneL, by = 1) %>% mutate( prop = as.numeric( x) / as.numeric( geneL)) %>% select( ensID, prop) 
  
  finalDF <- merge( finalDF, annDF, by = 1, all = T)
  colnames( finalDF)[ ncol( finalDF)] <- annName
}

write.table( finalDF, paste0( "hg19geneRegions_FinucaneGenomeAnn_", args[1], "kb.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
