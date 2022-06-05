## get the proportion of gene regions made up of the chromatin state annotations ##

#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb

######################
## set up workspace ##

library(tidyverse)

# read in gene ensembl IDs and symbols
hg19genes <- read.table( "./supportingFiles/hg19_ensembl_prot_linc_basicGeneInfo.txt", sep = "\t", header = T)

# read in the chromHMM and variance rank metric tissue information and make a new df with the relevant info
chromHMMinfo <- read.table( "./supportingFiles/chromHMMinfo.txt", sep = "\t", header = T)
varrankTissueInfo <- data.frame( varrankTissue = chromHMMinfo$rankMetric, chromHMMtissue = chromHMMinfo$chromHMMannotation)
varrankTissueInfo$chromHMMtissue[ which( varrankTissueInfo$varrankTissue == "crossTissue")]  <- c( "genome") ## rename to match old code

# read in the information with chromHMM states and the corresponding functional annotation categories
# note: this file can be edited to define chromHMM state annotation categories (e.g., promoter, enhancer, weak promoter) as one wishes
# refer to Vu & Ernst (2022) Genome Biology for the universal annotation chromatin state definitions; and Ernst & Kellis (2015) Nature Biotechnology for the 25-state model chromatin state definitions
chromHMMann <- read.table( "./supportingFiles/chromHMM_state_categories.txt", sep = "\t", header = T)


###############################################
## summarize chromHMM states in gene regions ##

for ( j in varrankTissueInfo$chromHMMtissue) {
  
  ## match the chromHMM state names to functional annotation categories
  
  if ( j == "genome") { # universal chromHMM  
    
    thisTissueAnn <- read.delim( paste0( "./hg19geneRegions_bpENCODEchromHMM_genome_", args[ 1], "kb.txt"), sep = "\t", header = F)
    colnames( thisTissueAnn) <- c( "chr", "start", "end", "ensID", "tissue", "ann", "bp")
    
    ann <- unique( thisTissueAnn$ann) # get the list of chromHMM states
    ann <- sub( "'", "", ann)
    
    # create the df to indicate which chromHMM state corresponds to the functional annotation
    annDF <- as.data.frame( matrix( unlist (strsplit( ann, "_", fixed = T)), ncol = 2, byrow = T))
    colnames( annDF) <- c( "num", "name")
    
    annDF <- cbind( annDF, chromHMM = ann)

    tmp <- chromHMMann[ which( chromHMMann$model == "fullStack"), ] # universal chromHMM states
    
    annDF <- merge( annDF, tmp, by.x = 1, by.y = 2)
    
  } else { # tissue-level chromHMM
    
    thisTissueAnn <- read.delim( file = paste0( "./hg19geneRegions_bpENCODEchromHMM_", j, "_", args[ 1], "kb.txt"), sep = "\t", header = F)
    colnames( thisTissueAnn) <- c( "chr", "start", "end", "ensID", "tissue", "ann", "bp")
    
    ann <- unique( thisTissueAnn$ann)
    ann <- sub( "'", "", ann)
    
    # create the df to indicate which chromHMM state corresponds to the functional annotation
    annDF <- as.data.frame( matrix( unlist (strsplit( ann, "_", fixed = T)), ncol = 2, byrow = T))
    colnames( annDF) <- c( "num", "name")
    
    annDF <- cbind( annDF, chromHMM = ann)
    
    tmp <- chromHMMann[ which( chromHMMann$model == "25state"), ] # 25-state model chromHMM states
  
    annDF <- merge( annDF, tmp, by.x = 1, by.y = 2)
  }
  
  ## get the proportion of each gene region comprised of each of the functional annotations from chromHMM
  
  finalDF <- unique( data.frame( ensID = hg19genes$ensID))
  
  for ( a in unique( annDF$ann)) {
    
    thisAnn <- thisTissueAnn[ which( thisTissueAnn$ann %in% annDF$chromHMM[ which( annDF$ann == a)]), ]
    
    geneL <- unique( data.frame( cbind( ensID = thisAnn$ensID, geneL = thisAnn$end - thisAnn$start)))
    
    tmpDF <- aggregate( thisAnn$bp, by = list( ensID = thisAnn$ensID), FUN = sum)
    tmpDF <- merge( tmpDF, geneL, by = 1) %>% mutate( prop = as.numeric( x) / as.numeric( geneL)) %>% dplyr::select( ensID, prop) 
    
    finalDF <- merge( finalDF, tmpDF, by = 1, all.x = T)
    
    colnames( finalDF)[ ncol( finalDF)] <- a
  }
  
  assign( paste0( "chromHMMcov.", j), finalDF)
  
}

save.image( file = paste0( "summarizeAnnotations_propGeneRegions_out_", args[ 1], "kb.RData"))

