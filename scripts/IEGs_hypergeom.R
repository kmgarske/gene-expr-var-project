## IEG enrichment in the tissue-level high-variance genes ## 

##########
## arguments to pass to script
##########

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb
## args[ 2] is the cross-tissue variance and mean rank metric file
## args[ 3] is the directory where all the tissue-level variance and mean rank metrics files are


#############################
## set up workspace (Argo) ##

library( tidyverse)

IEGs.lit <- read.table( "arner_IEGs_curated_clean.txt", sep = "\t", header = T) 

hg19genes <- read.table( "./supportingFiles/hg19_ensembl_prot_linc_basicGeneInfo.txt", sep = "\t", header = T)

# load the data (output from summarizeAnnotations_propGeneRegions.R)
load( file = paste0( "summarizeAnnotations_propGeneRegions_out_", args[ 1], "kb.RData"))
args = commandArgs( trailingOnly = TRUE)


########################
## get the enrichment ##

hyperDF <- NULL

for ( i in varrankTissueInfo$varrankTissue) {
  
  ## get the variance and mean ranks for the genes
  
  if ( i == "crossTissue") { # cross-tissue variance and mean ranks
    
    thisVarrank <- read.csv( file = args[ 2], header = T) %>% dplyr::select( Gene, mean, sd)
    colnames( thisVarrank)[ 3] <- c( "varrank")
    
  } else { # tissue-level variance and mean ranks
    
    thisVarrank <- read.csv( file = paste0( args[3], i, ".csv"), header = T) %>% dplyr::select( Gene, mean, sd)
    colnames( thisVarrank)[ 3] <- c( "varrank")
  }
  
  highVarGenes <- thisVarrank$Gene[ which( thisVarrank$varrank > quantile( thisVarrank$varrank, 0.95))]
  highVarGenes.sym <- unique( hg19genes$geneName[ which( hg19genes$ensID %in% highVarGenes)])
  
  expGenes <- thisVarrank$Gene
  expGenes.sym <- unique( hg19genes$geneName[ which( hg19genes$ensID %in% expGenes)])
  
  # total number of genes 
  total = length( expGenes.sym)
  
  # m total number of genes in category
  group1 = nrow( IEGs.lit)
  
  # k total number of high var genes
  group2 = length( highVarGenes.sym)
  
  # x overlap
  Overlap = length( which( IEGs.lit$Hs_symbol %in% highVarGenes.sym))
  
  hyperDF <- rbind( hyperDF, data.frame( tissue = i, enrich = ( Overlap / group2) / ( group1 / total), phyper = phyper( Overlap - 1, group1, total - group1, group2, lower.tail = F)))
}

hyperDF <- hyperDF %>% mutate( padj = p.adjust( phyper, method = "BH"))
write.table( hyperDF, "IEGs_hyperDF_res.txt", sep = "\t", quote = F, col.names = T, row.names = F)

