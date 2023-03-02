## run partial correlations and the wilcoxon rank sum test to associate ranks with the chromHMM annotations ##


#################################
## arguments to pass to script ##

args <- commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb
## args[ 2] is the cross-tissue variance and mean rank metric file
## args[ 3] is the directory where all the tissue-level variance and mean rank metrics files are


######################
## set up workspace ##

library(tidyverse)
library(ppcor)
library(gridExtra)
library(lemon)
library(RColorBrewer)
library(ggrepel)

# load the data (output from summarizeAnnotations_propGeneRegions.R)
load( file = paste0( "summarizeAnnotations_propGeneRegions_out_", args[ 1], "kb.RData"))

args <- commandArgs( trailingOnly = TRUE)


#############################################
## run the rank-state partial correlations ##

pcorDF.var <- NULL
pcorDF.mean <- NULL

for ( i in varrankTissueInfo$varrankTissue) {
  
  ## get the variance and mean ranks for the genes
  
  if ( i == "crossTissue") { # cross-tissue variance and mean ranks
    
    thisVarrank <- read.csv( file = args[ 2], header = T) %>% dplyr::select( Gene, mean, sd)
    colnames( thisVarrank)[ 3] <- c( "varrank")
    
  } else { # tissue-level variance and mean ranks
    
    thisVarrank <- read.csv( file = paste0( args[ 3], i, ".csv"), header = T) %>% dplyr::select( Gene, mean, sd)
    colnames( thisVarrank)[ 3] <- c( "varrank")
  }
  
  ## get the chromHMM functional annotation coverage in the gene regions and merge with gene ranks; run the partial  correlation
  
  for ( j in varrankTissueInfo$chromHMMtissue) {
    
    thisChromHMM <- get( paste0( "chromHMMcov.", j))
    
    finalDF <- merge( thisVarrank, thisChromHMM, by = 1)
    
    for ( k in 4:ncol( finalDF)) {
      
      tmpDF.tissue <- finalDF[ , c( 2, 3, k)]
      tmpDF.tissue <- tmpDF.tissue[ complete.cases( tmpDF.tissue), ]
      mypcor.tissue <- pcor( tmpDF.tissue, method = "spearman")
      
      pcorDF.var <- rbind( pcorDF.var, data.frame( tissue = i, chromHMMtissue = j, ann = colnames( finalDF)[ k], tissueAnnCorr = mypcor.tissue$estimate[ 2, 3], tissueAnnCorrP = mypcor.tissue$p.value[ 2, 3]))
      
      pcorDF.mean <- rbind( pcorDF.mean, data.frame( tissue = i, chromHMMtissue = j, ann = colnames( finalDF)[ k], tissueAnnCorr = mypcor.tissue$estimate[ 1, 3], tissueAnnCorrP = mypcor.tissue$p.value[ 1, 3]))
    }
  }
}


###################################
## prepare the data for plotting ##

pcorDF.var <- cbind( pcorDF.var, padj = p.adjust( pcorDF.var$tissueAnnCorrP, method = "BH"))
pcorDF.mean <- cbind( pcorDF.mean, padj = p.adjust( pcorDF.mean$tissueAnnCorrP, method = "BH"))

write.table( pcorDF.var, file = paste0( "geneRegion_chromHMMstateProp_varrankPCor_all_", args[1], "kb.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table( pcorDF.mean, file = paste0( "geneRegion_chromHMMstateProp_meanrankPCor_all_", args[1], "kb.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

## rename the crossTissue to across-study
varrankTissueInfo$varrankTissue[ which( varrankTissueInfo$varrankTissue == "crossTissue")] <- c( "across-study")

## add heatmap labels

pcorDF.var <- cbind( pcorDF.var, sig = NA)
pcorDF.var$sig[ which( pcorDF.var$padj >= 0.05)] <- c( "X")
pcorDF.var$sig[ which( pcorDF.var$padj < 0.05)] <- c( "")

pcorDF.mean <- cbind( pcorDF.mean, sig = NA)
pcorDF.mean$sig[ which( pcorDF.mean$padj >= 0.05)] <- c( "X")
pcorDF.mean$sig[ which( pcorDF.mean$padj < 0.05)] <- c( "")

## refactor and categorize as comparing with self, universal, or other tissue chromatin annotations

pcorDF.var$tissue[ which( pcorDF.var$tissue == "crossTissue")] <- c( "across-study")
pcorDF.mean$tissue[ which( pcorDF.mean$tissue == "crossTissue")] <- c( "across-study")

pcorDF.var$tissue <- factor( pcorDF.var$tissue, levels = c( "across-study", "blood", "breast", "colon", "fat", "liver", "lung", "neuron", "stomach")) 
pcorDF.mean$tissue <- factor( pcorDF.mean$tissue, levels = c( "across-study", "blood", "breast", "colon", "fat", "liver", "lung", "neuron", "stomach")) 

pcorDF.var <- cbind( pcorDF.var, type = NA)
pcorDF.mean <- cbind( pcorDF.mean, type = NA)

for ( i in varrankTissueInfo$varrankTissue) {
  
  pcorDF.var$type[ which( pcorDF.var$tissue == i & pcorDF.var$chromHMMtissue == varrankTissueInfo$chromHMMtissue[ which( varrankTissueInfo$varrankTissue == i)])]  <- c( "self")
  pcorDF.mean$type[ which( pcorDF.mean$tissue == i & pcorDF.mean$chromHMMtissue == varrankTissueInfo$chromHMMtissue[ which( varrankTissueInfo$varrankTissue == i)])]  <- c( "self")
  
  if ( i == "across-study") {
    
    pcorDF.var$type[ which( pcorDF.var$tissue == i & is.na( pcorDF.var$type))] <- c( "other")
    pcorDF.mean$type[ which( pcorDF.mean$tissue == i & is.na( pcorDF.mean$type))] <- c( "other")
  } else {

    pcorDF.var$type[ which( pcorDF.var$tissue == i & pcorDF.var$chromHMMtissue == "genome")] <- c( "universal")
    pcorDF.mean$type[ which( pcorDF.mean$tissue == i & pcorDF.mean$chromHMMtissue == "genome")] <- c( "universal")
    
    pcorDF.var$type[ which( pcorDF.var$tissue == i & is.na( pcorDF.var$type))] <- c( "other")
    pcorDF.mean$type[ which( pcorDF.mean$tissue == i & is.na( pcorDF.mean$type))] <- c( "other")
  }
}

## separate out the self-comparisons for the simple line plots
pcorDF.var.self <- pcorDF.var[ which( pcorDF.var$type == "self"), ]
pcorDF.mean.self <- pcorDF.mean[ which( pcorDF.mean$type == "self"), ]

## save this RData file so plotting the figure can be modified as needed
save.image( file = paste0( "pcor_rankVchromHMM", args[ 1], "_forHeatmap_Lineplot.RData"))