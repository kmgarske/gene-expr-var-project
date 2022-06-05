## run partial correlations and the wilcoxon rank sum test to associate ranks with the chromHMM annotations ##

#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb
## args[ 2] is the crosse-tissue variance and mean rank metric file
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
load( file = paste0( "summarizeAnnotations_propGeneRegions_out_", args[1], "kb.RData"))

args = commandArgs( trailingOnly = TRUE)


#############################################
## run the rank-state partial correlations ##

pcorDF.var <- NULL
pcorDF.mean <- NULL

for ( i in varrankTissueInfo$varrankTissue) {
  
  ## get the variance and mean ranks for the genes
  
  if ( i == "crossTissue") { # cross-tisue variance and mean ranks
    
    thisVarrank <- read.csv( file = args[2], header = T) %>% dplyr::select( Gene, mean, sd)
    colnames( thisVarrank)[ 3] <- c( "varrank")
    
  } else { # tissue-level variance and mean ranks
    
    thisVarrank <- read.csv( file = paste0( args[3], i, ".csv"), header = T) %>% dplyr::select( Gene, mean, sd)
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


##################################################################
## plot the heatmaps and dotplots with the partial correlations ##

pcorDF.var <- cbind( pcorDF.var, padj = p.adjust( pcorDF.var$tissueAnnCorrP, method = "BH"))
pcorDF.mean <- cbind( pcorDF.mean, padj = p.adjust( pcorDF.mean$tissueAnnCorrP, method = "BH"))

write.table( pcorDF.var, file = paste0( "geneRegion_chromHMMstateProp_varrankPCor_all_", args[1], "kb.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table( pcorDF.mean, file = paste0( "geneRegion_chromHMMstateProp_meanrankPCor_all_", args[1], "kb.txt"), sep = "\t", col.names = T, row.names = F, quote = F)

## add heatmap labels

pcorDF.var <- cbind( pcorDF.var, sig = NA)
pcorDF.var$sig[ which( pcorDF.var$padj >= 0.05)] <- c( "X")
pcorDF.var$sig[ which( pcorDF.var$padj < 0.05)] <- c( "")

pcorDF.mean <- cbind( pcorDF.mean, sig = NA)
pcorDF.mean$sig[ which( pcorDF.mean$padj >= 0.05)] <- c( "X")
pcorDF.mean$sig[ which( pcorDF.mean$padj < 0.05)] <- c( "")

## refactor and categorize as comparing with self, universal, or other tissue chromatin annotations

pcorDF.var$tissue <- factor( pcorDF.var$tissue, levels = c( "crossTissue", "blood", "breast", "colon", "fat", "liver", "lung", "neuron", "stomach")) 
pcorDF.mean$tissue <- factor( pcorDF.mean$tissue, levels = c( "crossTissue", "blood", "breast", "colon", "fat", "liver", "lung", "neuron", "stomach")) 

pcorDF.var <- cbind( pcorDF.var, type = NA)
pcorDF.mean <- cbind( pcorDF.mean, type = NA)

for ( i in varrankTissueInfo$varrankTissue) {
  
  pcorDF.var$type[ which( pcorDF.var$tissue == i & pcorDF.var$chromHMMtissue == varrankTissueInfo$chromHMMtissue[ which( varrankTissueInfo$varrankTissue == i)])]  <- c( "self")
  pcorDF.mean$type[ which( pcorDF.mean$tissue == i & pcorDF.mean$chromHMMtissue == varrankTissueInfo$chromHMMtissue[ which( varrankTissueInfo$varrankTissue == i)])]  <- c( "self")
  
  if ( i == "crossTissue") {
    
    pcorDF.var$type[ which( pcorDF.var$tissue == i & is.na( pcorDF.var$type))] <- c( "other")
    pcorDF.mean$type[ which( pcorDF.mean$tissue == i & is.na( pcorDF.mean$type))] <- c( "other")
  } else {

    pcorDF.var$type[ which( pcorDF.var$tissue == i & pcorDF.var$chromHMMtissue == "genome")] <- c( "universal")
    pcorDF.mean$type[ which( pcorDF.mean$tissue == i & pcorDF.mean$chromHMMtissue == "genome")] <- c( "universal")
    
    pcorDF.var$type[ which( pcorDF.var$tissue == i & is.na( pcorDF.var$type))] <- c( "other")
    pcorDF.mean$type[ which( pcorDF.mean$tissue == i & is.na( pcorDF.mean$type))] <- c( "other")
  }
}

## plot tissue- and universal-level correlation heatmap

pcorDF.var.self <- pcorDF.var[ which( pcorDF.var$type == "self"), ]
pcorDF.mean.self <- pcorDF.mean[ which( pcorDF.mean$type == "self"), ]

gg1 <- ggplot( pcorDF.var.self, aes( x = ann, y = tissue)) + geom_tile( aes( fill = tissueAnnCorr)) + geom_text( aes( label = sig)) + theme_bw() + scale_fill_gradient2( low = "blue", mid = "white", high = "red", limits = c( min( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)), max( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)))) + ylab( "Tissue of gene expression variance rank") + xlab( "Tissue-level chromatin state") + labs( fill = "rho") + theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1))

gg2 <- ggplot( pcorDF.mean.self, aes( x = ann, y = tissue)) + geom_tile( aes( fill = tissueAnnCorr)) + geom_text( aes( label = sig)) + theme_bw() + scale_fill_gradient2( low = "blue", mid = "white", high = "red", limits = c( min( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)), max( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)))) + ylab( "Tissue of gene expression mean rank") + xlab( "Tissue-level chromatin state") + labs( fill = "rho") + theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1))

pdf( file = paste0( "tissueLevel_var_meanRank_association_", args[ 1], "kb.pdf"), width = 12, height = 8)
grid.arrange( gg1, gg2, nrow = 1)
dev.off()

## plot all tissue, cross-tissue, and universal comparisons

pcorDF.var$sig[ which( pcorDF.var$sig == "")] <- c( "significant")
pcorDF.var$sig[ which( pcorDF.var$sig == "X")] <- c( "notSignificant")

pcorDF.mean$sig[ which( pcorDF.mean$sig == "")] <- c( "significant")
pcorDF.mean$sig[ which( pcorDF.mean$sig == "X")] <- c( "notSignificant")

pcorDF.var$type[ which( pcorDF.var$tissue == "crossTissue" & pcorDF.var$chromHMMtissue == "genome")] <- c( "universal")
pcorDF.mean$type[ which( pcorDF.mean$tissue == "crossTissue" & pcorDF.mean$chromHMMtissue == "genome")] <- c( "universal")


# took this from online to shift the legend into the empty space 
shift_legend2 <- function(p) {
  # ...
  # to grob
  gp <- ggplotGrob(p)
  facet.panels <- grep("^panel", gp[["layout"]][["name"]])
  empty.facet.panels <- sapply(facet.panels, function(i) "zeroGrob" %in% class(gp[["grobs"]][[i]]))
  empty.facet.panels <- facet.panels[empty.facet.panels]
  
  # establish name of empty panels
  empty.facet.panels <- gp[["layout"]][empty.facet.panels, ]
  names <- empty.facet.panels$name
  # example of names:
  #[1] "panel-3-2" "panel-3-3"
  
  # now we just need a simple call to reposition the legend
  # reposition_legend(p, 'center', panel=names)
  reposition_legend(p, 'center', panel=names)
}

p <- ggplot( pcorDF.var, aes( x = tissue, y = tissueAnnCorr, color = type, shape = sig)) + geom_point( size = 5, alpha = 0.5) + geom_hline( yintercept = 0, linetype = "dashed", color = "grey70") + facet_wrap( ~ ann) + theme_bw() + theme( legend.text = element_text( size = 16), axis.text.x = element_text( size = 12, angle = 90, hjust = 1, vjust = 1), legend.title = element_text( size = 16), axis.text.y = element_text( size = 12), axis.title = element_text( size = 14)) + xlab( "Tissue of gene expression variance rank") + ylab( "Correlation of variance rank with indicated chromatin state") + labs( color = "Tissue-level chromatin state\nrelative to tissue variance rank", shape = "Significance of correlation\n(Padj<0.05)") + ylim( c( min( c( pcorDF.var$tissueAnnCorr, pcorDF.mean$tissueAnnCorr)), max( c( pcorDF.var$tissueAnnCorr, pcorDF.mean$tissueAnnCorr)))) + scale_color_brewer( palette = "Dark2") + theme( legend.direction = "horizontal")

gg1 <- shift_legend2( p)
dev.off()

pdf( file = paste0( "tissueLevel_varrank_association_crossTissueCorr_", args[ 1], "kb.pdf"), width = 15, height = 15)
grid.arrange( gg1, ncol = 1)
dev.off()

q <- ggplot( pcorDF.mean, aes( x = tissue, y = tissueAnnCorr, color = type, shape = sig)) + geom_point( size = 5, alpha = 0.5) + geom_hline( yintercept = 0, linetype = "dashed", color = "grey70") + facet_wrap( ~ ann) + theme_bw() + theme( legend.text = element_text( size = 16), axis.text.x = element_text( size = 12, angle = 90, hjust = 1, vjust = 1), legend.title = element_text( size = 16), axis.text.y = element_text( size = 12), axis.title = element_text( size = 14)) + xlab( "Tissue of gene expression mean rank") + ylab( "Correlation of mean rank with indicated chromatin state") + labs( color = "Tissue-level chromatin state\nrelative to tissue mean rank", shape = "Significance of correlation\n(Padj<0.05)") + ylim( c( min( c( pcorDF.var$tissueAnnCorr, pcorDF.mean$tissueAnnCorr)), max( c( pcorDF.var$tissueAnnCorr, pcorDF.mean$tissueAnnCorr)))) + scale_color_brewer( palette = "Dark2") + theme( legend.direction = "horizontal")

gg2 <- shift_legend2( q)
dev.off()

pdf( file = paste0( "tissueLevel_meanrank_association_crossTissueCorr_", args[ 1], "kb.pdf"), width = 15, height = 15)
grid.arrange( gg2, ncol = 1)
dev.off()


####################################
## run the Wilcoxon rank sum test ##

wilcoxRes.var <- NULL
wilcoxRes.mean <- NULL

thisVarrank <- read.csv( args[ 2], header = T) %>% dplyr::select( Gene, mean, sd) ## cross-tissue ranks
colnames( thisVarrank)[ 3] <- c( "varrank")

thisVarrank <- thisVarrank %>% mutate( varrankBin = ntile( thisVarrank$varrank, n = 20))
thisVarrank <- thisVarrank %>% mutate( meanrankBin = ntile( thisVarrank$mean, n = 20))

thisChromHMM <- get( "chromHMMcov.genome")

finalDF <- merge( thisVarrank, thisChromHMM, by = 1)

for ( k in 6:ncol( finalDF)) {
  
  thisCat <- colnames( finalDF)[ k] ## annotation
  
  tmpDF.tissue <- finalDF[ , c( 4, k)] ## varrank
  tmpDF.tissue <- tmpDF.tissue[ complete.cases( tmpDF.tissue), ]
  
  medianL <- median( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 1), 2]) ## lowest 5%
  medianH <- median( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 20), 2]) ## highest 5%
  
  wilcoxRes.var <- rbind( wilcoxRes.var, data.frame( category = thisCat, 
                                                     medianL = medianL, 
                                                     medianH = medianH, 
                                                     medianDiffHsubtractL = medianH - medianL, 
                                                     wilcoxP = wilcox.test( tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 20), 2], 
                                                                            tmpDF.tissue[ which( tmpDF.tissue$varrankBin == 1), 2])$p.value
                                                     )
                          )
  
  tmpDF.tissue <- finalDF[ , c( 5, k)] ## meanrank
  tmpDF.tissue <- tmpDF.tissue[ complete.cases( tmpDF.tissue), ]
  
  medianL <- median( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 1), 2]) ## lowest 5%
  medianH <- median( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 20), 2]) ## highest 5%
  
  wilcoxRes.mean <- rbind( wilcoxRes.mean, data.frame( category = thisCat, 
                                                     medianL = medianL, 
                                                     medianH = medianH, 
                                                     medianDiffHsubtractL = medianH - medianL, 
                                                     wilcoxP = wilcox.test( tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 20), 2], 
                                                                            tmpDF.tissue[ which( tmpDF.tissue$meanrankBin == 1), 2])$p.value
                                                     )
                          )
  
  
}

wilcoxRes.var <- wilcoxRes.var %>% mutate( padj = p.adjust( wilcoxP, method = "BH")) %>% mutate( sig = NA)
wilcoxRes.var$sig[ which( wilcoxRes.var$padj < 0.05)] <- c( "yes")
wilcoxRes.var$sig[ which( wilcoxRes.var$padj >= 0.05)] <- c( "no")

wilcoxRes.mean <- wilcoxRes.mean %>% mutate( padj = p.adjust( wilcoxP, method = "BH")) %>% mutate( sig = NA)
wilcoxRes.mean$sig[ which( wilcoxRes.mean$padj < 0.05)] <- c( "yes")
wilcoxRes.mean$sig[ which( wilcoxRes.mean$padj >= 0.05)] <- c( "no")


##########################################
## plot the boxplot and dotplot figures ##

## dotplot showing the median value in the top and bottom 5%

