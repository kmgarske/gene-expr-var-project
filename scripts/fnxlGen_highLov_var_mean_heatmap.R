## plot the heatmap of chromHMM category


#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is distance around gene used in kb


######################
## set up workspace ##

library(tidyverse)
library(ggthemes)
library(gridExtra)
library(RColorBrewer)
library(ggrepel)

load( file = paste0( "pcor_rankVchromHMM", args[ 1], "_forHeatmap_Lineplot.RData"))

gg1 <- ggplot( pcorDF.var.self, aes( x = ann, y = tissue)) + geom_tile( aes( fill = tissueAnnCorr)) + geom_text( aes( label = sig)) + theme_bw() + scale_fill_gradient2( low = "blue", mid = "white", high = "red", limits = c( min( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)), max( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)))) + ylab( "Tissue of gene expression variance rank") + xlab( "Tissue-level chromatin state") + labs( fill = "rho") + theme_tufte( base_size = 20) + theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1))

gg2 <- ggplot( pcorDF.mean.self, aes( x = ann, y = tissue)) + geom_tile( aes( fill = tissueAnnCorr)) + geom_text( aes( label = sig)) + theme_bw() + scale_fill_gradient2( low = "blue", mid = "white", high = "red", limits = c( min( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)), max( c( pcorDF.var.self$tissueAnnCorr, pcorDF.mean.self$tissueAnnCorr)))) + ylab( "Tissue of gene expression mean rank") + xlab( "Tissue-level chromatin state") + labs( fill = "rho") + theme_tufte( base_size = 20) + theme( axis.text.x = element_text( angle = 45, vjust = 1, hjust = 1))

# pdf( file = paste0( "tissueLevel_var_meanRank_association_", args[ 1], "kb.pdf"), width = 12, height = 8)
pdf( file = paste0( "tissueLevel_var_meanRank_association_", args[ 1], "kb.pdf"), width = 15, height = 8)
grid.arrange( gg1, gg2, nrow = 1)
dev.off()


