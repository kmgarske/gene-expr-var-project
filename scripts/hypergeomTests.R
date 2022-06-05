## Run the hypergeometric enrichment test for gene categories in the top/bottom 5% of the gene expression variance/mean ranks ##


#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is the crosse-tissue variance and mean rank metric file

######################
## set up workspace ##

library(tidyverse)
library(readxl)
library(tidyverse)

# read in the variance and mean rank metrics
ranks <- read.csv( file = args[ 1], header = T) %>% dplyr::select( Gene, mean, sd)
colnames( ranks)[ 3] <- c( "varrank")

# get the information on whether gene product is secreted from Protein Atlas
secreted <- read.table( "secretedProteins.txt", sep = "\t", header = T)

# read in the housekeeping genes (Hounkpe et al. 2021 NAR)
HKgenes <- read.table( "HKgenes_ensID_transID.txt", sep = "\t", header = F)
colnames( HKgenes) <- c( "transID", "ensID")
HKgenes <- cbind( HKgenes, ubiquitous = "yes")

# read in the Ensembl gene info 
hg19geneInfo <- read.table( "./supportingFiles/hg19_ensembl_prot_linc_basicGeneInfo.txt", sep = "\t", header = T) %>% dplyr::select( ensID, geneName, entrezID, chr, start, end, TSS, geneType) %>% mutate( length = end - start)

# read in the pLI information
pLI <- as.data.frame( read_excel( "nature19057-s2/nature19057-SI Table 13.xlsx", sheet = c( "Gene Constraint"))) %>% select( gene, pLI)
write.table( pLI, "geneConstraint_pLI.txt", sep = "\t", quote = F, row.names = F, col.names = T)
pLI <- merge( pLI, hg19geneInfo[ , c( "ensID", "geneName")], by.x = 1, by.y = 2) %>% dplyr::select( ensID, pLI) %>% distinct()


####################
## merge the data ##

allData <- merge( ranks, hg19geneInfo, by = 1, all.x = T) %>% dplyr::select( Gene, varrank, mean, chr, start, end, length) %>% distinct()

allData <- merge( allData, secreted, by = 1, all.x = T) 

allData <- merge( allData, HKgenes[ , c( 2, 3)], by = 1, all.x = TRUE)
allData$ubiquitous[ which( is.na( allData$ubiquitous))] <- c( "no")

allData <- merge( allData, pLI, by = 1, all.x = T) %>% distinct()

allData <- allData[ order( allData$varrank), ]
write.table( allData, "geneClassInfo_forHypGeom.txt", sep = "\t", col.names = T, row.names = F, quote = F)


#################################
## run the hypergeometric test ##

# create a dataframe with the information needed to run the hypergeometric test
# note: if you want to add more gene classes you will need to update this with the name of the class, the column name in the allData DF, and what value corresponds to the gene being in that class
dataInfo <- data.frame( category = c( "secreted", "HK", "pLI"), colname = c( "isSecreted", "ubiquitous", "pLI"), yes = c( "yes", "yes", 0.9))

# run the hypergeometric test
# note: this code will need to be updated if you edit the dataInfo df
hyperDF <- NULL

for ( i in 1:nrow( dataInfo)) {
  
  # get this category of interesting features
  thisSet <- allData[ , c( 2, which( colnames( allData) %in% dataInfo$colname[ i]))]
  thisSet <- thisSet[ complete.cases( thisSet), ]
  
  # total number of genes 
  total = nrow( thisSet)
  
  # pLI is a threshold whereas the others are yes/no
  if ( dataInfo$category[ i] == "pLI") {
    
    # m total number of genes in category (n = number of genes not in category)
    group1 = length( which( thisSet[ , 2] >= dataInfo$yes[ i]))
    
    # k total number of top/bottom 5% variance genes
    enrichSelection <- round( nrow( thisSet) * 0.05)
    # enrichSelection <- round( nrow( thisSet) * 0.1)
    group2 <- enrichSelection
    
    # x overlap
    Overlap.low = length( which( thisSet[ c( 1:enrichSelection), 2] >= dataInfo$yes[ i]))
    Overlap.high = length( which( thisSet[ c( ( nrow( thisSet) - enrichSelection + 1):nrow( thisSet)), 2] >= dataInfo$yes[ i]))
    
    # get the enrichment information
    hyperDF <- rbind( hyperDF, data.frame( category = dataInfo$category[ i], variance = "low", enrich = ( Overlap.low / group2) / ( group1 / total), phyper = phyper( Overlap.low - 1, group1, total - group1, group2, lower.tail = F)))
    hyperDF <- rbind( hyperDF, data.frame( category = dataInfo$category[ i], variance = "high", enrich = ( Overlap.high / group2) / ( group1 / total), phyper = phyper( Overlap.high - 1, group1, total - group1, group2, lower.tail = F)))
    
  } else { 
    
    group1 = length( which( thisSet[ , 2] == dataInfo$yes[ i]))
    
    # k total number of top/bottom 5% variance genes
    enrichSelection <- round( nrow( thisSet) * 0.05)
    # enrichSelection <- round( nrow( thisSet) * 0.1)
    group2 <- enrichSelection
    
    Overlap.low = length( which( thisSet[ c( 1:enrichSelection), 2] == dataInfo$yes[ i]))
    Overlap.high = length( which( thisSet[ c( ( nrow( thisSet) - enrichSelection + 1):nrow( thisSet)), 2] == dataInfo$yes[ i]))
    
    # get the enrichment information
    hyperDF <- rbind( hyperDF, data.frame( category = dataInfo$category[ i], variance = "low", enrich = ( Overlap.low / group2) / ( group1 / total), phyper = phyper( Overlap.low - 1, group1, total - group1, group2, lower.tail = F)))
    hyperDF <- rbind( hyperDF, data.frame( category = dataInfo$category[ i], variance = "high", enrich = ( Overlap.high / group2) / ( group1 / total), phyper = phyper( Overlap.high - 1, group1, total - group1, group2, lower.tail = F)))
  }
}

# save the dataframe for plotting locally
write.table( hyperDF, "top_bottom_5pGeneExpVarHypGeomEnrich.txt", sep = "\t", col.names = T, row.names = F, quote = F)



