## get the downloaded excel file with lit-curated IEGs from Arner et al (2015) Science into proper format for data analysis


#################################
## arguments to pass to script ##

args = commandArgs( trailingOnly = TRUE)

## args[ 1] is the downloaded Excel file


######################
## set up workspace ##

library(tidyverse)
library(readxl)


################################
## read in and clean the file ##

IEGs.lit <- as.data.frame( read_xlsx( "arner_table_S5.xlsx", sheet = "Potential_IER_genes", skip = 1) %>% slice( 1:232))

write.table( IEGs.lit, "arner_IEGs_curated_clean.txt", sep = "\t", quote = F, col.names = T, row.names = F)
