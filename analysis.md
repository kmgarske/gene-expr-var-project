# 

## Download and pre-process the data 

### Secreted proteins

Get the Protein Atlas data and extract the secreted proteins

```bash

wget https://www.proteinatlas.org/download/proteinatlas.tsv.zip

zcat proteinatlas.tsv.zip | awk -F '\t' -v OFS="\t" '{print $3, $8}' | awk -F '\t' -v OFS="\t" '{ if ( $2 ~ /ecrete/) {$3 = "yes"} else {$3 = "no"}; print $1, $3}' | sed '1d' > er1

echo -e 'ensID\tisSecreted' | cat - er1 > secretedProteins.txt 

```

### pLI

```sh

wget https://static-content.springer.com/esm/art%3A10.1038%2Fnature19057/MediaObjects/41586_2016_BFnature19057_MOESM241_ESM.zip

unzip 41586_2016_BFnature19057_MOESM241_ESM.zip

```

### Houskeeping genes

```sh

wget https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv

```

### Genome annotations

```sh

wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baseline_v1.1_bedfiles.tgz

tar zxvf baseline_v1.1_bedfiles.tgz baseline_v1.1/DHS_peaks_Trynka.bed baseline_v1.1/Enhancer_Hoffman.bed baseline_v1.1/Repressed_Hoffman.bed baseline_v1.1/PromoterFlanking_Hoffman.bed baseline_v1.1/SuperEnhancer_Hnisz.bed baseline_v1.1/Transcribed_Hoffman.bed


```

### ChromHMM chromatin states

Get the chromHMM chromatin states for various tissues assessed in the current study

```bash

for i in `echo 'E062 E027 E106 E063 E066 E096 E010 E111'`

  do
  
    wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/${i}_25_imputed12marks_mnemonics.bed.gz
    
    gunzip ${i}_25_imputed12marks_mnemonics.bed.gz 
    
    tissue=`awk -v roadmap=${i} '$3 == roadmap {print $2}' chromHMMinfo.txt` 

    sed '1d' /Genomics/ayroleslab2/kristina/geneExpVarProject/annotationData/hg19_geneAnn/hg19_ensembl_prot_linc_basicGeneInfo.txt | awk -v OFS="\t" '{print "chr"$1,$2-10000,$3+10000,$4}' | sort -u | awk -v OFS="\t" '{if ($2 < 0) {$2 = 0}; print $0}' | bedtools intersect -wao -a - -b ${i}_25_imputed12marks_mnemonics.bed | awk -v OFS="\t" -v tissue=${tissue} '$8 != "." {print $1,$2,$3,$4,tissue,$8,$9}' > hg19geneRegions_bpENCODEchromHMM_${tissue}_10kb.txt
    
done

```

Get the universal chromatin states

```bash

wget https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg19_genome_100_segments.bed.gz

gunzip hg19_genome_100_segments.bed.gz

tissue="genome"

i="/Genomics/ayroleslab2/kristina/geneExpVarProject/annotationData/ENCODEdata/hg19_genome_100_segments.bed"

sed '1d' /Genomics/ayroleslab2/kristina/geneExpVarProject/annotationData/hg19_geneAnn/hg19_ensembl_prot_linc_basicGeneInfo.txt | awk -v OFS="\t" '{print "chr"$1,$2-10000,$3+10000,$4}' | sort -u | awk -v OFS="\t" '{if ($2 < 0) {$2 = 0}; print $0}' | bedtools intersect -wao -a - -b ${i} | awk -v OFS="\t" -v tissue=${tissue} '$8 != "." {print $1,$2,$3,$4,tissue,$8,$9}' > hg19geneRegions_bpENCODEchromHMM_${tissue}_10kb.txt

```


```R

myFile <- <path/to/suppTable13.xlsx>

library(readxl)
library(tidyverse)

pLI <- as.data.frame( read_excel( myFile, sheet = c( "Gene Constraint"))) %>% select( gene, pLI)

write.table( pLI, "geneConstraint_pLI.txt", sep = "\t", quote = F, row.names = F, col.names = T)

```
