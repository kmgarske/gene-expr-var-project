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

```bash

wget https://housekeeping.unicamp.br/Housekeeping_GenesHuman.csv

## match transcript to ensID 

awk -F ';' '{print $1}' Housekeeping_GenesHuman.csv | sed '1d' | sort > er1

sort -k2,2 ./supportingFiles/gencode.ensID.transID.txt | join -1 2 -2 1 -t $'\t' - er1 > HKgenes_ensID_transID.txt

```

### Gene regions

```bash

kb=10 # set how many kb you want around the gene 

outGeneLoc="hg19_geneLoc_extend${kb}kb.bed"

if [ ! -f ${outGeneLoc} ]

  then
  
    sed '1d' ./supportingFiles/hg19_ensembl_prot_linc_basicGeneInfo.txt | awk -v OFS="\t" -v kb=${kb} '{print "chr"$1,$2-(kb * 10^3),$3+(kb * 10^3),$4}' | sort -u > er1
    
    # check chromosome ends
    
    for j in {1..22}

     do

       awk -v chr=${j} '$1 == "chr"chr' er1 > er2
        
       end=`awk -v chr=${j} '$1 == chr {print $2}' ./supportingFiles/chromsizes.txt`
        
       awk -v OFS="\t" -v end=${end} '{if ($2 < 0) {$2 = 0}; if ($3 > end) {$3 = end}; print $0}' er2 >> ${outGeneLoc}
        
    done
    
fi

```

### Genome annotations

```bash

wget https://storage.googleapis.com/broad-alkesgroup-public/LDSCORE/baseline_v1.1_bedfiles.tgz

tar zxvf baseline_v1.1_bedfiles.tgz baseline_v1.1/DHS_peaks_Trynka.bed baseline_v1.1/Enhancer_Hoffman.bed baseline_v1.1/Repressed_Hoffman.bed baseline_v1.1/PromoterFlanking_Hoffman.bed baseline_v1.1/SuperEnhancer_Hnisz.bed baseline_v1.1/Transcribed_Hoffman.bed baseline_v1.1/TFBS_ENCODE.bed baseline_v1.1/CTCF_Hoffman.bed 

## overlap with gene regions

for i in `ls ./baseline_v1.1/*bed` 

  do
  
    annotation=`echo ${i} | sed 's@./baseline_v1.1/@@' | sed 's/.bed//'`
    
    bedtools intersect -wao -a ${outGeneLoc} -b ${i} | awk -v OFS="\t" -v ann=${annotation} '{print $1,$2,$3,$4,ann,$8}' > hg19geneRegions_bp_${annotation}_${kb}kb.txt

done

```

Summarize the genome annotation information into proportion of gene regions made up of each of the annotations

Make sure you have the following R packages and their dependencies installed: tidyverse

```bash

Rscript ./scripts/summarizeAnnotations_propGeneRegions_Finucane.R ${kb}

```

### ChromHMM chromatin states

Get the chromHMM chromatin states for various tissues assessed in the current study as well as the universal chromHMM states; overlap with gene regions

```bash

mkdir chromHMMdata

cd chromHMMdata/

## tissue-level

for i in `echo 'E062 E027 E106 E063 E066 E096 E010 E111'`

  do
  
    wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/${i}_25_imputed12marks_mnemonics.bed.gz
    
    gunzip ${i}_25_imputed12marks_mnemonics.bed.gz 
    
done

## universal

wget https://public.hoffman2.idre.ucla.edu/ernst/1G6UT/hg19_genome_100_segments.bed.gz

gunzip hg19_genome_100_segments.bed.gz

cd ../

###############################
## overlap with gene regions ##

## tissue-level

for i in `ls ./chromHMMdata/E*bed`

  do
  
    id=`echo ${i} | sed 's@./chromHMMdata/@@' | sed 's/_25_imputed12marks_mnemonics.bed//'`
    
    tissue=`awk -v roadmap=${id} -F '\t' '$4 == roadmap {print $3}' ./supportingFiles/chromHMMinfo.txt`
    
    bedtools intersect -wao -a ${outGeneLoc} -b ${i} | awk -v OFS="\t" -v tissue=${tissue} '$8 != "." {print $1,$2,$3,$4,tissue,$8,$9}' > hg19geneRegions_bpENCODEchromHMM_${tissue}_${kb}kb.txt
           
done

## universal

tissue="genome"

i="./chromHMMdata/hg19_genome_100_segments.bed"

bedtools intersect -wao -a ${outGeneLoc} -b ${i} | awk -v OFS="\t" -v tissue=${tissue} '$8 != "." {print $1,$2,$3,$4,tissue,$8,$9}' > hg19geneRegions_bpENCODEchromHMM_${tissue}_${kb}kb.txt

```

Summarize the chromHMM state information into proportion of gene regions made up of each of the state functional annotations

Make sure you have the following R packages and their dependencies installed: tidyverse

```bash

Rscript ./scripts/summarizeAnnotations_propGeneRegions_chromHMM.R ${kb}

```

## Test for enrichment of various gene classes in top/bottom 5% gene expression variance/mean rank

```bash

Rscript ./scripts/hypergeomTests.R <path/to/crossTissueRank/csv>

```

## Associate gene expression variance/mean rank metric with genome annotations

Run the partial correlation of gene expression variance/mean rank metrics with each genome (Finucane) or chromHMM annotation. Check if the genes in the top/bottom 5% of expression variance/mean exhibit differences in the proportion of gene regions made up of the annotations. 

Make sure you have the following R packages and their dependencies installed: tidyverse, ppcor

```bash

Rscript ./scripts/pcor_rankVchromHMM.R ${kb} <path/to/crossTissueRank/csv> <path/to/tissueLevelRank/directory>

```
