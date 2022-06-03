# 

## Download the data to be used for the enrichment

Get the Protein Atlas data and extract the secreted proteins

```bash

wget https://www.proteinatlas.org/download/proteinatlas.tsv.zip

zcat proteinatlas.tsv.zip | awk -F '\t' -v OFS="\t" '{print $3, $8}' | awk -F '\t' -v OFS="\t" '{ if ( $2 ~ /ecrete/) {$3 = "yes"} else {$3 = "no"}; print $1, $3}' | sed '1d' > er1

echo -e 'ensID\tisSecreted' | cat - er1 > secretedProteins.txt 

```

Get the chromHMM chromatin states for various tissues assessed in the current study

```bash

for i in `echo 'E062 E027 E106 E063 E066 E096 E010 E111'`

  do
  
    wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/imputed12marks/jointModel/final/${i}_25_imputed12marks_mnemonics.bed.gz

```

