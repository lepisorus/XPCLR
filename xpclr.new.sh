#!/bin/bash

##This script takes two command parameters: the population composition of high-low contrast and the chromosome number
##usage ./xpclrInput.sh MexLow_MexHigh.txt 4

KEEP=$1 #the high-low contrast population; this file listed all individuals in the two populations for contrast comparison 
CHR=$2 #the chromosome number
PREFIX=$(echo $1 | sed 's/\.txt//g')

#include SNPs non-singleton, biallelic, and <= 50% missing data  ## exactly the same filtering as in Hufford et al. 2012
vcftools --gzvcf /work/LAS/lwang/HMP3/HMP321/merged_flt_c${CHR}.vcf.gz --mac 2 --remove-indels --min-alleles 2 --max-alleles 2 --chr ${CHR} --max-missing 0.5 --keep teosinte.txt --recode --out teosinte.chr${CHR}
vcftools --gzvcf /work/LAS/lwang/HMP3/HMP321/merged_flt_c${CHR}.vcf.gz --mac 2 --remove-indels --min-alleles 2 --max-alleles 2 --chr ${CHR} --max-missing 0.5 --keep landrace.sorted.txt --recode --out landrace.chr${CHR}
vcftools --gzvcf /work/LAS/lwang/HMP3/HMP321/merged_flt_c${CHR}.vcf.gz --mac 2 --remove-indels --min-alleles 2 --max-alleles 2 --chr ${CHR} --max-missing 0.5 --keep improved.txt --recode --out improved.chr${CHR}

#get the list of sites, which are kept after the filtering in each population (teosinte, maize, improved lines)
cat <(less teosinte.chr${CHR}.recode.vcf | grep -v "^#" | awk '{print $1"\t"$2}') <(less landrace.chr${CHR}.recode.vcf | grep -v "^#" | awk '{print $1"\t"$2}') <(less improved.chr${CHR}.recode.vcf | grep -v "^#" | awk '{print $1"\t"$2}') | sort -n -k1,1 -n -k2,2 | uniq -c | sed 's/^[ \t]*//g' | awk -v OFS="\t" '$1==3 {print $2, $3}' | sort -n -k1,1 -n -k2,2 > chr${CHR}.sites

#a new vcf, keeping all the sites and individuals
vcftools --gzvcf /work/LAS/lwang/HMP3/HMP321/merged_flt_c${CHR}.vcf.gz --keep ${KEEP} --positions chr${CHR}.sites --recode --out ${PREFIX}.chr${CHR}

#prepare genotype file for XPCLR program
perl ../vcf2xpclrGeno.pl ${PREFIX}.chr${CHR}.recode.vcf

#convert to the genetic map
Rscript ../Phys2Genetic.R -i ${PREFIX}.chr${CHR}.snp -o ${PREFIX}.chr${CHR}.snp.geneticPos

##need to adjust the column number for High and Low Populations
cut -f1-4,51-108,137-144 ${PREFIX}.chr${CHR}.geno > improved.chr${CHR}.geno
cut -f5-50 ${PREFIX}.chr${CHR}.geno > landrace.chr${CHR}.geno
cut -f109-136 ${PREFIX}.chr${CHR}.geno > teosinte.chr${CHR}.geno

##The input genotype file is space delimited and with one space and line break at the end.
perl -p -i -e 's/\t/ /g' improved.chr${CHR}.geno
perl -p -i -e 's/\t/ /g' landrace.chr${CHR}.geno
perl -p -i -e 's/\t/ /g' teosinte.chr${CHR}.geno


#add one more space at end of each line

sed -i 's/$/ /g' improved.chr${CHR}.geno
sed -i 's/$/ /g' landrace.chr${CHR}.geno
sed -i 's/$/ /g' teosinte.chr${CHR}.geno

#scripts to run the XPCLR
/work/LAS/lwang/bin/XPCLR/XPCLR -xpclr landrace.chr${CHR}.geno teosinte.chr${CHR}.geno ${PREFIX}.chr${CHR}.snp.geneticPos dom.chr${CHR} -w1 0.0005 50 100 ${CHR} -p0 0.7
/work/LAS/lwang/bin/XPCLR/XPCLR -xpclr improved.chr${CHR}.geno landrace.chr${CHR}.geno ${PREFIX}.chr${CHR}.snp.geneticPos imp.chr${CHR} -w1 0.0005 50 100 ${CHR} -p0 0.7 


