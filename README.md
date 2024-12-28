# pando
scripts for the Pando study

### Estimating the number of invariant bases

## Step : extract depth from alignments

#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bamDepth
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rozennpineau@uchicago.edu

module load samtools/1.9
module load bcftools/1.16

cd /uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/pando_data_2022/all_bams/pando-friends/pando

ref=/uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/genome/Potrs01-genome.fa
bam_list=pando_bam.list

#calculate depth at each site
samtools depth -f $bam_list --reference $ref -q 20 -Q 30 > pando_depth.txt

#modify to bed format
awk '{printf "%s\t%s\t%s", $1, $2, $2+1; for (i=3; i<=NF; i++) { printf("\t%s", $i) } print "" }'  pando_depth.txt > pando_depth.bed


## Step : compare to variant file to remove variant sites 

#make bed from vcf
cd /uufs/chpc.utah.edu/common/home/u6028866/Pando/entire-clone/variants/pando_only_variants/vcf
grep -v '#' filtered2xHiCov_pando_only_variants.vcf | awk '{print $1"\t"$2-1"\t"$2}' > pando_unfiltered_variant_sites_chrom_pos.bed

#compare and filter

#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=bed_intersect
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rozennpineau@uchicago.edu

#goal : report all entries in $depth that are not in $bed

module load bedtools/2.28.0

depth=/uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/pando_data_2022/all_bams/pando-friends/pando/pando_depth.bed

bed=/uufs/chpc.utah.edu/common/home/u6028866/Pando/entire-clone/variants/pando_only_variants/vcf/pando_unfiltered_variant_sites_chrom_pos.bed

bedtools intersect -v -wa -a $depth -b $bed > /uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/pando_data_2022/all_bams/pando-friends/pando/invariant_sites_depth.bed

#-v #Only report those entries in A that have no overlap in B
#-wa report original entry in A (full line)


## Step : filter sites with mean depth < 2 and less then 80% of individuals with data (in R)
##(I ran the job on the cluster using sbatch because it was too much work for an interactive session)

#!/usr/bin/env Rscript
  
#goal : to read in the invariant depth file and filter for thos with mean depth >2 and with data in more than 80% of the individuals

data <- read.table("invariant_sites_depth.bed", sep = "\t", header = F)
dm<-as.matrix(data[,-c(1:3)])
mn<-apply(dm,1,mean)
prop<-apply(dm > 0,1,mean)
keep<-mn > 2 & prop > 0.8

export <- data.frame(chrom=as.factor(data[keep,1]),
                     pos1=as.numeric(unlist(data[keep,2])),
                     pos2=as.numeric(unlist(data[keep,3])))
dim(export)
#export region file
write.table(export, "region_file_pando_invariant_sites.bed", sep = "\t", col.names = F, row.names = F, quote=FALSE)


## Step : count invariant bases 


#!/bin/sh
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=gompert-kp
#SBATCH --partition=gompert-kp
#SBATCH --job-name=countBases
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rozennpineau@uchicago.edu

#goal : extract the regions from the fasta file that are invariant, and count the number of each base

##based on the bed file for the invariant sites and the fasta file, output the invariant bases from the fasta file

module load samtools/1.16

cd /uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/pando_data_2022/all_bams/pando-friends/pando/

#modify file format from bed to region file format --> chr:from-to 

awk '{print $1 ":" $2 "-" $2}' region_file_pando_invariant_sites.bed > region_file_pando_invariant_sites.region_file


ref=/uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/genome/Potrs01-genome.fa

bed=/uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/pando_data_2022/all_bams/pando-friends/pando/region_file_pando_invariant_sites.region_file

samtools faidx --region-file $bed -o invar_bases_pando.fa $ref

cat invar_bases_pando.fa | grep -v "^>" | sort | uniq -c

#2077338 A
#1231148 C
#1230613 G
#    116 N
#2075411 T

#--> a total of #6614510 invariant bases. 

# I did the same steps for the fine scane (subclone) files and obtained : 
#1928021 A
#1143053 C
#1142243 G
#     99 N
#1925426 T

#--> a total of #6138743 invariant bases. 
