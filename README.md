# Pando project
Scripts and notes for the Pando project. 

#  Alignments and variant calling

#  Identifying Pando samples

We calculated point estimates for each sample in the large scale dataset (i.e. the posterior mean genotype as a point estimate based on the genotype likelihood from bcftools and a binomial prior based on the allele frequency estimates from the vcf file, [vcf2gl.pl](https://github.com/rozenn-pineau/pando/blob/main/vcf2gl.pl) ).
We used principal component analysis (PCA) to ordinate the samples and k-means clustering (R kmeans function, with K=2) to label the different clusters of samples and further split the variant file into two files: the Pando variant file and the surrounding clones variant file, with 9 424 and 20 178 SNPs, respectively, see [identify_pando.R](https://github.com/rozenn-pineau/pando/blob/main/identify_pando.R). All files called in the script can be found in the supplementary material. 

#  Filtering for somatic mutations

### Step (1) : Initial quality filtering of the vcf files

### Step (2) : removing "common" mutations
To separate somatic mutations from the pool of genetic variants, we created a set of variants found in the neighboring clones and in 100 \textit{P. tremuloides} samples from the USA's Intermountain region (Colorado, Wyoming, Nevada, Idaho). We removed variants that were found in both the Pando clone samples and this comparative dataset, with the reasoning that common mutations may be germline in origin, or highly mutable sites. 

```
bcftools isec TO COMPLETE
```
Secondly, to minimize the effects of sequencing errors, we removed mutations that were found in only one sample. 


Summary table for the *replicate dataset* : 

| Step   | description                        | # mutations  | script |
| -------|:----------------------------------:| ---------:|  -------------------:|
| 1      | compare to "panel of normals"       | 332120   | [compare_vcfs_template.sh](https://github.com/rozenn-pineau/pando/blob/main/compare_vcfs_template.sh) | 
| 2      | compare to neighboring clones       | 320081   | [compare_vcfs_template.sh](https://github.com/rozenn-pineau/pando/blob/main/compare_vcfs_template.sh) | 
| 3      | filter crap (without p value tests) | 311551   | [filter_vcf_reps_1.pl](https://github.com/rozenn-pineau/pando/blob/main/vcfFilter_reps_1.pl) | 
| 4      | calculate genotype likelihood and AF | 4607    | [vcf2gl.pl](https://github.com/rozenn-pineau/pando/blob/main/vcf2gl.pl) |
| 5      | modify AF vector (anything > 0.7 is set to 0, anything between 0.5 and 0.7 is set to 0.5) | 4607     | [filter_af.R](https://github.com/rozenn-pineau/pando/blob/main/filter_af.R) |
| 6      | calculate point estimate and depth| 4607     | [gl2genet.pl](https://github.com/rozenn-pineau/pando/blob/main/gl2genet.pl) and [getDepth.sh](https://github.com/rozenn-pineau/pando/blob/main/getDepth.pl) |
| 7      | filter for individuals: mean depth > 4 | 4607   | [replicates_analysis_v3.R](https://github.com/rozenn-pineau/pando/blob/main/replicate_analysis_v3.Rmd) |
| 8      | binarize the point estimates and remove singletons | 4607   | [replicates_analysis_v3.R](https://github.com/rozenn-pineau/pando/blob/main/replicate_analysis_v3.Rmd)  |
| 9      | keep variants found in >2 samples per group and <=8 groups  | 536   | [replicates_analysis_v3.R](https://github.com/rozenn-pineau/pando/blob/main/replicate_analysis_v3.Rmd)  |
| 10     | filter vcf based on bool | 536   | [filter_vcf_based_on_bool.py](https://github.com/rozenn-pineau/pando/blob/main/filter_vcf_based_on_bool.py) |
| 11     | filter individuals with depth <4 | 536   | [replicates_analysis_v3.R](https://github.com/rozenn-pineau/pando/blob/main/replicate_analysis_v3.Rmd) |
| 12     | remove more crap using the p-value filters | 101   | [vcfFilter_536nps_80inds.pl](https://github.com/rozenn-pineau/pando/blob/main/vcfFilter_536nps_80inds.pl) |



### Step (1) : assessing our ability to recover rare mutations

 To assess our ability to consistently recover somatic mutations, we sequenced the same sample several times (12 samples sequenced 8 times each, from the same DNA extraction).

```
#!/usr/bin/Rscript
  
af <- read.table("af_101SNPs_80inds.txt", sep = "\t")

af_mod <- af
af_mod[af>0.7] <- 0
af_mod[af>0.5 & af<0.7] <- 0.5


write.table(af_mod, "afmod_101SNPs_80inds.txt", sep = "\t", col.names = F, row.names = F)
```


#  Estimating the number of invariant bases

### Step 1 : extract depth from alignments
```
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
```

### Step 2 : compare to variant file to remove variant sites 

#make bed from vcf
```
cd /uufs/chpc.utah.edu/common/home/u6028866/Pando/entire-clone/variants/pando_only_variants/vcf
grep -v '#' filtered2xHiCov_pando_only_variants.vcf | awk '{print $1"\t"$2-1"\t"$2}' > pando_unfiltered_variant_sites_chrom_pos.bed
```

#compare and filter
```
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
```

### Step 3 : filter sites with mean depth < 2 and less then 80% of individuals with data (in R)
##(I ran the job on the cluster using sbatch because it was too much work for an interactive session)
```
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
```

### Step 4 : count invariant bases 

```
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
```
```
#2077338 A
#1231148 C
#1230613 G
#116 N
#2075411 T

#--> a total of #6614510 invariant bases. 
```
### same steps for the fine-scale dataset :
```
#1928021 A
#1143053 C
#1142243 G
#99 N
#1925426 T

#--> a total of #6138743 invariant bases. 
```

#  Running BEAST2 

### Step 1 : Filtering for somatic mutations and converting to nexus format 

The nexus file was obtained by concatenating the set of somatic SNPs with binary coding of the presence of the homozygote genotype with one of the base pair (for example, "T"), a heterozygote with another base pair (for example, "A") and a missing site (no variant calling information for that site) with an "N".
See the scripts in R: 

[prep_data_for_beast_pando.Rmd](https://github.com/rozenn-pineau/pando/blob/main/prep_data_for_beast_pando.Rmd)


[prep_data_for_beast_finescale.Rmd](https://github.com/rozenn-pineau/pando/blob/main/prep_data_for_beast_finescale.Rmd)

### Step 2 : Specify model parameters using Beauti

We used the software package BEAST (version 2.7.5) to estimate the height of the phylogenetic tree for the Pando samples based on the accumulated somatic mutations; this was done on a coalescent Bayesian skyline model for effective population size. We chose the GTR nucleotide-substitution model to account for unequal substitutions rates between bases, starting with equal rates for all substitutions. We selected an optimized relaxed clock with a mean clock rate of 1. 

### Step 3 : Manually add the "invariant sites" line to the xml file

We manually added the following line to all xml after the data block, to take into account invariant bases:

```
    <data id='pando_3498snps_102inds_cstsites' spec='FilteredAlignment' filter='-' data='@pando_3498snps_102inds' constantSiteWeights='0 0 0 6614510'/>
```

### Step 4 : Run Beast 

sample code

```
#!/bin/sh 
#SBATCH --time=60:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --account=gompert
#SBATCH --partition=notchpeak
#SBATCH --job-name=beast2
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=rpineau3@gatech.edu

module load beast

beast -overwrite pando_3934snps_102inds.xml


```
All xmls used in the study are under the [xml folder](https://github.com/rozenn-pineau/pando/tree/main/xml). 


### Step 5 : Analyze .log and .trees files
We base our age estimation on the conversion of the phylogenetic tree height to an age in years. This was done in R based on the .log file. 

To visualize effective population size through time, we used Beast Traver v1.7.2 to reconstruct the Bayesian Skyline, under "Analysis", "Bayesian Skyline reconstruction" and saving the output data table. 


R scripts to plot Beast results : 





