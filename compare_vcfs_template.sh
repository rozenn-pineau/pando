!#/bin/bash
module load htslib

# friends
bgzip filtered_friends_optionIB_variants.vcf
bcftools index filtered_friends_optionIB_variants.vcf.gz 

# reps
bgzip filtered_reps_optionIB_variants.vcf 
bcftools index filtered_reps_optionIB_variants.vcf.gz 

# compare
bcftools isec -p filter_out_friends_reps filtered_reps_optionIB_variants.vcf.gz filtered_friends_optionIB_variants.vcf.gz
