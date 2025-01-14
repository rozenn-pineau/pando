#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script opens a vcf file and a boolean vector of the same length.
It reads every boolean value from the vector and keeps the "true" lines of the file.
"""

def filter_vcf_based_on_bool(input_file, input_vec,output_file):

    final = open(output_file,"w+")

    with open(input_file) as file, open(input_vec) as vec:
        # Go through the file line by line
        for line, boolean in zip(file,vec):
            line = line.strip()
            boolean = int(boolean)
            #print(boolean)
            if boolean == 1:
                final.write(f'{line}\n')

# Call the function here::wq

#filter_vcf_based_on_bool("/uufs/chpc.utah.edu/common/home/u6028866/Pando/replicate_analysis/variants/filtering_v2/noheader_filtered_ponfil_friendfil_replicates.vcf", 
#                        "/uufs/chpc.utah.edu/common/home/u6028866/Pando/replicate_analysis/variants/filtering_v2/bool_filter_replicate_vcf_536snps_today.txt",
#                        "/uufs/chpc.utah.edu/common/home/u6028866/Pando/replicate_analysis/variants/filtering_v2/noheader_bool_filtered_ponfil_friendfil_replicates.vcf")
