#!/usr/bin/perl

use warnings;
use strict;

# this program filters a vcf file based on overall sequence coverage, number of non-reference reads, number of alleles, and reverse orientation reads

# usage vcfFilter.pl infile.vcf
#
# change the marked variables below to adjust settings
#

#### stringency variables, edits as desired
my $minCoverage = 160; # minimum number of sequences; DP, 234 = 2X, N = 117 trees
my $minAltRds = 1; # minimum number of sequences with the alternative allele; AC
my $notFixed = 1.0; # removes loci fixed for alt; AF
my $bqrs = 0.005; # Mann-Whitney U test of Base Quality Bias; BQB
my $mqrs = 0.005; # Mann-Whitney U test of Mapping Quality Bias; MQB
my $rprs = 0.01; # Mann-Whitney U test of Read Position Bias; RPB
my $mq = 30; # minimum mapping quality; MQ
my $miss = 16; # maximum number of individuals with no data, 20% of N, 80% with data
my $d;

my @line;

my $in = shift(@ARGV);
open (IN, $in) or die "Could not read the infile = $in\n";
$in =~ m/^([a-zA-Z_0-9\-]+\.vcf)$/ or die "Failed to match the variant file\n";
open (OUT, "> filtered_$1") or die "Could not write the outfile\n";

my $flag = 0;
my $cnt = 0;

while (<IN>){
        chomp;
        $flag = 1;
        if (m/^\#/){ ## header row, always write
                $flag = 1;
        }
        elsif (m/^Pot/){ ## this is a sequence line, you migh need to edit this reg. expr.
                $flag = 1;
                $d = () = (m/0,0,0:0:0,0/g); ## for bcftools call
                if ($d >= $miss){
                        $flag = 0;
                        ##print "fail missing : ";
                }
                if (m/[ACTGN]\,[ACTGN]/){ ## two alternative alleles identified
                        $flag = 0;
                        #print "fail allele : ";
                }
                @line = split(/\s+/,$_);
                if(length($line[3]) > 1 or length($line[4]) > 1){
                        $flag = 0;
                        #print "fail INDEL : ";
                }
                m/DP=(\d+)/ or die "Syntax error, DP not found\n";
                if ($1 < $minCoverage){
                        $flag = 0;
                        #print "fail DP : ";
                }
## bcftools call version

                m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
                if(($1 + $2) < $minAltRds){
                        $flag = 0;
                }
                ## filter if forward and reverse reads
                m/DP4=\d+,\d+,(\d+),(\d+)/ or die "Syntax error DP4 not found\n";
                if(($1 > 0) and ($2 > 0)){
                        $flag = 0;
                }

                ## fitler out fixed
                m/AC=(\d+),*[0-9,]*;AN=(\d+)/ or die "Syntax error, AF not found: $_\n";
                if ($1==$2){
                        $flag = 0;
                #       print "fail AF : ";
                }

## bcftools call verions, these are p-values, use 0.01
                if(m/BQB=([0-9\-\.]*)/){
                        if ($1 < $bqrs){
                                $flag = 0;
#                               print "fail BQRS : ";
                        }
                }
                if(m/MQB=([0-9\-\.]*)/){
                        if ($1 < $mqrs){
                                $flag = 0;
#                               print "fail MQRS : ";
                        }
                }
                if(m/RPB=([0-9\-\.]*)/){
                        if ($1 < $rprs){
                                $flag = 0;
#                               print "fail RPRS : ";
                        }
                }
                if(m/MQ=([0-9\.]+)/){
                        if ($1 < $mq){
                                $flag = 0;
#                               print "fail MQ : ";
                        }
                }
                else{
                        $flag = 0;
                        print "faile no MQ : ";
                }
                if ($flag == 1){
                        $cnt++; ## this is a good SNV
                }
        }
        else{
                print "Warning, failed to match the chromosome or scaffold name regular expression for this line\n$_\n";
                $flag = 0;
        }
        if ($flag == 1){
                print OUT "$_\n";
        }
}
close (IN);
close (OUT);

print "Finished filtering $in\nRetained $cnt variable loci\n";
             
