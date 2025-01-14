#!/usr/bin/perl
# 
# This script takes in one vcf file and grabs the depth per variant per individual. 
#
#
#

my @line = ();
my $word;
my $nind = 0;
my $nloc = 0;
my $first = 1; ## first vcf file, get ids from here
my $out = "/uufs/chpc.utah.edu/common/home/u6028866/Pando/replicate_analysis/variants/filtering_v2/101snps_80inds_depth.txt";

open (OUT, "> $out") or die "Could not write the outfile\n";


foreach my $in (@ARGV){
        open (IN, $in) or die "Could not read the vcf file\n";
        while (<IN>){
                chomp;
                ## get individual ids
                if (m/^#CHROM/ & ($first == 1)){
                        @line = split(m/\s+/, $_);
                        foreach $word (@line){
                                if ($word =~ m/(P|p)/){
                                        push (@inds, $word);
                                        $nind++;
                                }
                        }
                        #print OUT "$nind $nloc\n";
                        $word = join (" ", @inds);
                        print OUT "$word\n";
                }
                ## read genetic data lines, write gl
                elsif ((m/^Potrs(\d+)\s+(\d+)/)){# && (!m/[AGCT],[AGCT]/)){
                #       $word = "$1".":$1";
                        $cnt1++;
                        #print "line is found";
                                $nloc++;
                        @line = split(m/\s+/, $_);
                        $i = 0;
                        foreach $word (@line){
                                if ($word =~ s/^\d\/\d\://){
                                        #print "$word\n" ;
                                        $word =~ /([0-9]+),([0-9]+),([0-9]+):([0-9]+):/ or die "failed match $word \n";
                                                        #print "$1, $2, $3\n" ;
                                        print OUT " $4";

                                        }
                                elsif ($word =~ m/\.\/\./){
                                        print OUT " 0";
                                }
                                }
                        print OUT "\n";
                        }
                }

        }

close (OUT);
close (OUT2);
print "Number of loci: $nloc; number of individuals $nind\n";
print "$cnt1 .... $cnt2\n";
