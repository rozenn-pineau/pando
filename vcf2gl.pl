#!/usr/bin/perl
#
# This scripts conversts a vcf file to a simpler format for downstream analysis. I am calling this format genetoype likelihood (gl). The first line lists: the number of individuals and loci. The next line has individual ids. This if followed by on line per SNP that gives the SNP id (scaffold, position) and the phred-scaled genotype likelihoods, three per individual. This version only prints a sub-set of loci with maf > $maf. 
#
# 16vi14, this script now prints the non-reference allele frequency for the retained loci to a file so that they can be used to estimate genotypes
#
# 20vii14, this now combines data from multiple files
# USAGE: perl subvcf2mpgl.pl maf *.vcf
#

my @line = ();
my $word;
my $nind = 0;
my $nloc = 0;
my $first = 1; ## first vcf file, get ids from here
my $out = "PS.gl";

open (OUT, "> $out") or die "Could not write the outfile\n";

if ($out =~ s/gl/txt/){
        open (OUT2, "> af_$out") or die "Count not write the alt. af file\n";
}

foreach my $in (@ARGV){
        open (IN, $in) or die "Could not read the vcf file\n";
        while (<IN>){
                chomp;
                ## get individual ids
                if (m/^#CHROM/ & ($first == 1)){
                        @line = split(m/\s+/, $_);
                        foreach $word (@line){
                                if ($word =~ m/(Pt|po)/){
                                        push (@inds, $word);
                                        $nind++;
                                }
                        }
                        #print OUT "$nind $nloc\n";
                        $word = join (" ", @inds);
                        #print OUT "$word\n";
                }
                ## read genetic data lines, write gl
                elsif ((m/^Potrs(\d+)\s+(\d+)/)){# && (!m/[AGCT],[AGCT]/)){
                        $word = "$1".":$1";
                        $cnt1++;
                        #print "$word \n";
                        if (m/AC=([0-9\.e\-]+);AN=([0-9\.e\-]+)/){ # here, AF= is being replaced by AF=AC/AN
                                #print "$1, $2\n" ;
                                $cnt2++;
                                #print $cnt2 . '\n';    
                                $palt = $1/$2;
                                #print "$palt\n";
                                if ($palt > 0.5){
                                        $p = 1 - $palt;
                                }
                                else {
                                        $p = $palt;
                                }
                                #print "$p\n";
                                if ($p >= $maf){ ## keep this locus
                                        $nloc++;
                                        print OUT "$word ";
                                        @line = split(m/\s+/, $_);
                                        $i = 0;
                                        foreach $word (@line){
                                                if ($word =~ s/^\d\/\d\://){
                                                        #print "$word\n" ;
                                                        $word =~ /([0-9]+),([0-9]+),([0-9]+)/ or die "failed match $word \n";
                                                        #print "$1, $2, $3\n" ;
                                                        print OUT " $1 $2 $3";

                                                }
                                                elsif ($word =~ m/\.\/\./){
                                                        print OUT " 0 0 0";
                                                }
                                        }
                                        print OUT "\n";
                                        print OUT2 "$palt\n"; ## print p before converting to maf
                                }
                        }
                }

        }
}
close (OUT);
close (OUT2);
print "Number of loci: $nloc; number of individuals $nind\n";
