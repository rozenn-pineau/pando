#!/usr/bin/perl
  
## This script returns a genotype matirx (locus by ind) with genotype
## means (point estimates) from a genotype likelihood file; a HW prior is used based on allele frequencies provided in a separate file

## USAGE: perl gl2genest.pl file.gl af_file.txt
use warnings;

$in = shift (@ARGV);
$af = shift (@ARGV);

## read in and store maf's
open (IN, $af) or die "read failed: $af\n";
while (<IN>){
        chomp;
        push (@af,$_);
}
close (IN);


## read through gl file and estimate genotypes
open (IN, $in) or die "read failed: $in\n";
$out = $in;
$out =~ s/gl$/txt/;
open (OUT, "> pntest_som_$out") or die;
while (<IN>){
    chomp;
    if (s/^\d+:\d+\s+//){ ## this line has genotype data, get rid of locus id
        $p = shift(@af); ## get alt. af. for this locus
        $prior[0] = 1-2*$p ; ## calculate prior probs as q^2, 2pq, p^2
        $prior[1] = 2*$p;
        $prior[2] = 0 ;
        @line = split(" ",$_);
        @gest = ();
        while (@line){
            $sum = 0;
            for $i (0..2){ ## three genotyple likelihoods for each individual
                $gl[$i] = shift(@line);
                $gl[$i] = (10 ** ($gl[$i]/-10)) * $prior[$i];
                $sum += $gl[$i];
            }
            $gest = 0;
            for $i (0..2){ ## normalize, and calculate mean genotype
                $gl[$i] = $gl[$i]/$sum;
                $gest += $i * $gl[$i];
            }
            $gest = sprintf("%.5f",$gest);
            push(@gest, $gest);
        }
        $gest = join(" ",@gest);
        print OUT "$gest\n";
    }
    else {
        print "failed to match $_\n";
    }
}
close (IN);
close (OUT);
