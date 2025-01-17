#!/usr/bin/perl
# 
# alignment with bwa mem 
#


use Parallel::ForkManager;
my $max = 80;
my $pm = Parallel::ForkManager->new($max);
my $genome = "/uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/genome/Potrs01-genome.fa";
my $out = "/uufs/chpc.utah.edu/common/home/u6028866/Pando/replicate_analysis/aligned/";

# /uufs/chpc.utah.edu/common/home/u6028866/Pando/subclone/genome/Potrs01-genome.fa
#
FILES:
foreach $file (@ARGV){
        $pm->start and next FILES; ## fork
        if ($file =~ m/^([A-Z0-9a-z\-\_]+)/){
                $ind = $1;
        }
        else {
                die "Failed to match $file\n";
        }
        system "bwa mem -t 1 -k 15 -r 1.3 -T 30 -R \'\@RG\\tID:Ptrem-"."$ind\\tPL:ILLUMINA\\tLB:Ptrem-"."$ind\\tSM:Ptrem-"."$ind"."\' $genome $file > $out/aln"."$ind".".sam\n";

        $pm->finish;
}

$pm->wait_all_children;
