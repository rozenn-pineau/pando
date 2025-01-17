#!/usr/bin/perl
# 
# filter phix by aligning to the phix reference, uses fork
#

system "module load bwa\n";
system "module load samtools\n";

use Parallel::ForkManager;
my $max = 12;
my $pm = Parallel::ForkManager->new($max);
my $phix = "/uufs/chpc.utah.edu/common/home/gompert-group4/data/phixSeqIndex/NC_001422.fasta";

FILES:
foreach $fq (@ARGV){
        $pm->start and next FILES; ## fork
        $fq =~ m/^([a-zA-Z_\-0-9]+)/ or die "failed to match $fq\n";
        $base = $1;
        print "Running alignment for $fq to phix\n";
        system "bwa aln -n 5 -l 20 -k 2 -t 4 $phix $fq -f $base.sai\n";
        system "bwa samse -f $base.sam -r \'\@RG\\tID:PHIX\' $phix $base.sai $fq\n";
        system "samtools view -S -b $base.sam > $base.bam\n";
        print "Removing reads mapped to phix for $base\n";
        system "samtools view -f4 $base.bam > $base.sam\n";
        system "samtools view -S -b $base.sam > $base.bam\n";
        print "Creating filtered fastq for $base\n";
        system "samtools bam2fq $base.bam > clean_$base.fastq\n";
        $pm->finish;
}

$pm->wait_all_children;
