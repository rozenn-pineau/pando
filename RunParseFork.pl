#!/usr/bin/perl
# 
# manage barcode parsing 
#


use Parallel::ForkManager;
my $max = 36;
my $pm = Parallel::ForkManager->new($max);

FILES:
foreach $fq (@ARGV){
        $pm->start and next FILES; ## fork
        $bc = "mod_barcode_rmp2.csv";
        print "Parsing $fq with $bc\n";
        system "perl parse_barcodes768.pl $bc $fq\n";
        $pm->finish;
}

$pm->wait_all_children
