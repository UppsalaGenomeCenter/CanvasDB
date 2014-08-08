#!/usr/bin/perl
use warnings;
use strict;

my $in_file = $ARGV[0];

open(IN_FILE,$in_file) or die "Can't open ".$in_file."\n";

while(!eof(IN_FILE)){

    my $line = <IN_FILE>;

    if($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$/){
		my $chr=$1;
		my $start=$2;
		my $end=$3;
		my $ref=$4;
		my $alt=$5;
		my $SIFT_score=$6;
		my $PolyPhen_score=$8;
		my $PhyloP_score=$10;
		my $LRT_score=$12;
		my $MutTaster_score=$14;
		my $GERP_score=$16;

		my $output="chr$chr|$end|$ref|$alt\t$SIFT_score\t$PolyPhen_score\t$PhyloP_score\t$LRT_score\t$MutTaster_score\t$GERP_score\n";
		print $output;
	}
}

## Close infile
close(IN_FILE);

