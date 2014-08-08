#!/usr/bin/perl
use warnings;
use strict;

sub reverse_complement {
        my $dna = shift;

	# reverse the DNA sequence
        my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}


my $in_file = $ARGV[0];

open(IN_FILE,$in_file) or die "Can't open ".$in_file."\n";

while(!eof(IN_FILE)){

    my $line = <IN_FILE>;

    if($line =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*$/){
		my $chr=$2;
		my $start=$3;
		my $end=$4;
		my $rsId=$5;
		my $score=$6;
		my $strand=$7;
		my $refNCBI=$8;
		my $refUCSC=$9;
		my $observed=$10;
		my $molType=$11;
		my $class=$12;

		if(!($class eq "single")){
			if($observed =~ /^([ACGT-]+)\/([ACGT-]+)$/){ ## If two alleles are observed
				my $alt="undefined";
				my $base1 = $1;
				my $base2 = $2;
				if($refUCSC eq $base1){
					$alt = $base2;
				}
				if($refUCSC eq $base2){
					$alt = $base1;
				}
				if($alt eq "undefined"){
					$base1 = reverse_complement($base1);
					$base2 = reverse_complement($base2);
					if($refUCSC eq $base1){
						$alt = $base2;
					}
					if($refUCSC eq $base2){
						$alt = $base1;
					}
				}

				if(!($alt eq "undefined")){
					print "$chr|$end|$refUCSC|$alt\t$rsId\n";
				}

			}
		}
	}
}

## Close infile
close(IN_FILE);

