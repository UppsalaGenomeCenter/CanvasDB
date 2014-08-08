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

		if($class eq "single"){
			if($observed =~ /^([ACGT])\/([ACGT])$/){ ## If two alleles are observed
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

			if($observed =~ /^([ACGT])\/([ACGT])\/([ACGT])$/){ ## If three alleles are observed
				my $alt1="undefined";
				my $alt2="undefined";

				my $base1 = $1;
				my $base2 = $2;
				my $base3 = $3;

				if($refUCSC eq $base1){
					$alt1 = $base2;
					$alt2 = $base3;
				}
				if($refUCSC eq $base2){
					$alt1 = $base1;
					$alt2 = $base3;
				}
				if($refUCSC eq $base3){
					$alt1 = $base1;
					$alt2 = $base2;
				}

				if($alt1 eq "undefined"){
					$base1 = reverse_complement($base1);
					$base2 = reverse_complement($base2);
					$base3 = reverse_complement($base3);

					if($refUCSC eq $base1){
						$alt1 = $base2;
						$alt2 = $base3;
					}
					if($refUCSC eq $base2){
						$alt1 = $base1;
						$alt2 = $base3;
					}
					if($refUCSC eq $base3){
						$alt1 = $base1;
						$alt2 = $base2;
					}
				}

				if(!($alt1 eq "undefined")){
					print "$chr|$end|$refUCSC|$alt1\t$rsId\n";
				}
				if(!($alt2 eq "undefined")){
					print "$chr|$end|$refUCSC|$alt2\t$rsId\n";
				}

			}

			if($observed =~ /^([ACGT])\/([ACGT])\/([ACGT])\/([ACGT])$/){ ## If all four alleles are observed
				my $alt1="undefined";
				my $alt2="undefined";
				my $alt3="undefined";

				my $base1 = $1;
				my $base2 = $2;
				my $base3 = $3;
				my $base4 = $4;

				if($refUCSC eq $base1){
					$alt1 = $base2;
					$alt2 = $base3;
					$alt3 = $base4;
				}
				if($refUCSC eq $base2){
					$alt1 = $base1;
					$alt2 = $base3;
					$alt3 = $base4;
				}
				if($refUCSC eq $base3){
					$alt1 = $base1;
					$alt2 = $base2;
					$alt3 = $base4;
				}
				if($refUCSC eq $base4){
					$alt1 = $base1;
					$alt2 = $base2;
					$alt3 = $base3;
				}

				if(!($alt1 eq "undefined")){
					print "$chr|$end|$refUCSC|$alt1\t$rsId\n";
				}
				if(!($alt2 eq "undefined")){
					print "$chr|$end|$refUCSC|$alt2\t$rsId\n";
				}
				if(!($alt3 eq "undefined")){
					print "$chr|$end|$refUCSC|$alt3\t$rsId\n";
				}
			}
		}
	}
}

## Close infile
close(IN_FILE);

