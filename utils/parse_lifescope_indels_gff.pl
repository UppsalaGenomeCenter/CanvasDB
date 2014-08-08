#!/usr/bin/perl
use warnings;
use strict;

my $in_file = $ARGV[0];

open(IN_FILE,$in_file) or die "Can't open ".$in_file."\n";

while(!eof(IN_FILE)){

    my $line = <IN_FILE>;

	while($line =~ /^#.*$/){
		$line = <IN_FILE>;
	}

    if($line =~ /^(\S+)\s+.*Tool\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*.*$/){
		my $chr=$1;
		my $type=$2;
		my $start=$3;
		my $end=$4;
		my $unk1=$5;
		my $unk2=$6;
		my $unk3=$7;
		my $indel_str=$8;

		my $type_out = "NA";
		my $len = "NA";
		my $ref = "NA";
		my $alt = "NA";
		my $heterozygous = "NA";
		my $ref_reads = "NA";
		my $alt_reads = "NA";
		my $total_reads = "NA";

		my $rangeMin = $start;
		my $rangeMax = $end;

		if($type eq "insertion_site"){
			$type_out = "ins";
			if($indel_str =~ /ins_len=([^;]+);.*reference=([^;]+);.*allele-call=([^;]+);.*reference-reads=([^;]+);.*context-variant-reads=([^;]+);.*total-reads=([^;]+);.*tight_chrom_pos=([^;]+);.*loose_chrom_pos=([^;]+);.*$/){
				$len = $1;
				$ref = $2;
				my $allele_call = $3;
				$ref_reads = $4;
				my $context_var_reads = $5;
				$total_reads = $6;
				my $tight_range = $7;
				my $loose_range = $8;

				$alt_reads = $context_var_reads;

				if($context_var_reads =~ /^(\d+)\,.*$/){
					$alt_reads = $1;
				}

				if($allele_call =~ /^([ACGT-]+)\/([ACGT-]+)$/){
					my $allele1 = $1;
					my $allele2 = $2;
					if($allele1 eq $allele2){
						if($allele1 ne "-"){
							$alt = $allele1;
							$heterozygous = 0;
						}
					}
					else{
						if($allele1 ne "-"){
							$alt = $allele1;
							$heterozygous = 1;
						}
						if($allele2 ne "-"){
							$alt = $allele2;
							$heterozygous = 1;
						}
					}
				}

				if($tight_range =~ /^(\d+)\/(\d+)$/){
					my $tight_start = $1;
					my $tight_end = $2;
					if($tight_start < $rangeMin){
						$rangeMin = $tight_start;
					}
					if($tight_end > $rangeMax){
						$rangeMax = $tight_end;
					}
				}

				if($loose_range =~ /^(\d+)\/(\d+)$/){
					my $loose_start = $1;
					my $loose_end = $2;
					if($loose_start < $rangeMin){
						$rangeMin = $loose_start;
					}
					if($loose_end > $rangeMax){
						$rangeMax = $loose_end;
					}
				}


			}

		}
		if($type eq "deletion"){
			$type_out = "del";

			if($indel_str =~ /del_len=([^;]+);.*reference=([^;]+);.*allele-call=([^;]+);.*reference-reads=([^;]+);.*context-variant-reads=([^;]+);.*total-reads=([^;]+);.*tight_chrom_pos=([^;]+);.*loose_chrom_pos=([^;]+);.*$/){
				$len = $1;
				$ref = $2;
				my $allele_call = $3;
				$ref_reads = $4;
				my $context_var_reads = $5;
				$total_reads = $6;
				my $tight_range = $7;
				my $loose_range = $8;

				$alt_reads = $context_var_reads;

				if($context_var_reads =~ /^(\d+)\,.*$/){
					$alt_reads = $1;
				}

				if($allele_call =~ /^([ACGT-]+)\/([ACGT-]+)$/){
					my $allele1 = $1;
					my $allele2 = $2;
					if($allele1 eq $allele2){
						if($allele1 eq "-"){
							$alt = $allele1;
							$heterozygous = 0;
						}
					}
					else{
						if($allele1 eq "-"){
							$alt = $allele1;
							$heterozygous = 1;
						}
						if($allele2 eq "-"){
							$alt = $allele2;
							$heterozygous = 1;
						}
					}
				}

				if($tight_range =~ /^(\d+)\/(\d+)$/){
					my $tight_start = $1;
					my $tight_end = $2;
					if($tight_start < $rangeMin){
						$rangeMin = $tight_start;
					}
					if($tight_end > $rangeMax){
						$rangeMax = $tight_end;
					}
				}

				if($loose_range =~ /^(\d+)\/(\d+)$/){
					my $loose_start = $1;
					my $loose_end = $2;
					if($loose_start < $rangeMin){
						$rangeMin = $loose_start;
					}
					if($loose_end > $rangeMax){
						$rangeMax = $loose_end;
					}
				}


			}

		}

		$start=$start-1;

		my $output="$chr\t$start\t$end\t$ref\t$alt\t$total_reads\t$ref_reads\t$alt_reads\t$heterozygous";
		print "$output\n";
	}
}

## Close infile
close(IN_FILE);

