#!/usr/bin/perl

use strict;

chdir("${OUTPUT_DIR}");

open IN, '<', "clean_inds_data_test_missing.missing" or die "Cannot open missing file \n";
open OUT, '>', "fail_diffmiss_data.txt";
while(<IN>){
	s/^\s+//;
	my @fields = split /\s+/, $_;
	unless($fields[0] eq 'CHR'){
		if($fields[4] < 0.00001){
			print OUT "$fields[1]\n";
		}
	}
}
