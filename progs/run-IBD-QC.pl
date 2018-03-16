#!/usr/bin/perl

use strict;

my %imiss;
my %removed;

chdir("/dc2/wzhang01/Adarsh/outputs/QC_train_t321-0105");

open IMISS, '<', "data_miss.imiss"
        or die "Cannot open genotypes file (pairwiseIBS.genome): $!\n";
print "Reading PLINK .imiss file data_miss.imiss\n";
while(<IMISS>){
	s/^\s+//;
    my @fields = split /\s+/, $_;
    $imiss{$fields[0]}{$fields[1]} = $fields[5];
}

open GENOME, '<', "pairwiseIBS.genome"
        or die "Cannot open genotypes file (pairwiseIBS.genome): $!\n";
open OUT, '>', "fail-IBD-QC.txt";
print "Reading PLINK .genome file pairwiseIBS.genome\n";
while(<GENOME>){
    s/^\s+//;
    my @fields = split /\s+/, $_;
 	if($fields[9] > 0.1875){
 		if($imiss{$fields[0]}{$fields[1]}>$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 		elsif($imiss{$fields[0]}{$fields[1]}<$imiss{$fields[2]}{$fields[3]}){
 			unless($removed{$fields[2]}{$fields[3]}){
 				print OUT "$fields[2] $fields[3]\n";
 				$removed{$fields[2]}{$fields[3]} = 1;
 			}
 		}
 		else{
 			unless($removed{$fields[0]}{$fields[1]}){
 				print OUT "$fields[0] $fields[1]\n";
 				$removed{$fields[0]}{$fields[1]} = 1;
 			}
 		}
 	}
}
    
	

