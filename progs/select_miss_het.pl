#!/usr/bin/perl -w

$cut_het_high=0.197538707024795;
$cut_het_low=0.172850696902036;

$cut_miss=0.01;

chdir("${OUTPUT_DIR}");

open(MISSFILE,"data_miss.imiss");
open(HETFILE,"data_het.het");
@all=<HETFILE>;
chomp(@all);
open(OUT,">fail_miss_het.txt");

$line=0;
while(<MISSFILE>){
    chomp($_);

    if($line>=1){
        chomp($_);
        @parts_miss=split(/\s+/,$_);
        $missing=$parts_miss[6];

        @parts_het=split(/\s+/,$all[$line]);
        $meanHet=sprintf("%.3f", ($parts_het[5]-$parts_het[3])/$parts_het[5]);

        if($missing>$cut_miss or $meanHet>$cut_het_high or $meanHet<$cut_het_low){
            print OUT $parts_miss[1],"\t",$parts_miss[2],"\t",$missing,"\t",$meanHet,"\n";
        }
    }


    ++$line;
}
