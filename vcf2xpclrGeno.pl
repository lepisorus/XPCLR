#!/usr/bin/perl
# USAGE:
# unix command line:
use strict;
use warnings;

my $infile = shift @ARGV;    # read file names from terminal input
open IN, '<', $infile or die "Cannot open '$infile' because: $!";

my $new_file = $infile;
$new_file =~ s/\.recode\.vcf//;

my $snp_file = $new_file."\.snp";    # give a output file name
open OUT1, '>', $snp_file or die "Cannot open '$snp_file' because: $!";
print OUT1 "SNPname\tchr\tPhysicalPosition\tRefAllele\tAltAllele\n";

my $geno_file = $new_file."\.geno";
open OUT2, '>', $geno_file or die "Cannot open '$geno_file' because: $!";


while (<IN>){
chomp;
next if /^#/;
next if /^[sPM]/;
my @tmpo=split("\t", $_);

##deal with the sample geno file
my $SNPname = $tmpo[0]."\_".$tmpo[1];
print OUT1 "$SNPname\t$tmpo[0]\t$tmpo[1]\t$tmpo[3]\t$tmpo[4]\n";
my $geno;

for (my $i=9; $i <= $#tmpo; $i++){
if ($tmpo[$i] =~ /^0\/0/) {
$geno = "0\t0";
}
elsif ($tmpo[$i] =~ /^1\/1/) {
$geno = "1\t1";
}
elsif ($tmpo[$i] =~ /^0\/1/) {
$geno = "0\t1";
}
elsif ($tmpo[$i] =~ /^\.\/\./) {
$geno = "9\t9";
}
else{
print "something wrong, please check your scripts";
}
print OUT2 "$geno\t";
}#for
print OUT2 "\n";
}#while

   
close IN;
close OUT1;
close OUT2;



