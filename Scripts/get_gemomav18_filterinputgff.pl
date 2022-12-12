#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Excludes additional columns in GFF field 8 other than ID and Parent
# usage: perl get_gemomav18_filterinputgff.pl input.gff output.gff

my ($line, $name, $nameout);

my $ingff = $ARGV[0];
my $outgff = $ARGV[1];

open (Results, ">", $outgff);

# Reading GFF
open (GFFfile , "<", $ARGV[0]); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
#	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/ || $subline[2] =~ /mRNA/){

		my $genename = "";
		if ($subline[8] =~ /(\S+Parent=[^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {die "ERROR in get_gemoma_gff.pl: It fails detecting Parent ID in $line\n";}

		print Results "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\t$genename\n";
		
	}
	elsif ($subline[2] =~ /gene/){

		my $genename = "";
		if ($subline[8] =~ /(ID=[^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1; 
		}
		else {die "ERROR in get_gemoma_gff.pl: It fails detecting ID in $line\n";}

		print Results "$subline[0]\t$subline[1]\t$subline[2]\t$subline[3]\t$subline[4]\t$subline[5]\t$subline[6]\t$subline[7]\t$genename\n";
		
	} 
	else {
		die "Error in get_gemoma_v18_filterinputgff.pl No mRNA, gene or CDS found in $line\nInput GFF: $ingff\n";
	}
}
close GFFfile;
close Results;




