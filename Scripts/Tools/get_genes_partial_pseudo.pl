#!/usr/bin/perl
use strict;
use warnings;

# Script to classify protein annotations as complete, partial or pseudogenes based on protein length comparing with the average length for a family

# usage: perl seqsfile.fasta length

die "\nUsage: insert the protein fasta and the required length (aa) to classify proteins as complete or partial, i.e.:\nperl get_genes_partial_pseudo.pl input.fasta 100\n\n" unless @ARGV == 2;

my $filtlen = "$ARGV[1]"; # Length to classify as complete or partial
my $outname = ""; # Name to add in out sequences


my $in = $ARGV[0];

my %fasta; my $name;

open (Fasta , "<", "$in"); 
while (<Fasta>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= $line;
	}

}
close Fasta;



open (Results, ">", "$in\_genes.fasta");
open (Resultspse, ">", "$in\_pseudogenes.fasta");
open (Resultspar, ">", "$in\_partial_genes.fasta");

foreach my $key (sort keys %fasta) {
	my $length3 = length($fasta{$key});

	if ($length3 >= $filtlen){
		my $middleseq = substr ($fasta{$key}, 5, -10);
		#if ($fasta{$key} =~ /X/){
		if ($middleseq =~ /X/){
			print Resultspse ">$key\n$fasta{$key}\n";
		} else {
			print Results ">$key\n$fasta{$key}\n";
		}
	} else {
		my $middleseq = substr ($fasta{$key}, 5, -10);
		#if ($fasta{$key} =~ /X/){
		if ($middleseq =~ /X/){
			print Resultspse ">$key\n$fasta{$key}\n";
		} else {
			print Resultspar ">$key\n$fasta{$key}\n";
		}
	}

}
close Results;
close Resultspse;
close Resultspar;

print "\nDone. Results have been saved in $in\_genes.fasta, $in\_partial_genes.fasta and $in\_pseudogenes.fasta\n\n";


