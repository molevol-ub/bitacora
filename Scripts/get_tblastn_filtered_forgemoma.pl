#!/usr/bin/perl
use strict;
use warnings;

# Filter tblastn hits in the same scaffold positions that cause Gemoma to fail

# usage: perl get_tblastn_filtered_forgemoma.pl outfmt6_file

my ($line, $name, $nameout);
my (%blast, %fasta, %hits);
my %hitpos;

system ("sort -V -k 1 $ARGV[0] > $ARGV[0]\_sorted.txt");

open (Results, ">", "$ARGV[0]\_filtered.txt");
open (Blastfile , "<", "$ARGV[0]\_sorted.txt"); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	my $scahit = $subline[1];
	my $starthit = "";	
	my $endhit = "";
	if ($subline[8] < $subline[9]){
		$starthit = $subline[8];
		$endhit = $subline[9];
	} elsif ($subline[8] > $subline[9]){
		$starthit = $subline[9];
		$endhit = $subline[8];
	} else {
		next;
	}

	my $sphitpos = " $starthit\_$endhit ";

	if (exists $hitpos{$scahit}){
		if ($hitpos{$scahit} =~ /$sphitpos/){
			next;
		}
	}
	
	$hitpos{$scahit} .= "$sphitpos";

	print Results "$line\n";
}

close Blastfile;
close Results;

