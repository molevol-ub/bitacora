#!/usr/bin/perl
use strict;
use warnings;

# usage: perl seqsfile.fasta

die "\nUsage: insert the protein fasta to be filtered and the minimum length (aa) required to retain sequences, i.e.:\nperl exclude_short_sequences_infasta.pl input.fasta 100\n\n" unless @ARGV == 2;

my $filtlen = "$ARGV[1]"; # Length to filter sequences
my $nametag = ""; # Name to add in out sequences


my $in = $ARGV[0];

my %fasta; my $name;

my $outname="";
if ($in =~ /(\S+)\.fa/){
	$outname = $1;
} else {
	$outname = $in;
}

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



open (Results, ">", "$outname\_shortseqsremoved.fasta");
foreach my $key (sort keys %fasta) {
	my $length3 = length($fasta{$key});
	next if ($length3 < $filtlen); # Length filter
	print Results ">$nametag"."$key\n$fasta{$key}\n";
}
close Results;

print "\nDone, output file has been saved in $outname\_shortseqsremoved.fasta\n\n";




