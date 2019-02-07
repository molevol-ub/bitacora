#!/usr/bin/perl
use strict;
use warnings;

# Obtain fasta sequences from a file

# V2: If proteins contain two putative separate proteins, it cuts into X proteins: name must be X_splitN

my ($name, $line, $nameout, $contig);
my %fasta;


my $file= "$ARGV[0]";
open (Fasta , "<", $file);
while (<Fasta>) {
chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	}
	else {
		$fasta{$name} .= $line;
	}
}
close Fasta;

$nameout=$ARGV[2];

open (Results, ">", "$nameout\_proteins.fasta");
open (Text, "<", $ARGV[1]);
while (<Text>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /(\S+)/){
		$contig = $1;
		if ($contig =~ /(\S+)\_split\d/){
			my $gene = $1;
			print Results ">"."$contig"."\n$fasta{$gene}\n";
		} else {
			print Results ">"."$contig"."\n$fasta{$contig}\n";
		}
		
	}
	else {
		die "ERROR en get_fasta_fromalist_v2.pl in $line from $ARGV[1]\n";
	}
}
close Text;

