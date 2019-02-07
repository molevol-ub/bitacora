#!/usr/bin/perl
use strict;
use warnings;

# Obtain sequences from a file

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
		print Results ">"."$contig"."\n$fasta{$contig}\n";
	}
	else {
		die "ERROR in get_fasta_fromalist.pl in $line from $ARGV[1]\n";
	}
}
close Text;

