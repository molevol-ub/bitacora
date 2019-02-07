#!/usr/bin/perl
use strict;
use warnings;

#Saco los GENES

my ($name, $line, $nameout, $contig, $namerep);
my %fasta;

$nameout=$ARGV[1];
open (Results, ">", "$nameout\_noduplicates.fasta");
my $file= "$ARGV[0]";
open (Fasta , "<", $file);
$name = "start";
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	}
	else {
		next if (exists $fasta{$name});
		$fasta{$name} = $line;
		print Results ">$name\n$line\n";
	}
}
close Fasta;
close Results;

