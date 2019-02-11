#!/usr/bin/perl
use strict;
use warnings;

#usage: perl get_fasta_trimmed.pl tabular_file_with_positions_totrim Originalfasta.fasta outname

my ($name, $line, $nameout, $contig, $frame, $name2);
my (%fasta, %inicio, %fin, %fastacutted, %fastanames);


my $file= "$ARGV[1]";
open (Fasta , "<", $file);
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		$fastanames{$name} = $line;
	}
	else {
		$fasta{$name} .= $line;
	}
}
close Fasta;

$nameout=$ARGV[2];

open (Results, ">", "$nameout\_proteins_trimmed.fasta");
open (Text, "<", $ARGV[0]);
while (<Text>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);
	$inicio{$subline[0]} = $subline[2];
	$fin{$subline[0]} = $subline[3];
}
close Text;

foreach my $key (sort keys %fasta){
	if (exists $inicio{$key}){
		my $longitud = $fin{$key} - $inicio{$key} +1;
		my $desplazamiento = "";
		if ($inicio{$key} == 0){
			$desplazamiento = $inicio{$key};
		}
		else {
			$desplazamiento = $inicio{$key} - 1;
		}
		my $cutseq = substr ($fasta{$key}, $desplazamiento, $longitud);

		print Results ">$key\n$cutseq\n";
	}
	else {
		print "$key was not found in the list with start end positions, something seems wrong. Debug\n";
	}

}
