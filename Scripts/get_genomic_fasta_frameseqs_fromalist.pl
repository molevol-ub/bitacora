#!/usr/bin/perl
use strict;
use warnings;

# Obtention of protein sequences from tblastn hits
#usage: perl get_genomic_fasta_frameseqs_fromalist.pl genome_translation(file.fasta.frame) IDtblastn_parsed_list.txt outname


my ($name, $line, $nameout, $contig);
my (%fasta, %frame);

# Reading genome translation
for (my $i = 1; $i <= 6; $i++){
	
	my $file= "$ARGV[0]$i";
	open (Fasta , "<", $file);
	while (<Fasta>) {
	chomp;
		$line = $_;
		next if ($line !~ /\S+/);
		if ($line =~ /^>(\S+)/) {
			$name = $1;
		}
		else {
			$fasta{$i}{$name} .= $line;
		}
	}
	close Fasta;

}

$nameout=$ARGV[2];

open (Results, ">", "$nameout\_genomic_exon_proteins.fasta");
open (Text, "<", $ARGV[1]);
while (<Text>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my $i="";	
	if ($line =~ /(\S+)\_exon\d+ Frame(\d)/){
		$i=$2;
		$contig = $1;
		my @subline = split (/ /, $line);
		my $key = $subline[0];
		my $seq = $fasta{$i}{$contig};

		my $longitud = $subline[4] - $subline[3] +1;
		my $desplazamiento = "";
		if ( $subline[3] == 0){
			$desplazamiento =  $subline[3];
		}
		else {
			$desplazamiento =  $subline[3] - 1;
		}
		my $cutseq = substr ($seq, $desplazamiento, $longitud);

		print Results ">"."$key"." Frame$i\n$cutseq\n";
	}
	else {
		die "ERROR in get_genomic_fasta_frameseqs_fromalist.pl: $nameout $line does not contain any frame and can not be codified\n";
	}
}
close Text;

