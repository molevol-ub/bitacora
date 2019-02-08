#!/usr/bin/perl
use strict;
use warnings;

# WARNING!! No tiene el filtro de longitud mínima que sí el de IRs ->>> SI LO HAGO, LAS QUE NO TENGAN LA LONGITUD MINIMA FUERA

#Usage: perl remove_seqs_with_multipledomains.pl GR_Ptepidariorum_newdata_allcutednr_cheliceratasearch.fasta ../../../Ptepidariorum_newdata/GR/GR_genomic_exons_proteins.fasta

my $explength = "60";

my ($name, $line, $contig);
my $allnames = "";
my %fasta; my %fastaexons;

my $file= "$ARGV[0]";
open (Results, ">", "$file\_nomultdomains.fasta");
open (Fasta , "<", $file);
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/){
		$name = $1;
	}
	else {
		$fasta{$name} .= $line;
	}
}
close Fasta;

open (Fasta , "<", $ARGV[1]);
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/){
		$name = $1;
	}
	else {
		$fastaexons{$name} .= $line;
	}
}
close Fasta;

foreach my $key (keys %fasta){ #LrecScaffold6808_2_3_2dom
	if ($key =~ /^(\S{4})(\S+)\_(\d+)dom/) {
		my $sp = $1;
		my $allname = $2;
		my $numdom = $3;
		my $explength2 = $explength * $numdom;
		my $seqlength = length $fasta{$key};
		if ($seqlength < $explength2){ # No lo divido si no tiene la longitud esperada 
			next if ($seqlength < $explength); # Tampoco la printo si no tiene la longitud minima deseada
			print Results ">$key\n$fasta{$key}\n";
			next;
		} 
		my @exons = split (/_/, $allname);
		my $seqname = "";
		foreach my $exon (@exons){
			if ($exon !~ /^\d+$/){
				$seqname = $exon;
			} else {
				my $exonname = "$seqname\_$exon";
				my $printseq = $fastaexons{$exonname};
				my $length3 = length ($fastaexons{$exonname});
				next if ($length3 < $explength); # Volver a filtrar por longitud para no colar cosas pequeñas
				print Results ">$sp";
				print Results "$seqname";
				print Results "_$exon\n$printseq\n";		
			}
		}
	
	} else {
		my $seqlength = length $fasta{$key};
		next if ($seqlength < $explength); # Tampoco la printo si no tiene la longitud minima deseada
		print Results ">$key\n$fasta{$key}\n";
	}
	
}

close Results;
