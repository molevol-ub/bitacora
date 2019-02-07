#!/usr/bin/perl
use strict; use warnings;


die "Usage: insert the nucleotide fasta to translate and OUTNAME\n" unless @ARGV == 2;

my %genetic_code = (
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => 'X',    # Stop
    'TAG' => 'X',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => 'X',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
#Codon degeneracy
    'GCN' => 'A',
    'CGN' => 'R',
    'MGR' => 'R',
    'AAY' => 'N',
    'GAY' => 'D',
    'TGY' => 'C',
    'CAR' => 'Q',
    'GAR' => 'E',
    'GGN' => 'G',
    'CAY' => 'H',
    'ATH' => 'I',
    'YTR' => 'L',
    'CTN' => 'L',
    'AAR' => 'K',
    'TTY' => 'F',
    'CCN' => 'P',
    'TCN' => 'S',
    'AGY' => 'S',
    'ACN' => 'T',
    'TAY' => 'Y',
    'GTN' => 'V',
    'TAR' => 'X',
    'TRA' => 'X',
);

open (Fasta, "<$ARGV[0]");

my $contig = ();
my %fasta = ();

while (<Fasta>) {
	chomp;
	if ($_ =~ />(\S+)/) {
		$contig = $1;
	} else {
		$fasta{$contig} .= $_;
	}
}

close Fasta;

my $out = $ARGV[1];

open (frame_1, ">$out");


my (%frame_1, %frame_2, %frame_3) = ();	# frame_1 --> 0, _2 +1, _3 +2 [forward]
my (%frame_4, %frame_5, %frame_6) = ();	# 4 0 5 -1 6 -2 [reverse]

foreach my $contig (sort keys %fasta) {
	print frame_1 ">$contig\n";
	my $length = length $fasta{$contig};
#ahora la forward
	my $seq = uc $fasta{$contig};
	#if ($seq =~ /N/){
	#	warn "Hay Ns en la sequencia: $contig'";
	#}
	for (my $i = 0; $i < ($length - 2); $i += 3) {
		my $codon = substr ($seq, $i, 3);
		if (exists($genetic_code{$codon})){
			print frame_1 "$genetic_code{$codon}";
		}
		else {
			print frame_1 "X";
		}
    }

	print frame_1 "\n";
}

