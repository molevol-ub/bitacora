#!/usr/bin/perl
use strict;
use warnings;


## Editar las variables que corresponden

my $name = $ARGV[0];


open (Chemcounts, ">", "$name\_Chemosensory_genescounts_summary.txt");
print Chemcounts "ChemosensoryFamily Genes_in_protseq Genes_in_genome Genes_in_genome_CDS Genes_in_genome_inProtein Total_genes Total_genes_nr Total_genes_withIPR GenesIPR>75% GenesIPR>50%\n";

my (%geneprot, %genegenome, %genegenomeprot, %genegenomeCDS, %totalgene, %iprgene, %ipr75gene, %ipr50gene, %totalgenenr);

open (File, "<", "$name\_Chemosensory_genescounts_proteins_withgenomic.txt");
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/ /, $line);
	next if ($subline[0] =~ /Chemosensory/);
	$geneprot{$subline[0]} = $subline[1];
	$genegenome{$subline[0]} = $subline[4] + $subline[7];
	$genegenomeCDS{$subline[0]} = $subline[5] + $subline[6];
	$genegenomeprot{$subline[0]} = $subline[8] + $subline[9];
	$totalgene{$subline[0]} = $subline[10];
	$totalgenenr{$subline[0]} = $subline[11];
	
}
close File;

open (File, "<", "$name\_Chemosensory_genescounts_final_proteins.txt");
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/ /, $line);
	next if ($subline[0] =~ /Chemosensory/);
	$iprgene{$subline[0]} = $subline[2];
	$ipr75gene{$subline[0]} = $subline[3];
	$ipr50gene{$subline[0]} = $subline[4];
}
close File;


foreach my $key (sort keys %geneprot){
	print Chemcounts"$key $geneprot{$key} $genegenome{$key} $genegenomeCDS{$key} $genegenomeprot{$key} $totalgene{$key} $totalgenenr{$key} $iprgene{$key} $ipr75gene{$key} $ipr50gene{$key}\n";
}

close Chemcounts;



