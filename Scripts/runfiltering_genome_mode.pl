#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);
#########################################################
#####
#
# Run an additional filtering step
#
#####

## Input variables

my $name = $ARGV[0];
my $chemdir = $ARGV[1];
my $flength = $ARGV[2];

my @chemosensory;
system("ls $chemdir\/*_db.fasta > Query_genes.txt");
open(File, "<", "Query_genes.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my $id = "";
	if ($line =~ /([^\/]+)\_db.fasta/){
		$id = $1;
	} else {
		die "Are you sure you renamed you Query DB to QUERY_db.fasta?\nCannot find QUERY_db.fasta in $line\n";
	}
	push (@chemosensory, $id);
}
close File;

print "Clustering highly identical sequences\n\n";

foreach my $chem (@chemosensory){
	print "Filtering $chem\n";

	#run filtering step
	system ("perl $dirname/Tools/exclude_similar_sequences_infasta_andgff.pl $chem\/$chem\_genomic_proteins_trimmed.fasta $chem\/$chem\_genomic_genes_trimmed.gff3 $flength 2>>BITACORAstd.err 2>BITACORAstd.err");

	# Add new gene count to summary table	
	my $nseqs = 0;
	open (Filefasta, "<", "$chem\/$chem\_genomic_proteins_trimmed_idseqsclustered.fasta");
	while (<Filefasta>){
		chomp;
		my $line = $_;
		next if ($line !~ /\S+/);
		if ($line =~ /^>\S+/){
			$nseqs++;
		} 
	}
	close Filefasta;

	system ("cp $name\_genecounts_genomic_proteins.txt $name\_genecounts_genomic_proteins.txt.tmp");

	open (Results, ">", "$name\_genecounts_genomic_proteins.txt");	
	open (Filesum, "<", "$name\_genecounts_genomic_proteins.txt.tmp");
	print Results "Gene\/Gene Family\tNumber of putative identified genes\tNumber of identified genes clustering identical sequences\tNumber of identified genes clustering highly identical sequences and proteins shorter than $flength aa\n";
	while (<Filesum>){
		chomp;
		my $line = $_;
		next if ($line !~ /\S+/);
		next if ($line =~ /clustering/);	
		if ($line =~ /^$chem\t/){
			print Results "$line\t$nseqs\n";
		} else {
			print Results "$line\n";
		}
	}
	close Filesum;	
	close Results;

	system ("rm $name\_genecounts_genomic_proteins.txt.tmp");

}
print "\n";



