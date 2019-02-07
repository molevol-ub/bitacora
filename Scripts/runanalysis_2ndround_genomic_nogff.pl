#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);
#########################################################
#####
#
# Run second analysis, identifying new proteins in the genome
#
#####

## Input variables

my $name = $ARGV[0];
my $chemdir = $ARGV[1];
my $genome = $ARGV[2];
my $translation = "$ARGV[2].frame";
my $evalue = $ARGV[3];
my $maxintron = $ARGV[4];
my $threads = $ARGV[5];

## Start

my @chemosensory;
system("ls $chemdir\/*db.fasta > Query_genes.txt");
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


# Indexing and translating genome
print "Indexing and translating genome file $genome\n";

my $indexed = 0;
my $translated = 0;
system ("ls $genome* > index_check.txt");
open (File, "<", "index_check.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /$genome\.nin/ || $line =~ /$genome\.nhr/ || $line =~ /$genome\.nsq/){
		$indexed++;
	}
	elsif ($line =~ /$genome\.frame\d/){
		$translated++;
	}

}
close File;
system("rm index_check.txt");

if ($translated == 6){
	print "Found translation for $genome\n";
} else {
	print "Translating genome... It may take a while\n";
	system ("perl $dirname/translating_seqs_v3.pl $genome");
	print "Translation completed for $genome\n";
}

if ($indexed == 3){
	print "Found a blast index for $genome\n";
} else {
	system ("makeblastdb -dbtype nucl -in $genome");
}

#

open (Chemcounts, ">", "$name\_genecounts_genomic_proteins.txt");
print Chemcounts"Gene\/Gene Family\tNumber of putative identified genes\tNumber of putative identified genes clustering identical sequences\n";


foreach my $chem (@chemosensory){
	print "\n----------------- Starting $chem Genomic Gene Identification and Annotation\n";

	system ("mkdir -p $chem");
	system ("mkdir -p $chem\/hmmer");

	# run tblastn
	print "Doing $chem tblastn\n";
	system ("tblastn -query $chemdir\/$chem\_db.fasta -db $genome -out $chem\/$name\_Vs$chem\_tblastn\.outfmt6 -evalue $evalue -num_threads $threads -outfmt \"6 std sframe qlen slen\"");

	# Parsing tblastn output
	print "Parsing $chem tblastn\n";

	#system ("perl $dirname/get_genomic_tblastn_parsed_newv_nogff.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $evalue");
	system ("perl $dirname/get_genomic_tblastn_parsed_newv_nogff_genomic_positions.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $evalue");


	###### Exiting if no novel genes are found

	my $nonew = 0;
	open (Filetblastn, "<", "$chem\/$chem"."tblastn_parsed_list_nogff_filtered.txt");
	while (<Filetblastn>){
		chomp;
		my $nlin = $_;
		next if ($nlin !~ /\S+/);
		$nonew++;
	}
	close Filetblastn;

	if ($nonew == 0){
		print "No $chem identified in genomic regions\n";

		print Chemcounts"$chem\t0\t0\n";

		print "----------------- $chem DONE\n\n";
		next;
	}

	#######

	### Script continues
	# Obtaining protein sequences, exons and genes

	system ("perl $dirname/get_genomic_fasta_frameseqs_fromalist.pl $translation $chem/$chem"."tblastn_parsed_list_nogff_filtered.txt $chem/$chem");
	system ("perl $dirname/get_genomic_genes_from_exons.pl $chem/$chem"."tblastn_parsed_list_nogff_filtered.txt $chem/$chem\_genomic_exon_proteins.fasta $chem/$chem $maxintron");

	# Filtering putative erroneus proteins with HMMER

	print "Doing $chem hmmer in newly identified genomic genes\n";
	system ("hmmsearch -o $chem\/hmmer\/$name\_genomic_genes\_$chem\.out --notextw --tblout $chem\/hmmer\/$name\_genomic_genes_$chem\.tblout --domtblout $chem\/hmmer\/$name\_genomic_genes_$chem\.domtblout $chemdir\/$chem\_db\.hmm $chem/$chem\_genomic_genes_proteins.fasta");
	system ("perl $dirname/get_genomic_hmmer_parsed.pl $chem\/hmmer\/$name\_genomic_genes_$chem\.domtblout $chem/$chem\_genomic_genes $evalue");
	system ("perl $dirname/get_fasta_fromalist.pl $chem/$chem\_genomic_genes_proteins.fasta $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem\_genomic_genes_hmmerparsed");
	system ("perl $dirname/get_genomic_fasta_cut.pl $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem\_genomic_genes_hmmerparsed_proteins.fasta $chem/$chem\_genomic_genes_hmmerparsed");


	open (File, "<", "$chem/$chem\_genomic_genes_hmmerparsed_proteins_cut.fasta");
	my %fasta_genomic;
	my $fastaname;
	my $genomic_count = 0;
	while (<File>) {
		chomp;
		my $line = $_;
		next if ($line !~ /\S+/);
		if ($line =~ /^>(\S+)/){
			$fastaname = $1;
			$genomic_count++;
		}
		else {
			$fasta_genomic{$fastaname} .= $line;
		}
	}
	close File;

	system("cp $chem/$chem\_genomic_genes_hmmerparsed_proteins_cut.fasta $chem/$chem\_genomic_proteins_cut.fasta");


	# Creating a final file with non-redundant sequences (i.e. isoforms, sequences from duplicated scaffolds...)

	open (NRfasta, ">", "$chem/$chem\_genomic_proteins_cut_nr.fasta");
	open (File, "<", "$chem/$chem\_genomic_proteins_cut.fasta");
	my $totalnrseqs = "0";
	my @sequences = ();
	while (<File>) {
		chomp;
		my $line = $_;
		if ($line =~ /^>(\S+)/){
			$fastaname = $1;
		}
		else {
			my $match = "0";
			foreach my $individualseq (@sequences){
				if ($individualseq =~ /$line/){
					$match++;
				}
				elsif ($line =~ /$individualseq/){
					$match++;
				}
			}
			next if ($match > 0);
			print NRfasta ">$fastaname\n$line\n";
			$totalnrseqs++;
			push (@sequences, $line);
		}
	}
	close File;

	print Chemcounts"$chem\t$genomic_count\t$totalnrseqs\n";


	## Generating a GFF for genomic genes

	system("perl $dirname/get_genomic_gff.pl $chem/$chem\_genomic_genes_hmmerparsed_proteins_cut.fasta $chem/$chem"."tblastn_parsed_list_genomic_positions_nogff_filtered.txt $chem/$chem $genome"); # Getting GFF3 from genomic sequences, although it is very raw and should be edited after manually filtering, or via Apollo
	system ("perl $dirname/get_genomic_gff_filtered_cut.pl $chem/$chem\_genomic_genes_unfiltered.gff3 $genome $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem");

	# Validating the obtained GFF3

	system ("blastp -query $chem/$chem\_genomic_genes_hmmerparsed_proteins_cut.fasta -subject $chem/$chem"."gffgenomiccut.pep.fasta -out $chem\/$chem\_genomic_protsVsGFF\_blastp\.outfmt6 -evalue $evalue -num_threads $threads -outfmt \"6 std qlen slen\"");
	system ("perl $dirname/confirm_GFF_proteins.pl $chem/$chem\_genomic_genes_hmmerparsed_proteins_cut.fasta $chem\/$chem\_genomic_protsVsGFF\_blastp\.outfmt6 $chem\/$chem\_genomic");


	# Creating a final GFF for all annotations
	system("perl $dirname/get_nr_gff_only1gff.pl $chem/$chem\_genomic_proteins_cut_nr.fasta $chem/$chem\_genomic_genes_cut.gff3 $chem/$chem");
	system("mv $chem/$chem\_genomic_and_annotated_genes_nr.gff3 $chem/$chem\_genomic_genes_cut_nr.gff3");

	print "----------------- $chem DONE\n\n";

}
close Chemcounts;

