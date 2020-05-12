#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);
#########################################################
#####
#
# Run first analysis, using blastp and hmmer without a GFF file
#
#####

## Input variables

my $name = $ARGV[0];
my $transcripts = $ARGV[1];
my $chemdir = $ARGV[2];
my $evalue = $ARGV[3];
my $threads = $ARGV[4];

## Start

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

# Indexing proteins
print "Indexing protein file $transcripts\n";

my $indexed = 0;
system ("ls $transcripts* > index_check.txt");
open (File, "<", "index_check.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /$transcripts\.pin/ || $line =~ /$transcripts\.phr/ || $line =~ /$transcripts\.psq/){
		$indexed++;
	}


}
close File;
system("rm index_check.txt");

if ($indexed == 3){
	print "Found a blast index for $transcripts\n";
} else {
	system ("makeblastdb -dbtype prot -in $transcripts");
}


open (Chemcounts, ">", "$name\_genecounts_annotated_proteins.txt");
print Chemcounts"Gene\/Gene Family\tNumber of Genes Identified\tAverage Length of Annotated Proteins\tAverage Length of Original Raw Proteins\n";

foreach my $chem (@chemosensory){
	print "\n----------------- Starting $chem Protein Identification and Annotation\n";	
#=head
	system ("mkdir -p $chem");
	# run blastp
	print "Doing $chem blastp\n";
	system ("cp $chemdir\/$chem\_db.fasta $chem\/$chem\_db.fasta");
	#system ("blastp -query $chemdir\/$chem\_db\.fasta -db $transcripts -out $chem\/$name\_Vs$chem\_blastp\.outfmt6 -evalue $evalue -num_threads $threads -outfmt \"6 std qlen slen\" "); # Deprecated
	system ("perl $dirname/run_parallel_blastp.pl $chem\/$chem\_db\.fasta $transcripts $threads $evalue $chem\/$name\_Vs$chem\_blastp\.outfmt6");

	# run hmmer
	system ("mkdir -p $chem\/hmmer");
	print "Doing $chem hmmer\n";
	system ("hmmsearch -o $chem\/hmmer\/$name\_$chem\.out --notextw --tblout $chem\/hmmer\/$name\_$chem\.tblout --domtblout $chem\/hmmer\/$name\_$chem\.domtblout $chemdir\/$chem\_db\.hmm $transcripts");
#=cut

	# Parsing blast and hmmer outputs

	print "Parsing $chem blastp and hmmer\n";
	system ("perl $dirname/get_blastp_parsed_newv2.pl $chem\/$name\_Vs$chem\_blastp\.outfmt6 $chem/$chem $evalue");
	system ("perl $dirname/get_hmmer_parsed_newv.pl $chem\/hmmer\/$name\_$chem\.domtblout $chem/$chem $evalue");
	
	# Combining results from blast and hmmer

	system ("cat $chem/$chem"."blastp_parsed_list.txt $chem/$chem"."hmmer_parsed_list.txt > $chem/$chem\_allsearches_list.txt");
	system ("perl $dirname/get_blast_hmmer_combined.pl $chem/$chem\_allsearches_list.txt $chem/$chem");

	# Obtaining raw original and trimming protein sequences

	system ("perl $dirname/get_fasta_fromalist_v2.pl $transcripts $chem/$chem\_combinedsearches_list.txt $chem/$chem"); 
	system ("perl $dirname/get_fasta_trimmed.pl $chem/$chem\_combinedsearches_list.txt $chem/$chem\_proteins.fasta $chem/$chem");


	# Counting numbers
	open (File, "<", "$chem/$chem\_combinedsearches_list.txt");
	my $count = "0"; my $cutlength = "0"; 
	while (<File>){
		chomp;
		my $line = $_;
		my @subline = split (/\s/, $line);
		$count++;
		my $gene = $subline[0];
		my $length = $subline[3] - $subline[2];
		$cutlength += $length;
	}
	close File;

	open (File, "<", "$chem/$chem\_proteins.fasta");
	my $protlength = "0";
	while (<File>) {
		chomp;
		my $line = $_;
		unless ($line =~ /^>/){
			my $plength = length $line;
			$protlength += $plength;
		}
	}
	close File;

	if ($count == 0){
		print Chemcounts "$chem 0 0 0\n";
	}
	else {
		my $medialength = int($cutlength/$count);
		my $medialengthprot = int ($protlength/$count);
		print Chemcounts"$chem\t$count\t$medialength\t$medialengthprot\n";
	}

	print "----------------- $chem DONE\n\n";

}
close Chemcounts;



