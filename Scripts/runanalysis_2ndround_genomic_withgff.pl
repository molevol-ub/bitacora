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
my $transcripts = $ARGV[1];
my $chemdir = $ARGV[2];
my $genome = $ARGV[3];
my $translation = "$ARGV[3].frame";
my $gfffile = $ARGV[4];
my $evalue = $ARGV[5];
my $maxintron = $ARGV[6];
my $threads = $ARGV[7];

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
	if ($line =~ /$genome\.nin/ || $line =~ /$genome\.nhr/ || $line =~ /$genome\.nsq/ || $line =~ /$genome\.\d+\.nsq/){
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

if ($indexed > 1){ # Changed from > 3, for big genomes there are multiple files
	print "Found a blast index for $genome\n";
} else {
	system ("makeblastdb -dbtype nucl -in $genome");
}

#


open (Chemcounts, ">", "$name\_genecounts_annotated_and_genomic_proteins.txt");
print Chemcounts"Gene\/Gene Family\tNumber of annotated genes Identified\tAnnotated genes with tblastn hits in genomic regions\tAnnotated genes with tblastn hits in genomic regions identified in protein search\tAnnotated genes with tblastn hits in genomic regions NOT identified in protein search\tNumber of putative not annotated genes\tNumber of not annotated genes with high identity with identified annotated genes (removed as putative assembly artifacts)\tNumber of not annotated genes with high identity with NOT identified annotated genes (removed as putative assembly artifacts)\tTotal number of identified genes (Annotated + Genomic)\tTotal number of identified genes clustering identical sequences\n";

open (Chemcounts2, ">", "$name\_genecounts_summary.txt");
print Chemcounts2"Gene\/Gene Family\tNumber of annotated genes Identified\tNumber of putative not annotated genes\tTotal number of identified genes (Annotated + Genomic)\tTotal number of identified genes clustering identical sequences\n";


foreach my $chem (@chemosensory){
	print "\n----------------- Starting $chem Genomic Gene Identification and Annotation\n";
	system ("cat $chemdir\/$chem\_db\.fasta $chem\/$chem\_proteins_trimmed.fasta > $chem\/$chem\_db_masannot.fasta");

	# run tblastn
	print "Doing $chem tblastn\n";
	#system ("tblastn -query $chem\/$chem\_db_masannot.fasta -db $genome -out $chem\/$name\_Vs$chem\_tblastn\.outfmt6 -outfmt \"6 std sframe qlen slen\" -evalue $evalue -num_threads $threads "); # Deprecated
	system ("perl $dirname/run_parallel_tblastn.pl $chem\/$chem\_db_masannot.fasta $genome $threads $evalue $chem\/$name\_Vs$chem\_tblastn\.outfmt6");

	# Parsing tblastn output
	print "Parsing $chem tblastn\n";

	#system ("perl $dirname/get_genomic_tblastn_parsed_newv_nogff.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $evalue");
	#system ("perl $dirname/get_genomic_tblastn_parsed_newv_gff.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $gfffile $evalue");
	system ("perl $dirname/get_genomic_tblastn_parsed_newv_gff_genomic_positions.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $gfffile $evalue");
	system ("perl $dirname/get_genomic_tblastn_parsed_newv_nogff_genomic_positions.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $evalue");


	###### Exiting if no novel genes are found

	my $nonew = 0;
	open (Filetblastn, "<", "$chem\/$chem"."tblastn_parsed_list.txt");
	while (<Filetblastn>){
		chomp;
		my $nlin = $_;
		next if ($nlin !~ /\S+/);
		$nonew++;
	}
	close Filetblastn;

	if ($nonew == 0){
		print "No $chem identified in unannotated genomic regions\n";

		system ("cp $chem/$chem\_proteins_trimmed.fasta $chem/$chem\_genomic_and_annotated_proteins_trimmed.fasta");
		open (NRfasta, ">", "$chem/$chem\_genomic_and_annotated_proteins_trimmed_nr.fasta");
		open (File, "<", "$chem/$chem\_genomic_and_annotated_proteins_trimmed.fasta");
		my $totalnrseqs = "0";
		my $totalseqs = "0";
		my %fasta_annot;
		my $fastaname;
		my @sequences = ();
		while (<File>) {
			chomp;
			my $line = $_;
			if ($line =~ /^>(\S+)/){
				$fastaname = $1;
				$totalseqs++;
			}
			else {
				$fasta_annot{$fastaname} .= $line;
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

		open (File, "<", "$chem/$chem"."tblastn_genescatched_list.txt");
		my $genescatchedgff = "0";
		my $geneannotcatched = "0";
		my $genenotannotcatched = "0";
		while (<File>) {
			chomp;
			my $line = $_;
			if ($line =~ /Totalgenescatched (\d+)/){
				$genescatchedgff = $1;
			}
			else {
				my @subline = split (" ", $line);
				if (exists $fasta_annot{$subline[0]}){
					$geneannotcatched++;
				}
				else {
					$genenotannotcatched++;
				}
			}
		}
		close File;

		print Chemcounts"$chem\t$totalseqs\t$genescatchedgff\t$geneannotcatched\t$genenotannotcatched\t0\t0\t0\t$totalseqs\t$totalnrseqs\n";
		print Chemcounts2 "$chem\t$totalseqs\t0\t$totalseqs\t$totalnrseqs\n";

		# Creating a final GFF for all annotations
		system("perl $dirname/get_nr_gff_only1gff.pl $chem/$chem\_genomic_and_annotated_proteins_trimmed_nr.fasta $chem/$chem\_annot_genes_trimmed.gff3 $chem/$chem");
		system("cat $chem/$chem\_annot_genes_trimmed.gff3 |  sed '/^END\tEND\tmRNA/d' > $chem/$chem\_genomic_and_annotated_genes.gff3");


		print "----------------- $chem DONE\n\n";
		next;
	}

	#######

	### Script continues
	# Obtaining protein sequences, exons and genes

	system ("perl $dirname/get_genomic_fasta_frameseqs_fromalist.pl $translation $chem/$chem"."tblastn_parsed_list.txt $chem/$chem");
	system ("perl $dirname/get_genomic_genes_from_exons.pl $chem/$chem"."tblastn_parsed_list.txt $chem/$chem\_genomic_exon_proteins.fasta $chem/$chem $maxintron");

	# Filtering putative erroneus proteins with HMMER

	print "Doing $chem hmmer in newly identified genomic genes\n";
	system ("hmmsearch -o $chem\/hmmer\/$name\_genomic_genes\_$chem\.out --notextw --tblout $chem\/hmmer\/$name\_genomic_genes_$chem\.tblout --domtblout $chem\/hmmer\/$name\_genomic_genes_$chem\.domtblout $chemdir\/$chem\_db\.hmm $chem/$chem\_genomic_genes_proteins.fasta");
	system ("perl $dirname/get_genomic_hmmer_parsed.pl $chem\/hmmer\/$name\_genomic_genes_$chem\.domtblout $chem/$chem\_genomic_genes $evalue");
	system ("perl $dirname/get_fasta_fromalist.pl $chem/$chem\_genomic_genes_proteins.fasta $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem\_genomic_genes_hmmerparsed");
	system ("perl $dirname/get_genomic_fasta_trimmed.pl $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem\_genomic_genes_hmmerparsed_proteins.fasta $chem/$chem\_genomic_genes_hmmerparsed");

	# Removing putative new sequences already annotated (i.e. duplicated regions due to assembly artifacts)

	#print "Doing $chem blast genomic vs anotated $chem\n";
	system ("blastp -query $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta -subject $chem/$chem\_proteins\_trimmed.fasta -evalue 1e-5 -max_target_seqs 1 -outfmt 6 -out $chem/$chem\_genomicVsanotated.outfmt6");

	open (File, "<", "$chem/$chem\_proteins\_trimmed.fasta");
	my %fasta_annot;
	my $fastaname;
	while (<File>) {
		chomp;
		my $line = $_;
		next if ($line !~ /\S+/);
		if ($line =~ /^>(\S+)/){
			$fastaname = $1;
		}
		else {
			$fasta_annot{$fastaname} .= $line;
		}
	}
	close File;

	open (File, "<", "$chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta");
	my %fasta_genomic;
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

	open (Blast, "<", "$chem/$chem\_genomicVsanotated.outfmt6");
	my %blast_relation;
	while (<Blast>){
		chomp;
		my $line = $_;
		my @subline = split (/\t/, $line);
		next if ($subline[2] < 95); # ID at least of 90% to consider it as a duplicate artifact
		if (exists $blast_relation{$subline[1]}){
			next if ($blast_relation{$subline[1]} =~ /$subline[0] /);
		}
		$blast_relation{$subline[1]}.="$subline[0] ";
	}
	close Blast;

	my $annotingenomic = 0;
	my $totalannotgenomic = 0;
	open (Results, ">", "$chem/$chem\_genomic_and_annotated_proteins_trimmed.fasta");
	foreach my $key (sort keys %fasta_annot){
		if (exists $blast_relation{$key}){
			print Results ">$key $blast_relation{$key}\n$fasta_annot{$key}\n";	

			my @subblast = split (/\s/, $blast_relation{$key});
			foreach my $deletion (@subblast){
				delete $fasta_genomic{$deletion};
				$annotingenomic++;
			}
		}
		else {
			print Results ">$key\n$fasta_annot{$key}\n";	
		}
	$totalannotgenomic++;
	}

	# Checking now with annotated proteins not annotated in phase 1
	open (ResultsNotannot, ">", "$chem/$chem\_genomic_not_$chem\_annotated\_proteins_trimmed.fasta");
	foreach my $key (sort keys %fasta_genomic){
		print ResultsNotannot ">$key\n$fasta_genomic{$key}\n";
	}
	close ResultsNotannot;
	#print "Doing $chem blast genomic proteins vs all anotated proteins\n";
	system ("blastp -query $chem/$chem\_genomic_not_$chem\_annotated\_proteins_trimmed.fasta -db $transcripts -evalue 1e-5 -max_target_seqs 1 -outfmt 6 -num_threads $threads -out $chem/$chem\_genomicnotAnnotVsAllproteome.outfmt6");
	#system ("perl $dirname/run_parallel_blastp.pl $chem/$chem\_genomic_not_$chem\_annotated\_proteins_trimmed.fasta $transcripts $threads 1e-5 $chem/$chem\_genomicnotAnnotVsAllproteome.outfmt6");
	open (Blast, "<", "$chem/$chem\_genomicnotAnnotVsAllproteome.outfmt6");
	my $deletedgenomic = 0;
	while (<Blast>){
		chomp;
		my $line = $_;
		my @subline = split (/\t/, $line);
		next if ($subline[2] < 95);
		if (exists $fasta_genomic{$subline[0]}){
			delete $fasta_genomic{$subline[0]};
			$deletedgenomic++;
		}
	}
	close Blast;

	# Printing results in final files
	foreach my $key (sort keys %fasta_genomic){
		print Results ">$key\n$fasta_genomic{$key}\n";
		$totalannotgenomic++;
	}

	close Results;

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

	open (File, "<", "$chem/$chem"."tblastn_genescatched_list.txt");
	my $genescatchedgff = "0";
	my $geneannotcatched = "0";
	my $genenotannotcatched = "0";
	while (<File>) {
		chomp;
		my $line = $_;
		if ($line =~ /Totalgenescatched (\d+)/){
			$genescatchedgff = $1;
		}
		else {
			my @subline = split (" ", $line);
			if (exists $fasta_annot{$subline[0]}){
				$geneannotcatched++;
			} elsif ($subline[0] =~ /(\S+)\-P(\S+)/){
				my $genenew = "$1"."-R$2"; ## Change P for R if protein is named as transcript
				if (exists $fasta_annot{$genenew}){
					$geneannotcatched++;
				} else {
					$genenotannotcatched++;
				}
			} else {
				$genenotannotcatched++;
			}
		}
	}
	close File;

	# Creating a final file with non-redundant sequences (i.e. isoforms, sequences from duplicated scaffolds...)

	open (NRfasta, ">", "$chem/$chem\_genomic_and_annotated_proteins_trimmed_nr.fasta");
	open (File, "<", "$chem/$chem\_genomic_and_annotated_proteins_trimmed.fasta");
	my $totalnrseqs = "0";
	my %nrfa;
	while (<File>) {
		chomp;
		my $line = $_;
		if ($line =~ /^>(\S+)/){
			$fastaname = $1;
		}
		else {
			my $match = "0";
			my $hit = "";
			my $hitl = "";
			foreach my $kseq (keys %nrfa){
				my $individualseq = $nrfa{$kseq};
				if ($individualseq =~ /$line/){
					$match++;
					$hit = $kseq;
					$hitl = length($individualseq);
				}
				elsif ($line =~ /$individualseq/){
					$match++;
					$hit = $kseq;
					$hitl = length($individualseq);
				}
			}

			if ($match > 0){
				my $seql = length($line);
				if ($seql <= $hitl){
					next;
				} elsif ($seql > $hitl){
					$nrfa{$fastaname} = $line;
					delete $nrfa{$hit};
				}
			} else {
				$totalnrseqs++;
				$nrfa{$fastaname} = $line;

			}

		}
	}
	close File;

	foreach my $kseq (sort keys %nrfa){
		print NRfasta ">$kseq\n$nrfa{$kseq}\n";
	}
	close NRfasta;
	

	my $finalgenomic = $genomic_count - $annotingenomic - $deletedgenomic;

	if ($count == 0){
		print Chemcounts "$chem\t0\t$genescatchedgff\t$geneannotcatched\t$genenotannotcatched\t$genomic_count\t$annotingenomic\t$deletedgenomic\t$totalannotgenomic\t$totalnrseqs\n";
		print Chemcounts2 "$chem\t0\t$finalgenomic\t$totalannotgenomic\t$totalnrseqs\n";
	}
	else {
		my $medialength = int($cutlength/$count);
		my $medialengthprot = int ($protlength/$count);
		print Chemcounts"$chem\t$count\t$genescatchedgff\t$geneannotcatched\t$genenotannotcatched\t$genomic_count\t$annotingenomic\t$deletedgenomic\t$totalannotgenomic\t$totalnrseqs\n";
		print Chemcounts2 "$chem\t$count\t$finalgenomic\t$totalannotgenomic\t$totalnrseqs\n";
	}

	## Generating a GFF for genomic genes

	system("perl $dirname/get_genomic_gff.pl $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta $chem/$chem"."tblastn_parsed_list_genomic_positions.txt $chem/$chem $genome"); # Getting GFF3 from genomic sequences, although it is very raw and should be edited after manually filtering, or via Apollo
	system ("perl $dirname/get_genomic_gff_filtered_trimmed.pl $chem/$chem\_genomic_genes_unfiltered.gff3 $genome $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem");

	# Validating the obtained GFF3

	system ("blastp -query $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta -subject $chem/$chem"."gffgenomictrimmed.pep.fasta -out $chem\/$chem\_genomic_protsVsGFF\_blastp\.outfmt6 -evalue $evalue -outfmt \"6 std qlen slen\"");
	system ("perl $dirname/confirm_GFF_proteins_withlength.pl $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta $chem\/$chem\_genomic_protsVsGFF\_blastp\.outfmt6 $chem\/$chem\_genomic");


	# Creating a final GFF for all annotations
	system("perl $dirname/get_nr_gff.pl $chem/$chem\_genomic_and_annotated_proteins_trimmed_nr.fasta $chem/$chem\_annot_genes_trimmed.gff3 $chem/$chem\_genomic_genes_trimmed.gff3 $chem/$chem\_genomic_and_annotated_genes_nr.gff3");
	system("perl $dirname/get_nr_gff.pl $chem/$chem\_genomic_and_annotated_proteins_trimmed.fasta $chem/$chem\_annot_genes_trimmed.gff3 $chem/$chem\_genomic_genes_trimmed.gff3 $chem/$chem\_genomic_and_annotated_genes.gff3");
	#system("cat $chem/$chem\_annot_genes_trimmed.gff3 $chem/$chem\_genomic_genes_trimmed.gff3 |  sed '/^END\tEND\tmRNA/d' > $chem/$chem\_genomic_and_annotated_genes.gff3");

	print "----------------- $chem DONE\n\n";

}
close Chemcounts;

