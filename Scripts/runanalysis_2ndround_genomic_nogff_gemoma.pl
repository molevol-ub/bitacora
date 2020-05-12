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
my $gemomap = $ARGV[6];

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

	# Filter db to avoid including the same sequence more than once
	system ("cp $chemdir\/$chem\_db.fasta $chem\/$chem\_db.fasta");
	system ("perl $dirname/get_fasta_uniqueseq.pl $chem\/$chem\_db.fasta $chem\/$chem\_db_filt.fasta");

	# run tblastn
	print "Doing $chem tblastn\n";
	#system ("tblastn -version");
	#system ("tblastn -query $chem\/$chem\_db_filt.fasta -db $genome -out $chem\/$name\_Vs$chem\_tblastn\.outfmt6 -outfmt \"6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles\" -evalue $evalue -num_threads $threads"); # Deprecated
	system ("perl $dirname/run_parallel_tblastn_gemoma.pl $chem\/$chem\_db_filt.fasta $genome $threads $evalue $chem\/$name\_Vs$chem\_tblastn\.outfmt6");

	# Parsing tblastn output
	print "Parsing $chem tblastn\n";

	# Step moved after filtering the tblastn file, it runs faster
	#system ("perl $dirname/get_genomic_tblastn_parsed_newv_nogff_genomic_positions_gemoma.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6 $chem/$chem $evalue");

	# Filter tblastn file to exclude hits in same scaffold positions that cause errors in gemoma
	system ("perl $dirname/get_tblastn_filtered_forgemoma.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6");

	# Parse filtered blast file
	system ("perl $dirname/get_genomic_tblastn_parsed_newv_nogff_genomic_positions_gemoma.pl $chem\/$name\_Vs$chem\_tblastn\.outfmt6_filtered.txt $chem/$chem $evalue");

	
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

	## Predicting new genes using GeMoMa
	print "Using GeMoMa to predict novel genes\n";

	system ("mkdir -p $chem/gemoma_outdir");

	# Detect GeMoMa version
	my $gemomaversion = "0";
	open (Gemvfile, "<", "GeMoMa_version.txt");
	while (<Gemvfile>){
		chomp;
		my $gemovline = $_;
		next if ($gemovline !~ /\S+/);
		$gemomaversion = $gemovline;
	}
	close Gemvfile;

	system ("java -jar $gemomap CLI GeMoMa s=$chem\/$name\_Vs$chem\_tblastn\.outfmt6_filtered.txt t=$genome c=$chem\/$chem\_db_filt.fasta outdir=$chem/gemoma_outdir p=10000 ct=0.2 > $chem\/gemoma.out 2> $chem\/gemoma.err");

	# Check if GeMoMa finished without errors
	my $gemoerr = "0";
	open (FileGemo, "<", "$chem\/gemoma.err");
	while (<FileGemo>){
		chomp;
		my $gemline = $_;
		next if ($gemline !~ /\S+/);
		if ($gemline =~ /Exception/){
			$gemoerr++;
		}

	}
	close FileGemo;
	if ($gemoerr > 0){
		die "ERROR in $dirname/runanalysis_2ndround_genomic_nogff_gemoma: GeMoMa died with error running GeMoMa\n";
	}

	## Filtering annotations

	system ("java -jar $gemomap CLI GAF g=$chem\/gemoma_outdir/predicted_annotation.gff outdir=$chem/gemoma_outdir f= > $chem\/gemoma.out 2> $chem\/gemoma.err");
	
	# Check if GeMoMa finished without errors
	$gemoerr = "0";
	open (FileGemo, "<", "$chem\/gemoma.err");
	while (<FileGemo>){
		chomp;
		my $gemline = $_;
		next if ($gemline !~ /\S+/);
		if ($gemline =~ /Exception/){
			$gemoerr++;
		}

	}
	close FileGemo;
	if ($gemoerr > 0){
		die "ERROR in $dirname/runanalysis_2ndround_genomic_nogff_gemoma: GeMoMa died with error runnning GAF\n";
	}


	# AnnotationFinalizer, the parameters differ between versions

	if ($gemomaversion > 163){
		system ("java -jar $gemomap CLI AnnotationFinalizer g=$genome a=$chem/gemoma_outdir/filtered_predictions.gff outdir=$chem/gemoma_outdir rename=NO > $chem\/gemoma.out 2> $chem\/gemoma.err");
	} else {
		system ("java -jar $gemomap CLI AnnotationFinalizer a=$chem/gemoma_outdir/filtered_predictions.gff outdir=$chem/gemoma_outdir rename=NO > $chem\/gemoma.out 2> $chem\/gemoma.err");	
	}
	
	# Check if GeMoMa finished without errors
	$gemoerr = "0";
	open (FileGemo, "<", "$chem\/gemoma.err");
	while (<FileGemo>){
		chomp;
		my $gemline = $_;
		next if ($gemline !~ /\S+/);
		if ($gemline =~ /Exception/){
			$gemoerr++;
		}

	}
	close FileGemo;
	if ($gemoerr > 0){
		die "ERROR in $dirname/runanalysis_2ndround_genomic_nogff_gemoma: GeMoMa died with error running AF\n";
	}


	## Extract novel annotated genes only and rename GFF

	system ("perl $dirname/get_gemoma_gff_genomemode.pl $chem\/gemoma_outdir/filtered_predictions.gff $chem\/gemoma_outdir/$chem $chem $genome");
	system ("cp $chem\/gemoma_outdir/$chem\_gemoma_genes.gff3 $chem/$chem\_genomic_genes.gff3");
	system ("cp $chem\/gemoma_outdir/$chem\_gemoma_genes.pep.fasta $chem/$chem\_genomic_genes_proteins.fasta");


	# Filtering putative erroneus proteins with HMMER

	print "Doing $chem hmmer in newly identified genomic genes\n";
	system ("hmmsearch -o $chem\/hmmer\/$name\_genomic_genes\_$chem\.out --notextw --tblout $chem\/hmmer\/$name\_genomic_genes_$chem\.tblout --domtblout $chem\/hmmer\/$name\_genomic_genes_$chem\.domtblout $chemdir\/$chem\_db\.hmm $chem/$chem\_genomic_genes_proteins.fasta");
	system ("perl $dirname/get_genomic_hmmer_parsed.pl $chem\/hmmer\/$name\_genomic_genes_$chem\.domtblout $chem/$chem\_genomic_genes $evalue");
	system ("perl $dirname/get_fasta_fromalist.pl $chem/$chem\_genomic_genes_proteins.fasta $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem\_genomic_genes_hmmerparsed");
	system ("perl $dirname/get_genomic_fasta_trimmed.pl $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem\_genomic_genes_hmmerparsed_proteins.fasta $chem/$chem\_genomic_genes_hmmerparsed");


	open (File, "<", "$chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta");
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

	system("cp $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta $chem/$chem\_genomic_proteins_trimmed.fasta");


	# Creating a final file with non-redundant sequences (i.e. isoforms, sequences from duplicated scaffolds...)

	open (NRfasta, ">", "$chem/$chem\_genomic_proteins_trimmed_nr.fasta");
	open (File, "<", "$chem/$chem\_genomic_proteins_trimmed.fasta");
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
	

	print Chemcounts"$chem\t$genomic_count\t$totalnrseqs\n";


	## Generating a GFF for genomic genes

	system ("perl $dirname/get_annot_genomic_genes_gff_v2.pl $chem/$chem\_genomic_genes.gff3 $genome $chem/$chem\_genomic_geneshmmer_parsed_list.txt $chem/$chem");


	# Validating the obtained GFF3

	system ("blastp -query $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta -subject $chem/$chem"."gffgenomictrimmed.pep.fasta -out $chem\/$chem\_genomic_protsVsGFF\_blastp\.outfmt6 -evalue $evalue -outfmt \"6 std qlen slen\"");
	system ("perl $dirname/confirm_GFF_proteins_withlength.pl $chem/$chem\_genomic_genes_hmmerparsed_proteins_trimmed.fasta $chem\/$chem\_genomic_protsVsGFF\_blastp\.outfmt6 $chem\/$chem\_genomic");


	# Creating a final GFF for all annotations
	system("perl $dirname/get_nr_gff_only1gff.pl $chem/$chem\_genomic_proteins_trimmed_nr.fasta $chem/$chem\_genomic_genes_trimmed.gff3 $chem/$chem");
	system("mv $chem/$chem\_genomic_and_annotated_genes_nr.gff3 $chem/$chem\_genomic_genes_trimmed_nr.gff3");

	print "----------------- $chem DONE\n\n";

}
close Chemcounts;

