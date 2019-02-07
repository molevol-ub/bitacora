#!/usr/bin/perl
use strict;
use warnings;

# Script to check if input data is OK and parseable by PAIP

# usage: perl Scripts/check_data.pl $gff $genome $proteins query_directory



my ($line, $name, $nameout);
my $genome = $ARGV[1];
my $proteome = $ARGV[2];
my $querydir = $ARGV[3];

print "----------------- Checking input data and prerequisites\n";

# Checking GFF

my %gffcds; my %gffgene;
open (GFFfile , "<", $ARGV[0]); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){  # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
			$genename = $1;

			## Extra for ensembl GFF3
			if ($genename =~ /transcript\:(\S+)/){
				$genename = $1;
			}
			if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
				$genename = "$1"."-P$2";
			}

		}
		else {die "ERROR in Scripts/check_data.pl: It does not recognize Parent ID in the GFF3 in: $line\n";}


		if ($subline[6] =~ /\+/){ # Forward strand
			#OK
		} elsif ($subline[6] =~ /\-/){
			#OK
		} else {die "ERROR in Scripts/check_data.pl: No forward/reverse in: $line\n";}


		if ($subline[7] =~ /\d+/){ # Forward strand
			#OK
		} else {die "ERROR in Scripts/check_data.pl: No frame in CDS in: $line\n";}


		if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
			push ( @{$gffcds{$genename}{$subline[3]}}, $line);
		}
		 else {
		 	push ( @{$gffcds{$genename}{$subline[4]}}, join ('\t', $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]));
		}
		
	}
	elsif ($subline[2] =~ /mRNA/ || $subline[2] =~ /transcript/){
		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){ # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
			$genename = $1;

			## Extra for ensembl GFF3
			if ($genename =~ /transcript\:(\S+)/){
				$genename = $1;
			}
			if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
				$genename = "$1"."-P$2"; ## Change P for R if protein is named as transcript
			}
		}
		else {die "ERROR in Scripts/check_data.pl: It does not recognize ID in the GFF3 in: $line\n";}

		if (exists $gffgene{$genename}){ # Control for duplicated genes in GFF3
			die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Gene $genename is duplicated in the GFF3, found duplicate in $line\nPlease, take a look into your GFF3 and delete duplicated genes\n";
		}


		if ($subline[6] =~ /\+/){ # Forward strand
			#OK
		} elsif ($subline[6] =~ /\-/){
			#OK
		} else {die "ERROR in Scripts/check_data.pl: No forward/reverse in: $line\n";}


		if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
			$gffgene{$genename} = $line;
		}
		 else {
		 	$gffgene{$genename} = join ('\t', $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]);
		}

	}
}
close GFFfile;

print "GFF parsed correctly\n";


#Checking if IDs from Fastas and GFF match

my $scafnames = "";
open (File , "<", $genome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ />(\S+)/){
		my $gene = $1;
		$scafnames .= "$gene\_\_";
	}	
}
close File;

open (File , "<", $proteome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ />(\S+)/){
		my $gene = $1;
		if (exists $gffgene{$gene}){
			#OK
			my @subline = split(/\t/, $gffgene{$gene});
			if ($scafnames =~ /$subline[0]\_\_/){
				#OK
			} else {
				die "ERROR in Scripts/check_data.pl: Genome sequence $subline[0] is not found in the genome fasta file\n";
			}

			if (exists $gffcds{$gene}){
				#OK
			} else {
				die "ERROR in Scripts/check_data.pl: Genome sequence $subline[0] is not found in the genome fasta file\n";
			}


		} else {die "ERROR in Scripts/check_data.pl: Protein gene $gene is not found in the GFF3\n";}
	}	
}
close File;

print "Genome and protein fasta parsed correctly\n";


# Check if hmmer is installed

system ("hmmsearch -h > testHMMER.out 2> testHMMER.err");
my $hmmerr = 0;
open (File, "<", "testHMMER.err");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$hmmerr++;
}
close File;
my $hmmok = 0;
open (File, "<", "testHMMER.out");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	if ($line =~ /domtblout/){
		$hmmok++;
	}
}
close File;

if ($hmmerr == 0 && $hmmok > 0){
	#OK
} else {
	die "ERROR in Scripts/check_data.pl: hmmsearch could not be found\nAre you sure you set the path to HMMER bin correctly and it is properly installed?\n";
}

print "HMMER installed correctly\n";

# Check if blast is installed

system ("blastp -help > testBLAST.out 2> testBLAST.err");
system ("blastn -help >> testBLAST.out 2>> testBLAST.err");
system ("makeblastdb -help >> testBLAST.out 2>> testBLAST.err");
my $blasterr = 0;
open (File, "<", "testBLAST.err");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	$blasterr++;
}
close File;
my $blastok = 0;
open (File, "<", "testBLAST.out");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);	
	if ($line =~ /outfmt/){
		$blastok++;
	}
}
close File;

if ($blasterr == 0 && $blastok > 0){
	#OK
} else {
	die "ERROR in Scripts/check_data.pl: BLAST binaries could not be found\nAre you sure you set the path to BLAST bin correctly and it is properly installed?\n";
}

print "BLAST installed correctly\n";

# Check if the query dir contains proper renamed fasta and HMM files

my $nfile = 0;
my $dbfiles = "";
system("ls $querydir\/*_db.fasta > Query_genes.txt");
open(File, "<", "Query_genes.txt");
while(<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my $id = "";
	if ($line =~ /([^\/]+)\_db.fasta/){
		$id = $1;
	} else {
		die "Are you sure you renamed you Query DB to QUERY_db.fasta? Avoid using special characters in QUERY as\nCannot find QUERY_db.fasta in $line\n";
	}

	# Checking now HMM file

	system ("wc -l $querydir\/$id\_db.hmm > testhmm.out 2> testhmm.err");

	$hmmerr = 0;
	open (File2, "<", "testhmm.err");
	while(<File2>){
		chomp;
		$line = $_;
		next if ($line !~ /\S+/);	
		$hmmerr++;
	}
	close File2;
	$hmmok = 0;
	open (File2, "<", "testhmm.out");
	while(<File2>){
		chomp;
		$line = $_;
		next if ($line !~ /\S+/);	
		my @subline = split (/\s/, $line);
		$hmmok = $subline[0];

	}
	close File2;

	if ($blasterr == 0 && $blastok > 0){
		#OK
	} elsif ($blasterr == 0 && $blastok == 0){
		die "ERROR in Scripts/check_data.pl: $querydir\/$id\_db.hmm is empty!\n";		
	} else {
		die "ERROR in Scripts/check_data.pl: It was not found the expected HMM file: $querydir\/$id\_db.hmm\n";
	}

	$nfile++;
	$dbfiles .= "$id ";

}
close File;

if ($nfile > 0){
	#OK
} else {
	die "ERROR in Scripts/check_data.pl: It was not found any file formatted correctly in $querydir\nFiles must be names as ID_db.fasta and ID_db.hmm as detailed in manual\n";
}

print "Query directory files found and named correctly\nFound $nfile query databases to search: $dbfiles\n";

print "Everything looks fine\n----------------- DONE\n\n";

system("rm testHMMER.err testHMMER.out testBLAST.err testBLAST.out testhmm.err testhmm.out");