#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use Readgff qw (readgff);
my $dirname = dirname(__FILE__);

# Script to check if input data is OK and parseable by PAIP

# usage: perl Scripts/check_data.pl $gff $genome $proteins query_directory



my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $genome = $ARGV[1];
my $proteome = $ARGV[2];
my $querydir = $ARGV[3];

print "\n----------------- Checking input data and prerequisites\n";

# Checking GFF

my ($gffgeneref, $gffcdsref, $gffscafcdsref) = &readgff($gff);
my %gffcds = %$gffcdsref; 
my %gffgene = %$gffgeneref;
my %gffscafcds = %$gffscafcdsref;


print "GFF parsed correctly\n";


#Checking if IDs from Fastas and GFF match

my $checkgenome = "0";
my $scafnames = "";
open (File , "<", $genome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		my $gene = $1;
		$scafnames .= "$gene\_\_";
		$checkgenome++;
	}	
}
close File;

if ($checkgenome == 0){
	die "ERROR in $dirname/check_data.pl: Genome file is not found or empty! $genome\n";
}

my $checkprot = "0";
open (File , "<", $proteome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ />(\S+)/){
		my $gene = $1;
		$checkprot++;
		if (exists $gffgene{$gene}){
			#OK
			my @subline = split(/\t/, $gffgene{$gene});
			if ($scafnames =~ /$subline[0]\_\_/){
				#OK
			} else {
				die "ERROR in $dirname/check_data.pl: Genome sequence $subline[0] from GFF is not found in the genome fasta file\n";
			}

			if (exists $gffcds{$gene}){
				#OK
			} else {
				die "ERROR in $dirname/check_data.pl: Gene $gene does not contain any CDS sequence in GFF file. Please check your GFF\n";
			}


		} elsif ($gene =~ /(\S+)\-P(\S+)/) {
			my $genenew = "$1"."-R$2"; ## Change P for R if protein is named as transcript

			if (exists $gffgene{$genenew}){
				#OK
				my @subline = split(/\t/, $gffgene{$genenew});
				if ($scafnames =~ /$subline[0]\_\_/){
					#OK
				} else {
					die "ERROR in $dirname/check_data.pl: Genome sequence $subline[0] from GFF is not found in the genome fasta file\n";
				}

				if (exists $gffcds{$genenew}){
					#OK
				} else {
					die "ERROR in $dirname/check_data.pl: Gene $gene does not contain any CDS sequence in GFF file. Please check your GFF\n";
				}
			}
			else {
				die "ERROR in $dirname/check_data.pl: Protein gene $gene nor $genenew is not found in the GFF3\nYou can run the following script to check which proteins are not annotated in the GFF, or contain different name in mRNA|transcript field in the GFF\nperl $dirname/Tools/get_proteins_notfound_ingff.pl $gff $proteome\nIf your GFF comes from NCBI, you propably need to reformat to assign protein ID in mRNA and CDS fields. For that, you can use the following Script and input the formatted GFF in PAIP\nperl $dirname/Tools/reformat_ncbi_gff.pl $gff\n";
			}

		} else {
			die "ERROR in $dirname/check_data.pl: Protein gene $gene is not found in the GFF3\nYou can run the following script to check which proteins are not annotated in the GFF, or contain different name in mRNA|transcript field as ID\nperl $dirname/Tools/get_proteins_notfound_ingff.pl $gff $proteome\n";
		}
	}	
}
close File;

if ($checkprot == 0){
	die "ERROR in $dirname/check_data.pl: Protein file is not found or empty! $proteome\n";
}

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
	die "ERROR in $dirname/check_data.pl: hmmsearch could not be found\nAre you sure you set the path to HMMER bin correctly and it is properly installed?\n";
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
	die "ERROR in $dirname/check_data.pl: BLAST binaries could not be found\nAre you sure you set the path to BLAST bin correctly and it is properly installed?\n";
}

print "BLAST installed correctly\n";

# Check if the query dir contains proper renamed fasta and HMM files

my $nfile = 0;
my $dbfiles = "";
system("ls $querydir\/*_db.fasta > Query_genes.txt");
open(File, "<", "Query_genes.txt");
while(<File>){
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my $id = "";
	if ($line =~ /([^\/]+)\_db.fasta/){
		$id = $1;
	} else {
		die "Are you sure you renamed you Query DB to QUERY_db.fasta? Avoid using special characters in QUERY as\nCannot find QUERY_db.fasta in $line\n";
	}

	system ("wc -l $line > testdb.out 2> testdb.err");

	my $dbout = 0;
	open (File2, "<", "testdb.out");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line !~ /\S+/);	
		if ($line2 =~ /(\d+) /){
			$dbout = $1;
		} else {
			die "ERROR in $dirname/check_data.pl: reading testdb.out in $line2\n";
		}
	}
	close File2;

	if ($dbout == 0){
		die "ERROR in $dirname/check_data.pl: $line is empty!\n";		
	}

	# Checking now HMM file

	system ("wc -l $querydir\/$id\_db.hmm > testhmm.out 2> testhmm.err");

	$hmmerr = 0;
	open (File2, "<", "testhmm.err");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		$hmmerr++;
	}
	close File2;
	$hmmok = 0;
	open (File2, "<", "testhmm.out");
	while(<File2>){
		chomp;
		my $line2 = $_;
		next if ($line2 !~ /\S+/);	
		#my @subline = split (/\s/, $line);
		#$hmmok = $subline[0];
		if ($line2 =~ /(\d+) /){
			$hmmok = $1;
		} else {
			die "ERROR in $dirname/check_data.pl: reading testhmm.out in $line2\n";
		}

	}
	close File2;

	if ($hmmerr == 0 && $hmmok > 0){
		#OK
	} elsif ($hmmerr == 0 && $hmmok == 0){
		die "ERROR in $dirname/check_data.pl: $querydir\/$id\_db.hmm is empty!\n";		
	} else {
		die "ERROR in $dirname/check_data.pl: It was not found the expected HMM file: $querydir\/$id\_db.hmm\n";
	}

	$nfile++;
	$dbfiles .= "$id ";

}
close File;

if ($nfile > 0){
	#OK
} else {
	die "ERROR in $dirname/check_data.pl: It was not found any file formatted correctly in $querydir\nFiles must be named as ID_db.fasta and ID_db.hmm as detailed in manual\n";
}

print "Query directory files found and named correctly in $querydir\nFound $nfile query databases to search: $dbfiles\n";

print "Everything looks fine\n----------------- DONE\n\n";

system("rm testHMMER.err testHMMER.out testBLAST.err testBLAST.out testhmm.err testhmm.out testdb.err testdb.out");



