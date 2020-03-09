#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
my $dirname = dirname(__FILE__);

# Script to check if input data is OK and parseable by BITACORA

# usage: perl Scripts/check_data_genome_mode.pl $genome query_directory



my ($line, $name, $nameout);
my $genome = $ARGV[0];
my $querydir = $ARGV[1];
my $gemoma = $ARGV[2];

print "\n----------------- Checking input data and prerequisites\n";


#Checking Fasta

my $checkgenome = "0";
open (File , "<", $genome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		$checkgenome++;
	}	
}
close File;

if ($checkgenome == 0){
	die "ERROR in $dirname/check_data.pl: Genome file is not found or empty! $genome\n";
}


print "Genome fasta parsed correctly\n";


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

# Check if GeMoMa is installed (if needed)

if ($gemoma =~ /T/){
	my $gpath = $ARGV[3];
	system ("java -jar $gpath CLI GeMoMa info > testGemoma.out 2> testGemoma.err");
	my $gemerr = 0;
	#open (File, "<", "testGemoma.err");
	#while(<File>){
	#	chomp;
	#	$line = $_;
	#	next if ($line !~ /\S+/);	
	#	if ($line =~ /Unable/){
	#		$gemerr++;
	#	}
	#}
	#close File;
	my $gemok = 0;
	open (File, "<", "testGemoma.out");
	while(<File>){
		chomp;
		$line = $_;
		next if ($line !~ /\S+/);	
		if ($line =~ /GeMoMa/){
			$gemok++;
		}
	}
	close File;

	if ($gemerr == 0 && $gemok > 0){
		#OK
	} else {
		die "ERROR in $dirname/check_data.pl: GeMoMa could not be found\nAre you sure you set the path to GeMoMa jar file correctly and java is installed in your system?\n";
	}

	print "GeMoMa installed correctly\n";
}


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

system("rm testHMMER.err testHMMER.out testBLAST.err testBLAST.out testhmm.err testhmm.out testdb.err testdb.out testGemoma.err testGemoma.out");



