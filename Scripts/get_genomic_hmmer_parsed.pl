#!/usr/bin/perl
use strict;
use warnings;

#usage:

#  perl get_genomic_hmmer_parsed.pl hmmer.domtblout out_name evalue

my $evalue = $ARGV[2];

my $hits = "";
my ($line, $name, $id);
my (%frame, %ini, %fin, %multipledom);

my $file= "$ARGV[0]";
open (Blastfile , "<", $file);
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	$line =~ s/\s+/\t/g;$line =~ s/\t+/\t/g;
	my @subline = split (/\t/, $line);
	if ($subline[6] <= $evalue){ # E-value
		next if ($subline[12] > $evalue);
		if ($hits !~ /$subline[0]\_/) {
			$ini{$subline[0]} = $subline[17];
			$fin{$subline[0]} = $subline[18];
			$hits .= $subline[0];
			$hits .= "_";
			$multipledom{$subline[0]} = "1";
		}
		else {
			if ($subline[17] < $ini{$subline[0]} ) {
				$ini{$subline[0]} = $subline[17];
			}
			if ($subline[18] > $fin{$subline[0]} ) {
				$fin{$subline[0]} = $subline[18];
			}
			$multipledom{$subline[0]}++;
		}

	}
}
close Blastfile;



open (Results , ">", "$ARGV[1]hmmer_parsed_list.txt");
foreach my $key (sort keys %ini) {
	print Results "$key - $ini{$key} $fin{$key} hmmer $multipledom{$key}\n";	
}







