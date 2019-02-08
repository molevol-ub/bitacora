#!/usr/bin/perl
use strict;
use warnings;
#use File::Basename;
#use lib dirname (__FILE__);

use FindBin;
use lib "$FindBin::Bin/../";

use Readgff qw (readgff);

# Script to find proteins not annotated in the GFF file and separate from those annotated

# usage: perl Scripts/get_proteins_notfound_ingff.pl $gff $proteins


my ($line, $name, $nameout);
my $gff = $ARGV[0];
my $proteome = $ARGV[1];


print "\nReading GFF\n";

# Checking GFF

my ($gffgeneref, $gffcdsref, $gffscafcdsref) = &readgff($gff);
my %gffcds = %$gffcdsref; 
my %gffgene = %$gffgeneref;
my %gffscafcds = %$gffscafcdsref;


print "GFF parsed correctly\n";


#Checking proteins included in GFF

open (Results , ">", "$proteome\_ingff.fasta"); 
open (Resultsb , ">", "$proteome\_nogff.fasta"); 
my $ok = "0";
my $notok = "0";
my (@okgene, @notokgene);

my %fasta;
open (File , "<", $proteome); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);

	if ($line =~ />(\S+)/){
		my $gene = $1;
		$name = $gene;

		if (exists $gffgene{$gene}){
			#OK
			$ok++;
			push (@okgene, $gene);

			if (exists $gffcds{$gene}){
				#OK
			} else {
				print "WARNING: Gene $gene does not contain any CDS sequence in GFF file. Please check your GFF\nThis Warning will cause ERROR in PAIP\n";
			}


		} elsif ($gene =~ /(\S+)\-P(\S+)/) {
			my $genenew = "$1"."-R$2";
			if (exists $gffgene{$genenew}){
				#OK
				$ok++;
				push (@okgene, $gene);

				if (exists $gffcds{$genenew}){
					#OK
				} else {
					print "WARNING: Gene $gene does not contain any CDS sequence in GFF file. Please check your GFF\nThis Warning will cause ERROR in PAIP\n";
				}
			} else {
				$notok++;
				push (@notokgene, $gene);				
			}

		} else {
			$notok++;
			push (@notokgene, $gene);
		}
	}
	else {
		$fasta{$name} .= $line;
	}	
}
close File;

foreach my $gene (@okgene){
	print Results ">$gene\n$fasta{$gene}\n";
}
foreach my $gene (@notokgene){
	print Resultsb ">$gene\n$fasta{$gene}\n";
}


close Results;
close Resultsb;

my $tot = $ok + $notok;

print "Identified $tot protein sequences in $protome\n";
print "Saved $ok proteins with GFF annotation in $proteome\_ingff.fasta\n";
print "Saved $notok proteins without GFF annotation in $proteome\_nogff.fasta\n";
print "Done\n\n";

