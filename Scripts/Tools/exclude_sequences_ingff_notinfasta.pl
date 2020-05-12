#!/usr/bin/perl
use strict;
use warnings;

# Script to exclude from a GFF all sequences not included in a fasta file (considering that protein and mRNA are named the same)

# usage: perl exclude_sequences_ingff_notinfasta.pl sequences.fasta file.gff3 outname.gff3
die "Script to exclude from a GFF all sequences not included in a fasta file\n\nperl exclude_sequences_ingff_notinfasta.pl sequences.fasta file.gff3 outname.gff3 \n\n" unless @ARGV == 3;


my ($line, $name, $nameout);
my $nrids = "";
my $nrgeneids = "";

open (Fasta , "<", $ARGV[0]); 
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/){
		my $genename = $1;

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		$nrids .= "$genename\_\_";

		my $gene = "";
		my $geneparent = "";
		if ($genename =~ /(\S+)(\_split\d+)/){
			$gene = $1;
			$geneparent = $genename;
		} else {
			$gene = $genename;
		}

		unless ($genename =~ /\S+\_split\d+/){
			if ($gene =~ /(\S+)\.\d+$/){
				$geneparent = $1;
			} elsif ($gene =~ /(\S+)\-R\S$/){
				$geneparent = $1;
			} elsif ($gene =~ /(\S+)\-P\S$/){
				$geneparent = $1;
			} elsif ($gene =~ /(\S+)\.t\d+$/){
				$geneparent = $1;
			} else {
			#	die "Can't find gene parent in $gene in $line\n";
				$geneparent = "gene_$gene";
			}
		}

		$nrgeneids .= "$geneparent\_\_";

	}
}
close Fasta;


open (Results, ">", "$ARGV[2]");

open (GFFfile , "<", $ARGV[1]); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {die "ERROR in exclude_sequences_ingff_notinfasta.pl: It fails detecting Parent ID in $line\n";}


		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		if ($nrids =~ /$genename\__/){
			print Results "$line\n";
		}

	}
	elsif ($subline[2] =~ /mRNA/){

		next if ($line =~ /^END/);

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "ERROR in exclude_sequences_ingff_notinfasta.pl: It fails detecting ID in $line\n";}

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		if ($nrids =~ /$genename\__/){
			print Results "$line\n";
		}

	}

	elsif ($subline[2] =~ /gene/){

		next if ($line =~ /^END/);

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "ERROR in exclude_sequences_ingff_notinfasta.pl: It fails detecting ID in $line\n";}

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		if ($nrgeneids =~ /$genename\__/ ){ # || $nrgeneids =~ /$genename/ 
			print Results "$line\n";
		}

	}
}
close GFFfile;

close Results;

print "\nDone, output file has been saved saved in $ARGV[2]\n\n";




