#!/usr/bin/perl
use strict;
use warnings;

# usage: perl get_nr_gff.pl seqs_nr.fasta annnotated.gff genomic.gff outname


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
			#	die "No encuentro gene parent en $gene en $line\n";
				$geneparent = "gene_$gene";
			}
		}

		$nrgeneids .= "$geneparent\_\_";

	}
}
close Fasta;


open (Results, ">", "$ARGV[3]");

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
		else {die "ERROR in get_nr_gff.pl: It fails detecting Parent ID in $line\n";}


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
		else {print "ERROR in get_nr_gff.pl: It fails detecting ID in $line\n";}

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
		else {print "ERROR in get_nr_gff.pl: It fails detecting ID in $line\n";}

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		if ($nrgeneids =~ /$genename\__/ || $nrgeneids =~ /$genename/ ){
			print Results "$line\n";
		}

	}
}
close GFFfile;


open (GFFfile , "<", $ARGV[2]); 
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
		else {die "ERROR in get_nr_gff.pl: It fails detecting Parent ID in $line\n";}


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
		else {print "ERROR in get_nr_gff.pl: It fails detecting ID in $line\n";}

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
		else {print "ERROR in get_nr_gff.pl: It fails detecting ID in $line\n";}

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		if ($nrgeneids =~ /$genename\__/ || $nrgeneids =~ /$genename/ ){
			print Results "$line\n";
		}

	}
}
close GFFfile;




close Results;

