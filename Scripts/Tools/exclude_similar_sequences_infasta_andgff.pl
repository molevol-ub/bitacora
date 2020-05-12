#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
my $dirname = dirname(__FILE__);

# Script to exlude similar sequences in a fasta file and the corresponding GFF3, retaining the longest sequence as representative (useful to cluster isoforms or proteins resulting from putative assembly artifacts)

die "\nUsage: insert the protein fasta and its GFF3 to be filtered, and the minimum length (aa) required to retain sequences. BLAST needs to be installed in the system\nperl exclude_similar_sequences_infasta_andgff.pl input.fasta input.gff3 100\n\nAdditional parameters:\nThe % identity used to filter similar sequences can be modified in line 14 of the script (default is 98)\nOutput name to be added in each sequence can be included in line 15\nThe number of threads to be used in blast search can be input in line 16 (default is 1)\n\n" unless @ARGV == 3;


my $filtlen = "$ARGV[2]"; # Length to filter sequences
my $ident = "98"; # Percent of identity to filter sequences
my $nametag = ""; # Name to add in out sequences
my $threads = "1"; # Threads to use in blastp search


my $in = $ARGV[0];

my %fasta; my $name;
my $nrids = "";
my $nrgeneids = "";

my $outname="";
if ($in =~ /(\S+)\.fa/){
	$outname = $1;
} else {
	$outname = $in;
}

open (Fasta , "<", "$in"); 
while (<Fasta>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= $line;
	}

}
close Fasta;

my %blast;

system ("makeblastdb -in $in -dbtype prot >/dev/null");
system ("blastp -query $in -db $in -out $in\.outfmt6 -evalue 1e-3 -num_threads $threads -outfmt \"6 std qlen slen\" ");

open (Blastfile , "<", "$in\.outfmt6"); 
while (<Blastfile>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > 1e-5); # Evalue filter
	push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[13], $subline[12], $subline[8], $subline[9], $subline[2]));
}
close Blastfile;

my $exclude = "";
open (Results2, ">", "$outname\_idseqsclustered_table.txt");
foreach my $key (sort keys %blast) {
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[4]*0.8);
		my $filtro2 = ($subline[3]*0.8);

		next if ($subline[2] < $filtro1 && $subline[2] < $filtro2); # Alignment length filter

		### New script

		my $escaped = $subline[0];
		$escaped =~ s/\|/\\\|/g; # Escape \|

		my $escapedk = $key;
		$escapedk =~ s/\|/\\\|/g; # Escape \|			

		next if ($key =~ /^$escaped$/);

		if ($subline[7] >= $ident){ # Identity filter
			my $length1 = length($fasta{$key});
			my $length2 = length($fasta{$subline[0]});	

			print Results2 "$key\t$subline[0]\t$subline[7]\t$subline[2]\t$length1\t$length2\n";

			if ($length1 > $length2){
				$exclude .= "$subline[0] ";
			} elsif ($length1 < $length2) {
				$exclude .= "$key ";
			} else { # same length proteins
				next if ($exclude =~ /$escapedk\s/);
				next if ($exclude =~ /$escaped\s/);
				$exclude .= "$subline[0] ";
				#print "$key\n"; # debugging
			}

		}
	}

}
close Results2;

open (Results, ">", "$outname\_idseqsclustered.fasta");
open (Results2, ">", "$outname\_idseqsclustered_renamed.fasta");
foreach my $key (sort keys %fasta) {
	my $length3 = length($fasta{$key});

	my $escapedk = $key;
	$escapedk =~ s/\|/\\\|/g; # Escape \|		
	next if ($exclude =~ /$escapedk\s/); # Exclude similar seqs

	next if ($length3 < $filtlen); # Length filter
	
	print Results ">$key\n$fasta{$key}\n";
	print Results2 ">$nametag"."$key\n$fasta{$key}\n";

	# Keep ids to filter the gff
		my $genename = $key;

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
close Results;
close Results2;



open (Results, ">", "$outname\_idseqsclustered.gff3");

open (GFFfile , "<", $ARGV[1]); 
while (<GFFfile>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {die "ERROR in exclude_similar_sequences_infasta_andgff.pl: It fails detecting Parent ID in $line\n";}


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
		else {print "ERROR in exclude_similar_sequences_infasta_andgff.pl: It fails detecting ID in $line\n";}

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
		else {print "ERROR in exclude_similar_sequences_infasta_andgff.pl: It fails detecting ID in $line\n";}

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

## Check overlapping genes
system ("perl $dirname/get_overlapping_genes_fromgff.pl $outname\_idseqsclustered.gff3");


print "\nDone, output files have been saved in $outname\_idseqsclustered.fasta, $outname\_idseqsclustered.gff3 and $outname\_idseqsclustered_table.txt\n\n";




