#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

#usage: perl get_genomic_gff.pl X_genomic_genes_hmmerparsed_proteins_cuted.fasta Xtblastn_parsed_list_genomic_positions.txt nameout Genome_fasta


my $genome = $ARGV[3];

my ($line, $name, $id);
my (%gffprot, %gffgen, %gffcuted);
my @annogenes;
my @genomegenes;

open (File, "<", $ARGV[0]);
while(<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/); #$name =~ /g/
	if ($line =~ /^>(\S+)/){
		$name = $1;
		push (@genomegenes, $name);		
	}
}

open (Results, ">", "$ARGV[2]\_genomic_genes_unfiltered.gff3");


open (File , "<", $ARGV[1]); #Genome positions genomic genes
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/ /, $line);
	$gffgen{$subline[0]} = $line;
}
close File;

foreach my $gene (@genomegenes){
	if ($gene =~ /(\S+?)_(exon\S+)/){ 
		my $scaf = $1;
		my $exons = $2;
		my @exon = split (/\_/, $exons);

		my $ini = 9999999999999999999999999999999999990;
		my $fin = 0;
		my $chain = "+";

		my @cdsline;
		my $initialcds = 0;
		foreach my $ex (@exon){
			my $gen = "$scaf\_$ex";
			next if ($ex =~ /dom/); # Ignore information about domains, not an exon to read
			if (exists $gffgen{$gen}) {
				my @col = split (/ /, $gffgen{$gen});
#lg1_ord2_scaf1805	AnnotGFF	mRNA	458567	471355	.	.	.	ID=g114.t1;blastphmmer;complete;Pos:80-477
#lg1_ord2_scaf1805	AnnotGFF	CDS	458567	459228	.	.	.	Parent=g114.t1;blastphmmer;complete;Pos:80-477
				#print Results "$scaf\tGenomicGFF\tCDS\t$col[3]\t$col[4]\t.\t.\t.\tParent=$gene;\n";
				#print Results2 "$scaf\tGenomicGFF\tCDS\t$col[3]\t$col[4]\t.\t.\t.\tParent=$gene;\n";

				if ($col[1] =~ /Frame4/ || $col[1] =~ /Frame5/ || $col[1] =~ /Frame6/){
					$chain = "-";
				}
				my $frame = "";
				if ($col[1] =~ /Frame1/ || $col[1] =~ /Frame4/){
					$frame = "0";
				} elsif ($col[1] =~ /Frame2/ || $col[1] =~ /Frame5/){
					$frame = "1";
				} elsif ($col[1] =~ /Frame3/ || $col[1] =~ /Frame6/){
					$frame = "2";
				}

				my $inicds = "";
				#if ($initialcds == 0){
				#			$inicds = $col[3] - $frame;
				#			$initialcds++;		
				#} else {
					$inicds = $col[3];
				#}


				push (@cdsline, join("\t", $scaf,"GenomicGFF","CDS",$inicds,$col[4],".",$chain,$frame,"ID=$gen;Parent=$gene;"));
				#push (@cdsline, join("\t", $scaf,"GenomicGFF","CDS",$inicds,$col[4],".",$chain,"0","ID=$gen;Parent=$gene;")); # Debugging, replacing frame to 0 to test

				if ($ini > $col[3]){
					$ini = $col[3];
				}
				if ($fin < $col[4]){
					$fin = $col[4];
				}				


			} else {die "ERROR in get_genomic_gff.pl: It does not exist the genomic gene $gene\n";}
		}

		print Results "$scaf\tGenomicGFF\tgene\t$ini\t$fin\t.\t$chain\t.\tID=$gene;\n";
		print Results "$scaf\tGenomicGFF\tmRNA\t$ini\t$fin\t.\t$chain\t.\tID=$gene;Parent=$gene;\n";
		foreach my $cds (@cdsline){
			my @cdssplit = split(/\t/, $cds);
			print Results "$cdssplit[0]\t$cdssplit[1]\t$cdssplit[2]\t$cdssplit[3]\t$cdssplit[4]\t$cdssplit[5]\t$cdssplit[6]\t$cdssplit[7]\t$cdssplit[8]\n";
		}


	} else {die "ERROR in get_genomic_gff.pl: It does not identify the genomic name in $gene\n";}

}
#print Results "END	END	mRNA	x	x	.	.	.	X\n";
close Results;

#system ("perl $dirname/gff2fasta_v2.pl $genome $ARGV[2]\_genomic_genes_unfiltered.gff3 $ARGV[2]gffgenomic");
#system ("perl $dirname/reconvert_fasta.pl $ARGV[2]gffgenomic\.pep.fasta");
system ("perl $dirname/gff2fasta_v3.pl $genome $ARGV[2]\_genomic_genes_unfiltered.gff3 $ARGV[2]gffgenomic");


