#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# perl $dirname/get_gemoma_gff_filternovel_v2.pl $chem\/gemoma_outdir/$chem\_gemoma_genes_novel.gff3 $chem\/gemoma_outdir/$chem\_overlap_CDS_novelgenes.txt $chem\/gemoma_outdir/$chem $chem $genome
# usage: perl get_gemoma_gff.pl $chem\/gemoma_outdir/final_annotation.gff $chem\/gemoma_outdir/comparison.tabular $chem\/gemoma_outdir/ $chem $genome

my ($line, $name, $nameout);

my $ingff = $ARGV[0];
my $intable = $ARGV[1];
my $outdir = $ARGV[2];
my $chem = $ARGV[3];
my $genome = $ARGV[4];

my $novelgenes = "";


# Reading GFF
my %gffcds; my %gffgene; my %gffgenel;
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
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1; 
		}
		else {die "ERROR in get_gemoma_gff_filternovel_v2.pl: It fails detecting Parent ID in $line\n";}


#		if ($genename =~ /(\S+)\_\d+DOM/){ ## DOM in upper case instead of dom (because of line 34)
#			$genename = $1;
#		} 

		if ($subline[4] > $subline[3]){ # Control to be sure that positions are ordered
			push ( @{$gffcds{$genename}{$subline[3]}}, $line);
		}
		 else {
		 	push ( @{$gffcds{$genename}{$subline[4]}}, join ("\t", $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]));
		}
		
	}
	elsif ($subline[2] =~ /mRNA/){

		next if ($line =~ /^END/);

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "ERROR in get_gemoma_gff_filternovel_v2.pl: It fails detecting ID in $line\n";}

#		if ($genename =~ /(\S+)\_\d+DOM/){ ## DOM in upper case instead of dom (because of line 34)
#			$genename = $1;
#		} 

		if (exists $gffgene{$genename}){ # Control for duplicated genes in GFF3
#			die "ERROR in get_gemoma_gff_filternovel_v2.pl: Gene $genename is duplicated in the GFF3, found duplicate in $line\nPlease, take a look into your GFF3 and delete duplicated genes\n";
		}

		if ($subline[4] > $subline[3]){ # Control to be sure that positions are ordered
			$gffgene{$genename} = $line;
		}
		 else {
		 	$gffgene{$genename} = join ("\t", $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]);
		}

	}

}
close GFFfile;


# Read table

#open (Resultsa, ">", "$outdir\_gemoma_genes_all.gff3");
open (Results, ">", "$outdir\_gemoma_genes_novel_v2.gff3");

my $genecount = "000";

open (File , "<", $intable); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subl = split (/\t/, $line);

	#my $gene = uc($subl[3]); ### uc added to avoid gemoma upper case or not in different versions
	my $gene = $subl[0];
	#my $geneparent = "gene_$gene";

	#$genecount++;
	#my $ngene = "$chem"."g$genecount\.t1";
	#my $ngeneparent = "$chem"."g$genecount";

	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});

		my $parentid = "";
		my $desc = "";
		if ($subl2[8] =~ /ID=[^;]+;Parent=([^;]+);(\S+)/){
			$desc = $2;
			$parentid = $1;
		} else {
			die "ERROR in get_gemoma_gff_filternovel_v2.pl: Could not find description after ID in $gffgene{$gene}\n";
		}

#		print Resultsa "$subl2[0]\tGeMoMa\tgene\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$ngeneparent;$desc\n";
#		print Resultsa "$subl2[0]\tGeMoMa\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$ngene;Parent=$ngeneparent;$desc\n";

		if ($subl[9] > 4){
#			print Results "$subl2[0]\tGeMoMa\tgene\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$ngeneparent;$desc\n";
#			print Results "$subl2[0]\tGeMoMa\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$ngene;Parent=$ngeneparent;$desc\n";

			##print Results "$subl2[0]\tGeMoMa\tgene\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$parentid;$desc\n";
			##print Results "$subl2[0]\tGeMoMa\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\t$subl2[8]\n";

			delete $gffgene{$gene};

		}

		my $n = 1;

		my @positions = keys (%{$gffcds{$gene}});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "No forward/reverse in $gffgene{$gene}\n";}

		foreach my $posit (@posorder){
			foreach my $cds (@{$gffcds{$gene}{$posit}}) {
				my @subl3 = split (/\t/, $cds);
				#print Resultsa "$subl3[0]\tGeMoMa\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$ngene\.CDS$n;Parent=$ngene;\n";		
				if ($subl[9] > 4){
					#print Results "$subl3[0]\tGeMoMa\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$ngene\.CDS$n;Parent=$ngene;\n";
					##print Results "$subl3[0]\tGeMoMa\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\t$subl3[8]\n";
				}
				$n++;
			}
		}

	} else {
		#die "ERROR in get_gemoma_gff_filternovel_v2.pl: It fails finding the gene in GFF $gene\n";
		# It is ok, because a gene can appear multiple times and therefore already deleted from the hash
	}


#	if ($subl[9] < 4){
#		my $ngene = $subl[0];
#		$ngene =~ s/\.//g;
#		$novelgenes .= "$ngene ";
#	}
}
close File;


foreach my $gene (sort keys %gffgene) {

	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});

		my $parentid = "";
		my $desc = "";
		if ($subl2[8] =~ /ID=[^;]+;Parent=([^;]+);(\S+)/){
			$desc = $2;
			$parentid = $1;
		} else {
			die "ERROR in get_gemoma_gff_filternovel_v2.pl: Could not find description after ID in $gffgene{$gene}\n";
		}

		print Results "$subl2[0]\tGeMoMa\tgene\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$parentid;$desc\n";
		print Results "$subl2[0]\tGeMoMa\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\t$subl2[8]\n";


		my $n = 1;

		my @positions = keys (%{$gffcds{$gene}});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "No forward/reverse in $gffgene{$gene}\n";}

		foreach my $posit (@posorder){
			foreach my $cds (@{$gffcds{$gene}{$posit}}) {
				my @subl3 = split (/\t/, $cds);
				print Results "$subl3[0]\tGeMoMa\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\t$subl3[8]\n";
				$n++;
			}
		}

	} else {die "ERROR in get_gemoma_gff_filternovel_v2.pl: It fails finding the gene in GFF $gene\n";}

}
close File;


close Results;
#close Resultsa;



#Codify proteins and CDS from the generated GFFs

#system ("perl $dirname/gff2fasta_v3.pl $genome $outdir\_gemoma_genes_all.gff3 $outdir\_gemoma_genes_all");
system ("perl $dirname/gff2fasta_v3.pl $genome $outdir\_gemoma_genes_novel_v2.gff3 $outdir\_gemoma_genes_novel_v2");




