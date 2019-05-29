#!/usr/bin/perl
use strict;
use warnings;

# Script to reformat a NCBI GFF to be parsable by BITACORA. It searchs for the protein name and renames all genes, mRNA and CDS in the GFF

# Usage: perl reformat_ncbi_gff.pl gfffile

my $inputgff = $ARGV[0];

my %gffcds; my %gffgene;
my @genes;
my %protid;
my $dups = "0";
my $dupsid = "";
my $keep = "";

my $errprotid = "0";

open (GFFfile , "<", $inputgff) 
	or die "ERROR in reformat_ncbi_gff.pl:Could not open file '$inputgff' $!";
while (<GFFfile>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
		#if ($subline[8] =~ /(transcript_id|Transcript) \"([^"]+)\"/){  # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
			$genename = $1;
				## Extra for ensembl GFF3
			if ($genename =~ /ranscript\:(\S+)/){
				$genename = $1;
			}
			#if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
			#	$genename = "$1"."-P$2"; ## Change P for R if protein is named as transcript
			#}

		} elsif ($subline[8] =~ /(transcript_id|Transcript)."([^"]+)"/){  # Reading GTF formats.
			$genename = $2;

		}
		else {die "ERROR in reformat_ncbi_gff.pl: It does not recognize Parent ID in the GFF3 in: $line\n";}


		if ($subline[6] =~ /\+/){ # Forward strand
			#OK
		} elsif ($subline[6] =~ /\-/){
				#OK
		} else {die "ERROR in reformat_ncbi_gff.pl: No forward/reverse in: $line\n";}


		if ($subline[7] =~ /\d+/){ # Forward strand
			#OK
		} else {die "ERROR in reformat_ncbi_gff.pl: No frame in CDS in: $line\n";}

		## Protein ID from NCBI
		my $proteinid = "none";
		if ($subline[8] =~ /rotein_id=([^;]+)/){
			$proteinid = $1;
		}
		else {
			print "WARNING in reformat_ncbi_gff.pl: It does not find ProteinID (protein_id=) in the GFF3 in:\n$line\nAre you sure the GFF is from NCBI? Expected format is as:\nNW_018367575.1  Gnomon  CDS     3420    3965    .       -       0       ID=cds0;Parent=rna0;Dbxref=GeneID:107443680,Genbank:XP_015913113.1;Name=XP_015913113.1;gbkey=CDS;gene=LOC107443680;product=protein maternal effect lethal 26-like;protein_id=XP_015913113.1\nStepping this protein...\n";
			$errprotid = 0;
		}

		$protid{$genename}=$proteinid;

		if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
			push ( @{$gffcds{$genename}{$subline[3]}}, $line);	
		}
		 else {
		 	push ( @{$gffcds{$genename}{$subline[4]}}, join ("\t", $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]));
		}
			
	}
#	elsif ($subline[2] =~ /mRNA/ || $subline[2] =~ /transcript/){
	elsif ($subline[2] =~ /mRNA/){	# NCBI uses transcript different than mRNA
		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /(transcript_id|Transcript) \"([^"]+)\"/){  # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
		$genename = $1;

		## Extra for ensembl GFF3
		if ($genename =~ /ranscript\:(\S+)/){
			$genename = $1;
			}
			#if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
			#	$genename = "$1"."-P$2"; ## Change P for R if protein is named as transcript
			#}
		} elsif ($subline[8] =~ /(transcript_id|Transcript)."([^"]+)"/){  # Reading GTF formats.
			$genename = $2;

		}
		else {die "ERROR in reformat_ncbi_gff.pl: It does not recognize ID in the GFF3 in: $line\n";}

		if (exists $gffgene{$genename}){ # Control for duplicated genes in GFF3
			die "ERROR in reformat_ncbi_gff.pl: Gene $genename is duplicated in the GFF3, found duplicate in $line\nPlease, take a look into your GFF3 and delete duplicated genes\n";
		}


		if ($subline[6] =~ /\+/){ # Forward strand
			#OK
		} elsif ($subline[6] =~ /\-/){
			#OK
		} else {die "ERROR in reformat_ncbi_gff.pl: No forward/reverse in: $line\n";}


		if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
			$gffgene{$genename} = $line;
		}
		 else {
		 	$gffgene{$genename} = join ("\t", $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]);
		}

		#if ($keep =~ /$genename\_\_/){
		#	$dups++;
		#	$dupsid .= "$genename ";
		#} else {
		#	$keep .= "$genename\_\_";
			push (@genes, $genename);
		#}

	}
}
close GFFfile;


open (Results, ">", "$inputgff\_reformatted.gff3");
foreach my $gene (@genes){
	my @subl = split (/\t/, $gffgene{$gene});
	if ($subl[8] =~ /Parent=[^;]+;(.*);/){
		my $rest = $1;
		#my $geneid = $protid{$gene};
		my $geneid = "";
		if (exists $protid{$gene}){
			$geneid = $protid{$gene};
			next if ($geneid =~ /none/);
		} else {
			print "ERROR in reformat_ncbi_gff.pl: a protein ID cannot be found for gene $gene   Please exclude this gene from GFF\n";
			next;
		}

		if ($keep =~ /$geneid  /){
			$dups++;
			$dupsid .= "$geneid ";
			next;
		} else {
			$keep .= "$geneid  ";
		}

		my $parentid = "";
		if ($geneid =~ /(\S+)\.\d+$/){
			$parentid = $1;
		} else {
			$parentid = $geneid;
		}
		print Results "$subl[0]\t$subl[1]\tgene\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\t$subl[7]\tID=$geneid;$rest\n";
		print Results "$subl[0]\t$subl[1]\tmRNA\t$subl[3]\t$subl[4]\t$subl[5]\t$subl[6]\t$subl[7]\tID=$geneid;Parent=$parentid;$rest\n";

		my @positions = keys (%{$gffcds{$gene}});
		my @posorder;
		my $n = 1;

		if ($subl[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in reformat_ncbi_gff.pl: No forward/reverse in $gffgene{$gene}\n";}

		my $firstcds = 1;
		foreach my $posit (@posorder){
			my $cdsnum = @{$gffcds{$gene}{$posit}};
			foreach my $cds (@{$gffcds{$gene}{$posit}}) {
				my @subl3 = split (/\t/, $cds);
				if ($subl3[8] =~ /Parent=[^;]+;(.*)/){
					my $rest2 = $1;
					my $chain = $subl3[6];

					if ($chain =~ /\+/ && $firstcds == 1){
						my $newpos1 = $subl3[3] + $subl3[7];
						print Results "$subl3[0]\t$subl3[1]\t$subl3[2]\t$newpos1\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$geneid\.CDS$n;Parent=$geneid;$rest2\n";	
					} 

					elsif ($chain =~ /\-/ && $firstcds == 1){
						my $newpos1 = $subl3[4] - $subl3[7];
						print Results "$subl3[0]\t$subl3[1]\t$subl3[2]\t$subl3[3]\t$newpos1\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$geneid\.CDS$n;Parent=$geneid;$rest2\n";	
					} 

					else {
						print Results "$subl3[0]\t$subl3[1]\t$subl3[2]\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$geneid\.CDS$n;Parent=$geneid;$rest2\n";		
					}
					
					$firstcds++;
					$n++;
				}
				else {die "ERROR in reformat_ncbi_gff.pl: No Parent ID found in $cds\n$subl3[8]\n";}
			}
		}

	} else {
		die "ERROR in reformat_ncbi_gff.pl: It does not recognize Parent in the GFF3 in: $gffgene{$gene}\n";
	}


}
close Results;

if ($dups > 0){
	print "Found $dups duplicated genes in the GFF\n$dupsid\n";
}

if ($errprotid > 0){
	print "\nFound $errprotid protein in the GFF3 without protein_id that were ignored\n";
}


