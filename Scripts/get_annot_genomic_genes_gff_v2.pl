#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use Readgff qw (readgff);
my $dirname = dirname(__FILE__);

# usage: perl Scripts/get_annot_genes_gff_v2.pl $gff $chem/$chem\_combinedsearches_list.txt $chem/$chem

# V2: If there are hits in non-overlapping positions in the protein (i.e. two fused genes by the structural annotation software), it splits the protein (flag _split)


my ($line, $name, $nameout);
my $genome = $ARGV[1];
my $gff = $ARGV[0];

# Read GFF
my ($gffgeneref, $gffcdsref, $gffscafcdsref) = &readgff($gff);
my %gffcds = %$gffcdsref; 
my %gffgene = %$gffgeneref;
my %gffscafcds = %$gffscafcdsref;


#Opening identified proteins and writting results


open (Results, ">", "$ARGV[3]\_genomic_genes.gff3");
open (File , "<", $ARGV[2]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);

	my $gene = "";
	my $geneparent = "";
	if ($subline[0] =~ /(\S+)(\_split\d+)/){
		$gene = $1;
		$geneparent = $gene;
	} else {
		$gene = $subline[0];
	}

	unless ($subline[0] =~ /\S+\_split\d+/){
		if ($gene =~ /(\S+)\.\d+$/){
			$geneparent = $1;
		} elsif ($gene =~ /(\S+)\-R\S$/){
			$geneparent = $1;
		} elsif ($gene =~ /(\S+)\-P\S$/){
			$geneparent = $1;
		} elsif ($gene =~ /(\S+)\.t\d+$/){
			$geneparent = $1;
		} else {
		#	die "cannot find gene parent in $gene in $line\n";
			$geneparent = "gene_$gene";
		}
	}

	my $genera = "";
	if ($gene =~ /(\S+)\-P(\S+)/) {
		$genera = "$1"."-R$2";
	}

	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});
		print Results "$subl2[0]\tGenomicGFF\tgene\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$subl2[0]\tGenomicGFF\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$subline[0];Parent=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";

##old
##		my $n = 1;
##		foreach my $cds (@{$gffcds{$gene}}) {
##			my @subl3 = split (/\t/, $cds);
##			print Results "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$subline[0]\.CDS$n;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";		
##			$n++;
##		}
###
###new - ordering CDS
		my $n = 1;

		my @positions = keys (%{$gffcds{$gene}});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: No forward/reverse in $gffgene{$gene}\n";}

		foreach my $posit (@posorder){
			foreach my $cds (@{$gffcds{$gene}{$posit}}) {
				my @subl3 = split (/\t/, $cds);
				print Results "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$subline[0]\.CDS$n;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";		
				$n++;
			}
		}
###

	} elsif (exists $gffgene{$genera}){
		my @subl2 = split (/\t/, $gffgene{$genera});
		print Results "$subl2[0]\tGenomicGFF\tgene\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$subl2[0]\tGenomicGFF\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$subline[0];Parent=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";

		my $n = 1;

		my @positions = keys (%{$gffcds{$genera}});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: No forward/reverse in $gffgene{$gene}\n";}

		foreach my $posit (@posorder){
			foreach my $cds (@{$gffcds{$gene}{$posit}}) {
				my @subl3 = split (/\t/, $cds);
				print Results "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$subline[0]\.CDS$n;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";		
				$n++;
			}
		}

	} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: Protein gene $subline[0] is not found in the GFF3\n";}
	
}
close File;
#print Results "END	END	mRNA	x	x	.	.	.	X\n";
close Results;

# It prints now the genes trimmed according to hit positions

open (Results, ">", "$ARGV[3]\_genomic_genes_trimmed.gff3");

open (File , "<", $ARGV[2]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);

	my $gene = "";
	my $geneparent = "";
	if ($subline[0] =~ /(\S+)(\_split\d+)/){
		$gene = $1;
		$geneparent = $gene;
	} else {
		$gene = $subline[0];
	}

	unless ($subline[0] =~ /\S+\_split\d+/){
		if ($gene =~ /(\S+)\.\d+$/){
			$geneparent = $1;
		} elsif ($gene =~ /(\S+)\-R\S$/){
			$geneparent = $1;
		} elsif ($gene =~ /(\S+)\-P\S$/){
			$geneparent = $1;
		} elsif ($gene =~ /(\S+)\.t\d+$/){
			$geneparent = $1;
		} else {
		#	die "cannot find gene parent in $gene in $line\n";
			$geneparent = "gene_$gene";
		}
	}

	my $genera = "";
	if ($gene =~ /(\S+)\-P(\S+)/) {
		$genera = "$1"."-R$2";
	}


	## Determine if the protein name is -PA or -RA
	if (exists $gffgene{$gene}){
		$gene = $gene;
	} elsif (exists $gffgene{$genera}){
		$gene = $genera;
	} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: Protein gene $gene $subline[0] is not found in the GFF3\n";}


	my $printcds;
	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});

		my $inigene = 0; my $endgene = 0; my $endgenereached = 0; my $inicuted = 0; my $endcuted = 0;
		my $prelength = 0; my $precdslength = 0;
		my $firstcds = $subline[2]; ## To avoid modifying first cds if it starts from first position
		my $ncds = 1;

		my $ifix = 0;


		my @positions = keys (%{$gffcds{$gene}});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: No forward/reverse in $gffgene{$gene}\n";}



		################ FORWARD STRAND

		if ($subl2[6] =~ /\+/){
		foreach my $posit (@posorder){
		foreach my $cds (@{$gffcds{$gene}{$posit}}) {
			my @subl3 = split (/\t/, $cds);
			my $cdsid = "";
			if ($subl3[8] =~ /ID=([^;]+)/){
				$cdsid = $1;
			} else {
				$cdsid = "$subline[0]"."CDS$ncds";
				$ncds++;
			#	die "Cannot find proper name for cds\n";
			}

			$ifix += 2; 
			my $nifix = int($ifix/3)*3;
			
			if ($inicuted == 0){
				$inicuted = $subl3[3] + ($subline[2]*3) - $precdslength - $nifix + $subl3[7] -3; # $subl3[7] para el frame que codificaría coincida
			} else {
				$inicuted = $subl3[3];
			}

			if ($endcuted == 0){
				$endcuted = $subl3[3] + ($subline[3]*3) - $prelength - $precdslength - $nifix +3;
			} elsif ($endgenereached == 0){
				my $diflength = $subline[3] - $subline[2]; 
				$endcuted = $subl3[3] + ($diflength*3) - $prelength;
			} elsif ($endgenereached > 0){
				next; ## Finalizo el bucle
			}

			if ($inigene == 0){
				$inigene = $inicuted;
			} 
			if ($endgenereached == 0){
				$endgene = $endcuted + $subl3[7];
			}


			if ($inicuted < $subl3[3]){ # en el caso de restar demasiado
				$inicuted = $subl3[3] + $subl3[7];
			}


			if ($inicuted >= $subl3[3] && $endcuted <= $subl3[4]){
				$endgenereached++;
				next if ($endcuted <= $inicuted); ## check gene structure
				$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$inicuted\t$endcuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#if ($inigene == 0 || $inigene )
			} 
			elsif ($endcuted <= $subl3[3]){
				next;
			}
			elsif ($inicuted >= $subl3[4]){
				$inigene = 0;
				$inicuted = 0;
				$endcuted = 0;
				$precdslength += (int(($subl3[4] - $subl3[3])/3))*3;
				next;
			}
			elsif ($inicuted >= $subl3[3] && $endcuted >= $subl3[4]){
				if ($firstcds == 1){ ## not to modify first exon starting in alignment from first position
					$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
					$firstcds++;
				} else {
					$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$inicuted\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				}
				###$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$inicuted\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n"; # Pruebo a ver si sacando ese cds completo queda bien: NADA, sigue quedando mal la proteína, y algunas son extra largas con otros dominios en el mismo CDS anotado así que nada
				$prelength += (int(($subl3[4] - $inicuted)/3))*3;
			}
			#elsif ($inicuted < $subl3[3]){
			#	die "$inicuted menor que $subl3[3] en  get_annot_genes_gff.pl\n$line\n$cds\n";
			#}
			else {
				die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: Debugging control: Incongruence observed in (some situation not considered):\n$line\n$cds\n";
			}

		}
		}
		}


		################ REVERSE STRAND

		elsif ($subl2[6] =~ /\-/){
		foreach my $posit (@posorder){
		foreach my $cds (@{$gffcds{$gene}{$posit}}) {
			my @subl3 = split (/\t/, $cds);
			my $cdsid = "";
			if ($subl3[8] =~ /ID=([^;]+)/){
				$cdsid = $1;
			} else {
				$cdsid = "$subline[0]"."CDS$ncds";
				$ncds++;
			#	die "No encuentro nombre correcto cds";
			}


			$ifix += 2; 
			my $nifix = int($ifix/3)*3;
			
			if ($inicuted == 0){
				$inicuted = $subl3[4] - ($subline[2]*3) + $precdslength + $nifix - $subl3[7] +6; # $subl3[7] para el frame que codificaría coincida
			} else {
				$inicuted = $subl3[4];
			}

			if ($endcuted == 0){
				$endcuted = $subl3[4] - ($subline[3]*3) + $prelength + $precdslength + $nifix;
			} elsif ($endgenereached == 0){
				my $diflength = $subline[3] - $subline[2]; 
				$endcuted = $subl3[4] - ($diflength*3) + $prelength;
			} elsif ($endgenereached > 0){
				next; ## Finalizo el bucle
			}

			if ($endgene == 0){
				$endgene = $inicuted;
			} 
			if ($endgenereached == 0){
				$inigene = $endcuted;
				if ($inigene < 1){ # Avoid negative values when trimming genes
					$inigene = 1;
				}
			}


			if ($inicuted > $subl3[4]){ # en el caso de restar demasiado
				$inicuted = $subl3[4] - $subl3[7] ;
			}


			if ($inicuted <= $subl3[4] && $endcuted >= $subl3[3]){
				$endgenereached++;
				next if ($endcuted >= $inicuted); ## check gene structure
				$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$endcuted\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#if ($inigene == 0 || $inigene )
			} 
			elsif ($endcuted >= $subl3[4]){
				next;
			}
			elsif ($inicuted <= $subl3[3]){
				$endgene = 0;
				$inicuted = 0;
				$endcuted = 0;
				$precdslength += (int(($subl3[4] - $subl3[3])/3))*3;
				next;
			}
			elsif ($inicuted <= $subl3[4] && $endcuted <= $subl3[3]){
				if ($firstcds == 1){ #### not to modify first exon starting in alignment from first position
					$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
					$firstcds++;
				} else {
					$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				}
				###$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n"; # Pruebo a ver si sacando ese cds completo queda bien: NADA, sigue quedando mal la proteína, y algunas son extra largas con otros dominios en el mismo CDS anotado así que nada
				$prelength += (int(($inicuted - $subl3[3])/3))*3;
			}
			#elsif ($inicuted < $subl3[3]){
			#	die "$inicuted menor que $subl3[3] en  get_annot_genes_gff.pl\n$line\n$cds\n";
			#}
			else {
				die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: Debugging control: Incongruence observed in (some situation not considered):\n$line\n$cds\ninicuted:$inicuted\nendcuted:$endcuted\n";
			}

		}
		}
		} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: Cant't find strand in GFF3: $gffgene{$gene}\n";}


		print Results "$subl2[0]\tGenomicGFF\tgene\t$inigene\t$endgene\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$subl2[0]\tGenomicGFF\tmRNA\t$inigene\t$endgene\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$subline[0];Parent=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$printcds";

	} else {die "ERROR in $dirname/get_annot_genomic_genes_gff_v2.pl: Protein gene $gene $subline[0] is not found in the GFF3\n";}
	
}
close File;
#print Results "END	END	mRNA	x	x	.	.	.	X\n";
close Results;


#Encode proteins and CDS from the generated GFFs

system ("perl $dirname/gff2fasta_v3.pl $genome $ARGV[3]\_genomic_genes_trimmed.gff3 $ARGV[3]gffgenomictrimmed");





