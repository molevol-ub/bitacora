#!/usr/bin/perl
use strict;
use warnings;

# usage: perl Scripts/get_annot_genes_gff_v2.pl $gff $chem/$chem\_combinedsearches_list.txt $chem/$chem

# V2: If there are hits in non-overlapping positions in the protein (i.e. two fused genes by the structural annotation software), it splits the protein (flag _split)



my ($line, $name, $nameout);
my $genome = $ARGV[1];

# Read GFF
my %gffcds; my %gffgene;
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
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){  # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
			$genename = $1;

			## Extra for ensembl GFF3
			if ($genename =~ /transcript\:(\S+)/){
				$genename = $1;
			}
			if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
				$genename = "$1"."-P$2";
			}

		}
		else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: It does not recognize Parent ID in the GFF3 in: $line\n";}
		if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
			push ( @{$gffcds{$genename}{$subline[3]}}, $line);
		}
		 else {
		 	push ( @{$gffcds{$genename}{$subline[4]}}, join ('\t', $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]));
		}
		
	}
	elsif ($subline[2] =~ /mRNA/ || $subline[2] =~ /transcript/){
		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){ # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
			$genename = $1;

			## Extra for ensembl GFF3
			if ($genename =~ /transcript\:(\S+)/){
				$genename = $1;
			}
			if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
				$genename = "$1"."-P$2";
			}
		}
		else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: It does not recognize ID in the GFF3 in: $line\n";}

		if (exists $gffgene{$genename}){ # Control for duplicated genes in GFF3
			die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Gene $genename is duplicated in the GFF3, found duplicate in $line\nPlease, take a look into your GFF3 and delete duplicated genes\n";
		}

		if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
			$gffgene{$genename} = $line;
		}
		 else {
		 	$gffgene{$genename} = join ('\t', $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]);
		}

	}
}
close GFFfile;


#Opening identified proteins and writting results


open (Results, ">", "$ARGV[3]\_annot_genes.gff3");
open (File , "<", $ARGV[2]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);

	my $gene = "";
	my $geneparent = "";
	if ($subline[0] =~ /(\S+)\_split\d/){
		$gene = $1;
	} else {
		$gene = $subline[0];
	}

	if ($gene =~ /(\S+)\./){
		$geneparent = $1;
	} else {
	#	die "No encuentro gene parent en $gene en $line\n";
		$geneparent = $gene;
	}

	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});
		print Results "$subl2[0]\tAnnotGFF\tmRNA\t$subl2[3]\t$subl2[4]\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$subline[0];Parent=$geneparent;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";

##old
##		my $n = 1;
##		foreach my $cds (@{$gffcds{$gene}}) {
##			my @subl3 = split (/\t/, $cds);
##			print Results "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$subline[0]\.CDS$n;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";		
##			$n++;
#		}
###
###new - ordering CDS
		my $n = 1;

		my @positions = keys ($gffcds{$gene});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: No forward/reverse in $gffgene{$gene}\n";}

		foreach my $posit (@posorder){
			foreach my $cds (@{$gffcds{$gene}{$posit}}) {
				my @subl3 = split (/\t/, $cds);
				print Results "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$subline[0]\.CDS$n;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";		
				$n++;
			}
		}
###



	} else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Protein gene $subline[0] is not found in the GFF3\n";}
	
}
close File;
print Results "END	END	mRNA	x	x	.	.	.	X";
close Results;

# It prints now the genes cut according to hitted positions

open (Results, ">", "$ARGV[3]\_annot_genes_cut.gff3");

open (File , "<", $ARGV[2]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);

	my $gene = "";
	my $geneparent = "";
	if ($subline[0] =~ /(\S+)\_split\d/){
		$gene = $1;
	} else {
		$gene = $subline[0];
	}

	if ($gene =~ /(\S+)\./){
		$geneparent = $1;
	} else {
	#	die "No encuentro gene parent en $gene en $line\n";
		$geneparent = $gene;
	}

	my $printcds;
	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});

		my $inigene = 0; my $endgene = 0; my $endgenereached = 0; my $inicuted = 0; my $endcuted = 0;
		my $prelength = 0; my $precdslength = 0;
		my $ncds = 1;

		my $ifix = 0;


		my @positions = keys ($gffcds{$gene});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: No forward/reverse in $gffgene{$gene}\n";}



		################ FORWARD STRAND

		if ($subl2[6] =~ /\+/){
		foreach my $posit (@posorder){
		foreach my $cds (@{$gffcds{$gene}{$posit}}) {
			my @subl3 = split (/\t/, $cds);
			my $cdsid = "";
			if ($subl3[8] =~ /ID=([^;]+)/){
				$cdsid = $1;
			} else {
				$cdsid = "$subline[0]"."$ncds";
				$ncds++;
			#	die "No encuentro nombre correcto cds";
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
				$inicuted = $subl3[3];
			}


			if ($inicuted >= $subl3[3] && $endcuted <= $subl3[4]){
				$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$inicuted\t$endcuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				$endgenereached++;
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
				$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$inicuted\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n"; # Pruebo a ver si sacando ese cds completo queda bien: NADA, sigue quedando mal la proteína, y algunas son extra largas con otros dominios en el mismo CDS anotado así que nada
				$prelength += (int(($subl3[4] - $inicuted)/3))*3;
			}
			#elsif ($inicuted < $subl3[3]){
			#	die "$inicuted menor que $subl3[3] en  get_annot_genes_gff.pl\n$line\n$cds\n";
			#}
			else {
				die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Debugging control: Incongruence observed in (some situation not considered):\n$line\n$cds\n";
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
			if ($subl3[8] =~ /ID=(\S+);Parent/){
				$cdsid = $1;
			} else {
				$cdsid = "$subline[0]"."$ncds";
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
			}


			if ($inicuted > $subl3[4]){ # en el caso de restar demasiado
				$inicuted = $subl3[4] - $subl3[7] ;
			}


			if ($inicuted <= $subl3[4] && $endcuted >= $subl3[3]){
				$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$endcuted\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				$endgenereached++;
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
				$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n"; # Pruebo a ver si sacando ese cds completo queda bien: NADA, sigue quedando mal la proteína, y algunas son extra largas con otros dominios en el mismo CDS anotado así que nada
				$prelength += (int(($inicuted - $subl3[3])/3))*3;
			}
			#elsif ($inicuted < $subl3[3]){
			#	die "$inicuted menor que $subl3[3] en  get_annot_genes_gff.pl\n$line\n$cds\n";
			#}
			else {
				die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Debugging control: Incongruence observed in (some situation not considered):\n$line\n$cds\ninicuted:$inicuted\nendcuted:$endcuted\n";
			}

		}
		}
		} else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Cant't find strand in GFF3: $gffgene{$gene}\n";}


		print Results "$subl2[0]\tAnnotGFF\tmRNA\t$inigene\t$endgene\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$printcds";

	} else {die "ERROR in Scripts/get_annot_genes_gff_v2.pl: Protein gene $subline[0] is not found in the GFF3\n\n";}
	
}
close File;
print Results "END	END	mRNA	x	x	.	.	.	X";
close Results;


#Codify proteins and CDS from the generated GFFs
system ("perl Scripts/gff2fasta_v2.pl $genome $ARGV[3]\_annot_genes_cut.gff3 $ARGV[3]gffcut");
system ("perl Scripts/reconvert_fasta.pl $ARGV[3]gffcut\.pep.fasta");
#system ("perl gff2fasta_v2.pl $genome $ARGV[3]\_annot_genes.gff3 $ARGV[3]"."gff");




