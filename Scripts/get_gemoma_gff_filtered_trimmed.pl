#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# usage: perl get_genomic_gff_filtered_trimmed.pl X_genomic_genes_unfiltered.gff3 $genome X_genomic_geneshmmer_parsed_list.txt outname


my ($line, $name, $nameout);
my $genome = $ARGV[1];

# Reading GFF
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
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: It fails detecting Parent ID in $line\n";}


		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

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
		else {print "ERROR in get_gemoma_gff_filtered_trimmed.pl: It fails detecting ID in $line\n";}

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 

		if (exists $gffgene{$genename}){ # Control for duplicated genes in GFF3
			die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Gene $genename is duplicated in the GFF3, found duplicate in $line\nPlease, take a look into your GFF3 and delete duplicated genes\n";
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


# Opening hmmer file

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

	$gene = $subline[0];

	if ($gene =~ /(\S+)\.t/){
		$geneparent = "$1";
	} else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Cannot detect gene parent in $subline[0] in $line\n";}


	my $printcds;
	if (exists $gffgene{$gene}){
		my @subl2 = split (/\t/, $gffgene{$gene});

		my $inigene = 0; my $endgene = 0; my $endgenereached = 0; my $inicuted = 0; my $endcuted = 0;
		my $prelength = 0; my $precdslength = 0;

		my $ifix = 0;


		my @positions = keys (%{$gffcds{$gene}});
		my @posorder;

		if ($subl2[6] =~ /\+/){ # Forward strand
			@posorder = sort { $a <=> $b } @positions;
		} elsif ($subl2[6] =~ /\-/){
			@posorder = sort { $b <=> $a } @positions;
		} else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: No forward/reverse in $gffgene{$gene}\n";}



		################ FORWARD STRAND

		if ($subl2[6] =~ /\+/){
		foreach my $posit (@posorder){
		foreach my $cds (@{$gffcds{$gene}{$posit}}) {
			my @subl3 = split (/\t/, $cds);
			my $cdsid = "";
			if ($subl3[8] =~ /ID=([^;]+)/){
				$cdsid = $1;
			} else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Can't find CDS name in $cds\n";}

			$ifix += 2; 
			my $nifix = int($ifix/3)*3;
			
			if ($inicuted == 0){
				$inicuted = $subl3[3] + ($subline[2]*3) - $precdslength - $nifix; # $subl3[7] para el frame que codificaría coincida
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
				$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$inicuted\t$endcuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
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
				$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$inicuted\t$subl3[4]\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n"; # Pruebo a ver si sacando ese cds completo queda bien: NADA, sigue quedando mal la proteína, y algunas son extra largas con otros dominios en el mismo CDS anotado así que nada
				$prelength += (int(($subl3[4] - $inicuted)/3))*3;
			}
			#elsif ($inicuted < $subl3[3]){
			#	die "$inicuted menor que $subl3[3] en  get_annot_genes_gff.pl\n$line\n$cds\n";
			#}
			else {
				die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Debugging control: Incongruence observed in (some situation not considered):\n$line\n$cds\ninicuted:$inicuted\nendcuted:$endcuted\n";
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
			} else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Can't find CDs name in $cds\n";}

			# Hacer bien gff con reverse:
			# restar pos. inicio recortada = $subl3[4] - $subline[2]
			# pos. final recortada = $subl3[4] - $subline[3]

			$ifix += 2; 
			my $nifix = int($ifix/3)*3;
			
			if ($inicuted == 0){
				$inicuted = $subl3[4] - ($subline[2]*3) + $precdslength + $nifix +6; # $subl3[7] para el frame que codificaría coincida
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
				#$inicuted = $subl3[4] - $subl3[7] ; # salen mal las que empiezan el hmm en pos:1-X
				$inicuted = $subl3[4];
			}


			if ($inicuted <= $subl3[4] && $endcuted >= $subl3[3]){
				$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$endcuted\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
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
				$printcds .= "$subl3[0]\tGenomicGFF\tCDS\t$subl3[3]\t$inicuted\t$subl3[5]\t$subl3[6]\t$subl3[7]\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				#$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$subl3[3]\t$subl3[4]\t$subl3[5]\t$subl3[6]\t0\tID=$cdsid;Parent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n"; # Pruebo a ver si sacando ese cds completo queda bien: NADA, sigue quedando mal la proteína, y algunas son extra largas con otros dominios en el mismo CDS anotado así que nada
				$prelength += (int(($inicuted - $subl3[3])/3))*3;
			}
			#elsif ($inicuted < $subl3[3]){
			#	die "$inicuted menor que $subl3[3] en  get_annot_genes_gff.pl\n$line\n$cds\n";
			#}
			else {
				die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Debugging control: Incongruence observed in (some situation not considered):\n$line\n$cds\ninicuted:$inicuted\nendcuted:$endcuted\n";
			}

		}
		}
		} else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Cant't find strand in GFF3: $gffgene{$gene}\n";}

		my $desc = "";
		if ($subl2[8] =~ /Parent=[^;]+;(\S+)/){
			$desc = $1;
		} else {
			die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Cannot find description after parent in $line\n";
		}


		print Results "$subl2[0]\tGenomicGFF\tgene\t$inigene\t$endgene\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$geneparent;$desc\;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$subl2[0]\tGenomicGFF\tmRNA\t$inigene\t$endgene\t$subl2[5]\t$subl2[6]\t$subl2[7]\tID=$subline[0];Parent=$geneparent;$desc\;$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$printcds";

	} else {die "ERROR in get_gemoma_gff_filtered_trimmed.pl: Protein gene $gene $subline[0] is not found in the GFF3\n";}
	
}
close File;
#print Results "END	END	mRNA	x	x	.	.	.	X\n";
close Results;


#Codify proteins and CDS from the generated GFFs

system ("perl $dirname/gff2fasta_v3.pl $genome $ARGV[3]\_genomic_genes_trimmed.gff3 $ARGV[3]gffgenomictrimmed");



