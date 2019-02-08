#!/usr/bin/perl
use strict;
use warnings;

# usage: perl Scripts/get_annot_genes_gff.pl $gff $chem/$chem\_combinedsearches_list.txt $chem/$chem

my ($line, $name, $nameout);

# Abro el gff y lo guardo en memoria
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
		if ($subline[8] =~ /Parent=(\S+);/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "falla al pillar el nombre geneid del gff\n";}
		if ($subline[3] < $subline[4]){
			push ( @{$gffcds{$genename}}, join ("\t" , $subline[0], $subline[3], $subline[4])) ;
		}
		elsif ($subline[3] > $subline[4]){
			push ( @{$gffcds{$genename}}, join ("\t" , $subline[0], $subline[4], $subline[3])) ;
		}
	}
	elsif ($subline[2] =~ /mRNA/){
		my $genename = "";
		if ($subline[8] =~ /ID=(\S+);Parent/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "falla al pillar el nombre geneid del gff\n";}
		if ($subline[3] < $subline[4]){
			$gffgene{$genename} = join ("\t" , $subline[0], $subline[3], $subline[4]) ;
		}
		elsif ($subline[3] > $subline[4]){
			$gffgene{$genename} = join ("\t" , $subline[0], $subline[4], $subline[3]) ;
		}
	}
}
close GFFfile;

#Abriendo resultados búsqueda y printando


open (Results, ">", "$ARGV[2]\_annot_genes.gff3");
open (File , "<", $ARGV[1]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);
	if (exists $gffgene{$subline[0]}){
		my @subl2 = split (/\t/, $gffgene{$subline[0]});
		print Results "$subl2[0]\tAnnotGFF\tmRNA\t$subl2[1]\t$subl2[2]\t\.\t\.\t\.\tID=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";

		foreach my $cds (@{$gffcds{$subline[0]}}) {
			my @subl3 = split (/\t/, $cds);
			print Results "$subl3[0]\tAnnotGFF\tCDS\t$subl3[1]\t$subl3[2]\t\.\t\.\t\.\tParent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";		
		}
	} else {die "No encuentra en el gff $subline[0]\n";}
	
}
close File;
close Results;

#Ahora printo las posiciones recortadas

open (Results, ">", "$ARGV[2]\_annot_genes_cuted.gff3");

open (File , "<", $ARGV[1]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);
	my $printcds;
	if (exists $gffgene{$subline[0]}){
		my @subl2 = split (/\t/, $gffgene{$subline[0]});

		my $inigene = 0; my $endgene = 0; my $endgenereached = 0; my $inicuted = 0; my $endcuted = 0;
		my $prelength = 0; my $precdslength = 0;
		foreach my $cds (@{$gffcds{$subline[0]}}) {
			my @subl3 = split (/\t/, $cds);
			
			if ($inicuted == 0){
				$inicuted = $subl3[1] + ($subline[2]*3) - 3 - $precdslength;
			} else {
				$inicuted = $subl3[1];
			}

			if ($endcuted == 0){
				$endcuted = $subl3[1] + ($subline[3]*3) - 3 - $prelength - $precdslength;
			} elsif ($endgenereached == 0){
				my $diflength = $subline[3] - $subline[2]; 
				$endcuted = $subl3[1] + ($diflength*3) - $prelength;
			}

			if ($inigene == 0){
				$inigene = $inicuted;
			} 
			if ($endgenereached == 0){
				$endgene = $endcuted;
			}


			if ($inicuted >= $subl3[1] && $endcuted <= $subl3[2]){
				$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$inicuted\t$endcuted\t\.\t\.\t\.\tParent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				$endgenereached++;
				#if ($inigene == 0 || $inigene )
			} 
			elsif ($endcuted <= $subl3[1]){
				next;
			}
			elsif ($inicuted >= $subl3[2]){
				$inigene = 0;
				$inicuted = 0;
				$endcuted = 0;
				$precdslength += ($subl3[2] - $subl3[1]);
				next;
			}
			elsif ($inicuted >= $subl3[1] && $endcuted >= $subl3[2]){
				$printcds .= "$subl3[0]\tAnnotGFF\tCDS\t$inicuted\t$subl3[2]\t\.\t\.\t\.\tParent=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
				$prelength += ($subl3[2] - $inicuted);
			}
			else {
				die "No considero alguna opción en get_annot_genes_gff.pl\n";
			}

		}

		print Results "$subl2[0]\tAnnotGFF\tmRNA\t$inigene\t$endgene\t\.\t\.\t\.\tID=$subline[0];$subline[4];$subline[1];Pos:$subline[2]-$subline[3]\n";
		print Results "$printcds";

	} else {die "No encuentra en el gff $subline[0]\n";}
	
}
close File;
close Results;





