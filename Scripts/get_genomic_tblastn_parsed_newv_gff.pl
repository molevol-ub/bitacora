#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use lib dirname (__FILE__);
use Readgff qw (readgff);

# usage: perl get_tblastn_parsed_newv_gff.pl outfmt6_file outname gff_file evalue

my ($line, $name, $nameout);
my (%blast, %fasta, %hits);

my $evalue = $ARGV[3];
my $gff = $ARGV[2];

# Read GFF
my ($gffgeneref, $gffcdsref, $gffscafcdsref) = &readgff($gff);
my %gffcds = %$gffcdsref; 
my %gffgene = %$gffgeneref;
my %gffscafcds = %$gffscafcdsref;

# Opening blast output
my %countsgenescatched;

open (Blastfile , "<", $ARGV[0]); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > $evalue); # E-value filter
	my $frame = "";

	## tblastn positions in GFF annotation are deleted. Gene hitted is counted
	
	my $key = $subline[1];
	my $ini = "";
	my $fin = "";
	if ($subline[8] < $subline[9]){
		$ini = $subline[8];
		$fin = $subline[9];
	}
	elsif ($subline[8] > $subline[9]){
		$ini = $subline[9];
		$fin = $subline[8];

	}
	my $cdspillado = "0";
	if (exists $gffscafcds{$key}){
		GFFBUCLE: foreach my $gffgene (keys %{$gffscafcds{$key}}) {
			foreach my $positions (@{$gffscafcds{$key}{$gffgene}}) {
				my @subposition = split (/\t/, $positions);
				if ($ini <= $subposition[0] && $fin >= $subposition[1]){
					$countsgenescatched{$gffgene}++;
					$cdspillado++;
					last GFFBUCLE;
				}
				elsif ($ini <= $subposition[0] && $fin >= $subposition[0]){
					$countsgenescatched{$gffgene}++;
					$cdspillado++;
					last GFFBUCLE;
				}
				elsif ($ini <= $subposition[1] && $fin >= $subposition[1]){
					$countsgenescatched{$gffgene}++;
					$cdspillado++;
					last GFFBUCLE;
				}
				elsif ($ini >= $subposition[0] && $fin <= $subposition[1]){
					$countsgenescatched{$gffgene}++;
					$cdspillado++;
					last GFFBUCLE;
				}
			}
		}
	}
	#else {print "Warning, scaffold $key does not have any annotation\n";}
	next if ($cdspillado > 0);

	# Save hits not annotated in the GFF

	if ($subline[12] =~ /-1/){
		$frame = 4;
		my $left = $subline[14] - $subline[8];
		my $right = $subline[14] - $subline[9];
		push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[14], $subline[13], $frame, $left, $right, $subline[2]));
	}
	elsif ($subline[12] =~ /-2/){
		$frame = 5;
		my $left = $subline[14] - $subline[8];
		my $right = $subline[14] - $subline[9];
		push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[14], $subline[13], $frame, $left, $right, $subline[2]));
	}
	elsif ($subline[12] =~ /-3/){
		$frame = 6;
		my $left = $subline[14] - $subline[8];
		my $right = $subline[14] - $subline[9];
		push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[14], $subline[13], $frame, $left, $right, $subline[2]));
	}
	else {
		push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[14], $subline[13], $subline[12], $subline[8], $subline[9], $subline[2]));
	}	
}
close Blastfile;


foreach my $key (sort keys %blast) {
	my $blasthit2="";
	my $frame = "";
	my $ini = "";
	my $fin = "";
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		#my $filtro1 = ($subline[4]*2)/3;
		#my $filtro2 = ($subline[3]*0.8)/3;# Divido entre 3 porque son nucleótidos y estoy comparando con la longitud de alineamiento en proteina
		
		if ($subline[2] < 17) { # Si el alineamiento es menor a 50 nt (min length exon), que al menos tenga el 80% o más de identidad para quitar falsos positivos
			next if ($subline[8] < 80);
		}

		# ELIMINO EL FILTRO DE LONGITUD PORQUE PERDERÉ EXONES PEQUEÑOS
		#if ($subline[2] >= $filtro1 ) { # Filtrar el blast, o que cubra 2/3 del subject
		$frame = $subline[5];
		$ini = int($subline[6]/3);
		$fin = int($subline[7]/3);
		push (@{$hits{$key}}, join("\t",$frame,$ini, $fin));
		next; 
		#}
	}

}

# Parsing hits to extend overlapping positions

my %hits_parsed;
foreach my $key (sort keys %hits) {
	foreach my $position (@{$hits{$key}}){
		my @subpos = split (/\t/, $position);
		my $ini = $subpos[1];
		my $frame = $subpos[0];

		my $cadena = "";
		if ($frame == 1 || $frame == 2 || $frame == 3){
			$cadena = "F";
		} elsif ($frame == 4 || $frame == 5 || $frame == 6){
			$cadena = "R";
		} else {die "ERROR in get_genomic_tblastn_parsed_newv_gff.pl: Can't find frame in $position\n";}

		my $fin = $subpos[2];
		my $n = 0;
		if (exists $hits_parsed{$key}) {
			my $i = 0;
			foreach my $hit (@{$hits_parsed{$key}}){
				my @subhit = split (/\t/, $hit);

				if ($cadena =~ /F/){ # Filter according to chain
					if ($subhit[0] == 4 || $subhit[0] == 5 || $subhit[0] == 6){
						$n++;	
						$i++;
						next;
					}
				} elsif ($cadena =~ /R/){
					if ($subhit[0] == 1 || $subhit[0] == 2 || $subhit[0] == 3){
						$n++;	
						$i++;
						next;
					}
				} else {die "ERROR in get_genomic_tblastn_parsed_newv_gff.pl: Can't find frame in $position\n";}

				if ($ini <= $subhit[1] && $fin >= $subhit[2]){
					@{$hits_parsed{$key}}[$i] = join("\t",$frame,$ini, $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $fin <= $subhit[2] && $fin >= ($subhit[1]-10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($ini >= $subhit[1] && $fin >= $subhit[2] && $ini <= ($subhit[2]+10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$subhit[1], $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $ini >= ($subhit[1] - 10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($fin >= $subhit[2] && $fin <= ($subhit[2] + 10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$subhit[1], $fin);
					$n = 0;
					last;		
				}		
				elsif ($ini > $subhit[2] || $fin < $subhit[1]) {
					$n++;
				}
				elsif ($ini > $subhit[1] && $fin < $subhit[2]){
					$n = 0;
					last;
				}
				else {
					# Debugging
					print "Debugging error: Script: get_genomic_tblastn_parsed_newv_gff.pl There is some extra situation\n";
				}	
				$i++;		 
			}
			if ($n > 0) {
				push (@{$hits_parsed{$key}}, join("\t",$frame,$ini, $fin));
			}
		}
		else{
			push (@{$hits_parsed{$key}}, join("\t",$frame,$ini, $fin));
		}
		
	}
}

# Second round of parsing

my %hits_parsed_2nd;
foreach my $key (sort keys %hits_parsed) {
	foreach my $position (@{$hits_parsed{$key}}){
		my @subpos = split (/\t/, $position);
		my $ini = $subpos[1];
		my $frame = $subpos[0];

		my $cadena = "";
		if ($frame == 1 || $frame == 2 || $frame == 3){
			$cadena = "F";
		} elsif ($frame == 4 || $frame == 5 || $frame == 6){
			$cadena = "R";
		} else {die "ERROR in get_genomic_tblastn_parsed_newv_gff.pl: Can't find frame in $position\n";}

		my $fin = $subpos[2];
		my $n = 0;
		if (exists $hits_parsed_2nd{$key}) {
			my $i = 0;
			foreach my $hit (@{$hits_parsed_2nd{$key}}){
				my @subhit = split (/\t/, $hit);

				if ($cadena =~ /F/){ # Filter according to chain
					if ($subhit[0] == 4 || $subhit[0] == 5 || $subhit[0] == 6){
						$n++;	
						$i++;
						next;
					}
				} elsif ($cadena =~ /R/){
					if ($subhit[0] == 1 || $subhit[0] == 2 || $subhit[0] == 3){
						$n++;	
						$i++;
						next;
					}
				} else {die "ERROR in get_genomic_tblastn_parsed_newv_gff.pl: Can't find frame in $position\n";}

				if ($ini <= $subhit[1] && $fin >= $subhit[2]){
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$frame,$ini, $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $fin <= $subhit[2] && $fin >= ($subhit[1]-10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($ini >= $subhit[1] && $fin >= $subhit[2] && $ini <= ($subhit[2]+10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$subhit[1], $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $ini >= ($subhit[1] - 10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($fin >= $subhit[2] && $fin <= ($subhit[2] + 10)){
					unless ($frame == $subhit[0]){ #Nuevo filtro para no alargar si los frames son distintos
						$n = 0;
						last;
					}
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$subhit[1], $fin);
					$n = 0;
					last;		
				}		
				elsif ($ini > $subhit[2] || $fin < $subhit[1]) {
					$n++;
				}
				elsif ($ini > $subhit[1] && $fin < $subhit[2]){
					$n = 0;
					last;
				}
				else {
					# Debugging
					print "Debugging error: Script: get_genomic_tblastn_parsed_newv_gff.pl There is some extra situation\n";
				}	
				$i++;		 
			}
			if ($n > 0) {
				push (@{$hits_parsed_2nd{$key}}, join("\t",$frame,$ini, $fin));
			}
		}
		else{
			push (@{$hits_parsed_2nd{$key}}, join("\t",$frame,$ini, $fin));
		}
		
	}
}


open (Results, ">", "$ARGV[1]tblastn_parsed_list.txt");
foreach my $key (sort keys %hits_parsed_2nd) {
	my $i = 0;
	foreach my $hit (@{$hits_parsed_2nd{$key}}){
		$i++;
		my @subhit = split (/\t/, $hit);
		print Results "$key\_exon$i Frame$subhit[0] complete $subhit[1] $subhit[2] tblastn\n";
	}
}
close Results;

open (Results, ">", "$ARGV[1]tblastn_genescatched_list.txt");
my $totalgenescatched = scalar (keys %countsgenescatched);
print Results "Totalgenescatched $totalgenescatched\n";
foreach my $key (sort keys %countsgenescatched) {
	print Results "$key $countsgenescatched{$key}\n";
}
close Results;



