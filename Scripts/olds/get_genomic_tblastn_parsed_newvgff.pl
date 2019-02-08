#!/usr/bin/perl
use strict;
use warnings;

# usage: perl get_tblastn_parsed outfmt6_file nombre_out 

my ($line, $name, $nameout);
my (%blast, %fasta, %hits);

# Abro el gff y lo guardo en memoria
my %gff;
open (GFFfile , "<", $ARGV[2]); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);
	next if ($subline[2] !~ /CDS/);
	my $genename = "";
	if ($subline[8] =~ /Parent=([^;]+)/){
	#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
		$genename = $1;
	}
	else {die "falla al pillar el nombre geneid del gff $line\n";}
	if ($subline[3] < $subline[4]){
		push ( @{$gff{$subline[0]}{$genename}}, join ("\t" , $subline[3], $subline[4])) ;
	}
	elsif ($subline[3] > $subline[4]){
		push ( @{$gff{$subline[0]}{$genename}}, join ("\t" , $subline[4], $subline[3])) ;
	}
}
close GFFfile;

#Abriendo resultados Blast
my %countsgenescatched;

open (Blastfile , "<", $ARGV[0]); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > 1e-5); # Filtro por evalue
	my $frame = "";

	## NUEVO FILTRO, SI LA POSICIÓN PILLADA ES UN CDS ANOTADO, NO ME LA QUEDO Y CUENTO EL HIT EN EL GEN QUE HA DADO
	
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
	if (exists $gff{$key}){
		GFFBUCLE: foreach my $gffgene (keys %{$gff{$key}}) {
			foreach my $positions (@{$gff{$key}{$gffgene}}) {
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
#	else {print "Warning, el scaffold no está en el GFF o se llama distinto: scaffold $key\n";}
	next if ($cdspillado > 0);

	# Procedo ya a guardar los hits

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

#Parseando los hits para quitar posiciones solapantes

my %hits_parsed;
foreach my $key (sort keys %hits) {
	foreach my $position (@{$hits{$key}}){
		my @subpos = split (/\t/, $position);
		my $ini = $subpos[1];
		my $frame = $subpos[0];
		my $fin = $subpos[2];
		my $n = 0;
		if (exists $hits_parsed{$key}) {
			my $i = 0;
			foreach my $hit (@{$hits_parsed{$key}}){
				my @subhit = split (/\t/, $hit);
				if ($ini <= $subhit[1] && $fin >= $subhit[2]){
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$ini, $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $fin <= $subhit[2] && $fin >= ($subhit[1]-10)){
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($ini >= $subhit[1] && $fin >= $subhit[2] && $ini <= ($subhit[2]+10)){
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$subhit[1], $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $ini >= ($subhit[1] - 10)){
					@{$hits_parsed{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($fin >= $subhit[2] && $fin <= ($subhit[2] + 10)){
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
					print "Script: get_genomic_tblastn_parsed_newv.pl hay alguna situación no contemplada";
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

#Segunda ronda de parseo
my %hits_parsed_2nd;
foreach my $key (sort keys %hits_parsed) {
	foreach my $position (@{$hits_parsed{$key}}){
		my @subpos = split (/\t/, $position);
		my $ini = $subpos[1];
		my $frame = $subpos[0];
		my $fin = $subpos[2];
		my $n = 0;
		if (exists $hits_parsed_2nd{$key}) {
			my $i = 0;
			foreach my $hit (@{$hits_parsed_2nd{$key}}){
				my @subhit = split (/\t/, $hit);
				if ($ini <= $subhit[1] && $fin >= $subhit[2]){
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$ini, $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $fin <= $subhit[2] && $fin >= ($subhit[1]-10)){
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($ini >= $subhit[1] && $fin >= $subhit[2] && $ini <= ($subhit[2]+10)){
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$subhit[1], $fin);
					$n = 0;
					last;	
				}
				elsif ($ini <= $subhit[1] && $ini >= ($subhit[1] - 10)){
					@{$hits_parsed_2nd{$key}}[$i] = join("\t",$subhit[0],$ini, $subhit[2]);
					$n = 0;
					last;	
				}
				elsif ($fin >= $subhit[2] && $fin <= ($subhit[2] + 10)){
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
					print "Script: get_genomic_tblastn_parsed_newv.pl hay alguna situación no contemplada";
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
		print Results "$key\_$i Frame$subhit[0] complete $subhit[1] $subhit[2] tblastn\n";
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



