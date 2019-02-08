#!/usr/bin/perl
use strict;
use warnings;

#Antes printaba solo la mejor sin solapar posiciones, ahora se queda con las mejores posiciones

##Abre el blast y saca la secuencias con hit blast parseado AÑADIENDO EL FRAME
# usage: perl get_tblastn_parsed outfmt6_file nombre_out 

my ($line, $name, $nameout);
my (%blast, %fasta);


#Abriendo resultados Blast

open (Blastfile , "<", $ARGV[0]); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > 1e-5); # Filtro por evalue
	push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[13], $subline[12], $subline[8], $subline[9], $subline[2]));
}
close Blastfile;

open (Results, ">", "$ARGV[1]blastp_parsed_list.txt");
foreach my $key (sort keys %blast) {
	my $blasthit2="";
	my $hitlvl = "0";
	my $ini = "999999999999999999";
	my $fin = "1";
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[4]*2)/3;
		my $filtro2 = ($subline[3]*0.8);

		if ($subline[2] < 51) { # Si el alineamiento es menor a 50 aa, que al menos tenga el 80% o más de identidad para quitar falsos positivos
			next if ($subline[7] < 80);
		}

		if ($subline[2] >= $filtro1 ) { # Filtrar el blast, o que cubra 2/3 del subject o sea más de un 80% del query
			$hitlvl = "2";
			if ($ini > int($subline[5])){
				$ini = int($subline[5]);
			}
			if ($fin < int($subline[6])){
				$fin = int($subline[6]);
			}
			#last; # Que acabe el bucle porque ya ha encontrado el mejor hit
		}
		elsif ($subline[2] >= $filtro2){
			if ($ini > int($subline[5])){
				$ini = int($subline[5]);
			}
			if ($fin < int($subline[6])){
				$fin = int($subline[6]);
			}
			next if ($hitlvl == 2);
			$hitlvl = "1";
		}
	}
	if ($hitlvl == 1){
		print Results "$key fragment $ini $fin blastp\n";
	}
	elsif ($hitlvl == 2){
		print Results "$key complete $ini $fin blastp\n";
	}
}
close Results;




