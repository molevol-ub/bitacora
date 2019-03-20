#!/usr/bin/perl
use strict;
use warnings;

# Newv: Bug fixed: It kept the best blast hit positions; Now it combines all blast hits
# V2: If there are hits in non-overlapping positions in the protein (i.e. two fused genes by the structural annotation software), it splits the protein (flag _split)

# usage: perl get_blastp_parsed_newv2.pl blast_outfmt6_file(-outfmt "6 std qlen slen") Output_name E-value

my ($line, $name, $nameout);
my (%blast, %fasta);

my $evalue = $ARGV[2];


#Opening blast results

open (Blastfile , "<", $ARGV[0]); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > $evalue); # Evalue filtering
	push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[13], $subline[12], $subline[8], $subline[9], $subline[2]));
}
close Blastfile;

open (Results, ">", "$ARGV[1]blastp_parsed_list.txt");
foreach my $key (sort keys %blast) {
	my $blasthit2="";
	my $hitlvl = "0";
	my (@ini, @fin);
	$ini[0] = "99999999999999999999";
	$fin[0] = "1";
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[4]*2)/3;
		my $filtro2 = ($subline[3]*0.8);

		if ($subline[2] < 51) { # If the alignment is lower than 50aa, and smaller than 2/3 QUERY length protein used, it should contain a similarity higher than 80%. If not it is remove as false positive (small domains hitting non-related proteins)
			unless ($subline[2] >= $filtro1){
				next if ($subline[7] < 80);
			}
		}

		if ($subline[2] >= $filtro1 || $subline[2] >= $filtro2) { # BLAST filtering: alignment covering 2/3 of subject, or 80% of query
			$hitlvl = "2";
			my $n = 0;
			my $extrahit = 0;
			foreach my $i (@ini){
				my $f = $fin[$n];
				my $f2 = $f + 10;
				my $i2 = $i - 10;

				if ($i >= int($subline[5]) && $f <= int($subline[6])){
					$ini[$n] = int($subline[5]);
					$fin[$n] = int($subline[6]);
				}
				elsif ($i <= int($subline[5]) && $f < int($subline[6]) && $f2 >= int($subline[5])){
					$fin[$n] = int($subline[6]);
				}
				elsif ($i > int($subline[5]) && $f >= int($subline[6]) && $i2 <= int($subline[6])){
					$ini[$n] = int($subline[5]);
				}
				elsif ($f2 < int($subline[5]) || $i2 > int($subline[6])) {
					$extrahit++;
				}

				$n++;
			}

			if ($extrahit >= $n) {
					$ini[$n] = int($subline[5]);
					$fin[$n] = int($subline[6]);				
			}
		}
	}

	if ($hitlvl == 2){
		my $hits = scalar(@ini);

		if ($hits == 1){
			print Results "$key annot $ini[0] $fin[0] blastp\n";
		}
		else {
			my $n = 0;
			foreach my $i (@ini){
				my $f = $fin[$n];

				my $nn= $n+1;
				print Results "$key\_split$nn annot $ini[$n] $fin[$n] blastp\n";
				$n++;

			}
		}
	}

}
close Results;




