#!/usr/bin/perl
use strict;
use warnings;

# usage: perl get_hmmer_parsed_newv.pl hmmer_domtblout_outfile Output_name E-value

# newv: If there are hits in non-overlapping positions in the protein (i.e. two fused genes by the structural annotation software), it splits the protein (flag _split)

my ($line, $name, $nameout);
my (%blast, %fasta, %multipledom);

my $evalue = $ARGV[2];

#Opening HMMER results

open (Blastfile , "<", $ARGV[0]); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	$line =~ s/\s+/\t/g;$line =~ s/\t+/\t/g;
	my @subline = split (/\t/, $line);
	if ($subline[6] <= $evalue){ # E-value 
		next if ($subline[12] > $evalue); # E-value specific domain region
		push (@{$blast{$subline[0]}}, join("\t",$subline[0], $subline[17], $subline[18]));
		$multipledom{$subline[0]}++;
	}
}
close Blastfile;

open (Results, ">", "$ARGV[1]hmmer_parsed_list.txt");
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

		$hitlvl = "2";
		my $n = 0;
		my $extrahit = 0;
		foreach my $i (@ini){
			my $f = $fin[$n];
			my $f2 = $f + 10;
			my $i2 = $i - 10;

			if ($i >= int($subline[1]) && $f <= int($subline[2])){
				$ini[$n] = int($subline[1]);
				$fin[$n] = int($subline[2]);
			}
			elsif ($i <= int($subline[1]) && $f < int($subline[2]) && $f2 >= int($subline[1])){
				$fin[$n] = int($subline[2]);
			}
			elsif ($i > int($subline[1]) && $f >= int($subline[2]) && $i2 <= int($subline[2])){
				$ini[$n] = int($subline[1]);
			}
			elsif ($f2 < int($subline[1]) || $i2 > int($subline[2])) {
				$extrahit++;
			}

			$n++;
		}

		if ($extrahit >= $n) {
				$ini[$n] = int($subline[1]);
				$fin[$n] = int($subline[2]);				
		}

	}

	if ($hitlvl == 2){
		my $hits = scalar(@ini);

		if ($hits == 1){
			print Results "$key annot $ini[0] $fin[0] hmmer $multipledom{$key}\n";
		}
		else {
			my $n = 0;
			foreach my $i (@ini){
				my $f = $fin[$n];

				my $nn= $n+1;
				print Results "$key\_split$nn annot $ini[$n] $fin[$n] hmmer $multipledom{$key}\n";
				$n++;

			}
		}
	}


}
close Results;












