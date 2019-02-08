#!/usr/bin/perl
use strict;
use warnings;

# Usage: 	perl Scripts/get_blast_hmmer_combined.pl outname\_allsearches_list.txt outname


open (File, "<", $ARGV[0]);
my %genes;
while (<File>){
	chomp;
	my $line = $_;
	my @subline = split(/\s/, $line);
	my $gene = "";
	if ($subline[0] =~ /(\S+)\_split\d/){
		$gene = $1;
	} else {
		$gene = $subline[0];
	}
	push (@{$genes{$gene}}, $line);
}
close File;
open (Results, ">", "$ARGV[1]\_combinedsearches_list.txt");

foreach my $gene (sort keys %genes){
	my $fragcomp = "-";	
	my (@ini, @fin);
	$ini[0] = "999999999999999999";
	$fin[0] = "1"; 
	my $method = "";
	foreach my $list (@{$genes{$gene}}){
		my @subline = split(/\s/, $list);

		if ($subline[4] =~ /blastp/){
			unless ($fragcomp =~ /complete/){
				$fragcomp = $subline[1];
			}
			unless ($method =~ /blastp/){
				$method .= "blastp";
			}
		} elsif ($subline[4] =~ /hmmer/){
			unless ($method =~ /hmmer/){
				$method .= "hmmer";
			}
		}

		my $n = 0;
		my $extrahit = 0;
		foreach my $i (@ini){
			my $f = $fin[$n];
			my $f2 = $f + 10;
			my $i2 = $i - 10;
				if ($i >= int($subline[2]) && $f <= int($subline[3])){
				$ini[$n] = int($subline[2]);
				$fin[$n] = int($subline[3]);
			}
			elsif ($i <= int($subline[2]) && $f < int($subline[3]) && $f2 >= int($subline[2])){
				$fin[$n] = int($subline[3]);
			}
			elsif ($i > int($subline[2]) && $f >= int($subline[3]) && $i2 <= int($subline[3])){
				$ini[$n] = int($subline[2]);
			}
			elsif ($f2 < int($subline[2]) || $i2 > int($subline[3])) {
				$extrahit++;
			}
			$n++;
		}
		if ($extrahit >= $n) {
				$ini[$n] = int($subline[2]);
				$fin[$n] = int($subline[3]);				
		}

	}

	my $hits = scalar(@ini);
	if ($hits == 1){
		print Results "$gene $fragcomp $ini[0] $fin[0] $method\n";
	}
	else {
		my $n = 0;
		foreach my $i (@ini){
			my $f = $fin[$n];
			my $nn = $n+1;
			print Results "$gene\_split$nn $fragcomp $ini[$n] $fin[$n] $method\n";
			$n++;
		}
	}

}
close Results;







