#!/usr/bin/perl
use strict;
use warnings;


# Script to confirm that translated proteins from GFF match with those in proteome file using also gene length coincidence


my ($line, $name, $nameout);
my (%blast, %fasta);
my @fastaid;



open (Blastfile , "<", $ARGV[1]); 
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > 1e-5); # E-value filter
	push (@{$blast{$subline[0]}}, join("\t",$subline[1], $subline[2], $subline[4], $subline[10], $subline[12], $subline[13], $subline[3]));
}
close Blastfile;

open (Fasta , "<", $ARGV[0]); 
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/){
		$name = $1;
		push (@fastaid, $name);
	} else {
		next;
	}
}
close Fasta;


open (Results, ">", "$ARGV[2]\_prots_VsGFF_badannot_list.txt");
open (Results2, ">", "$ARGV[2]\_prots_VsGFF_goodannot_list.txt");
foreach my $key (@fastaid) {
	my $good = 0;

	my $nodomid = "";
	if ($key =~ /(\S+)\_(\d+)dom/){
		$nodomid = $1;
	}

	#$key =~ s/\./\\\./g;
	#$key =~ s/\|/\\\|/g;
	my $skey = $key;
	$skey =~ s/\|/\\\|/g;	

	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		if ($subline[0] =~ /^$skey$/){
			if ($subline[1] > 89){ ## ID filter

				my $klength = $subline[4];
				my $slength = $subline[5];
				my $alnlength = $subline[6];
				my $expected = $slength * 0.9;

				if ($alnlength > $expected){ ### Length filter
					$good = 1;
				}
			}
		} elsif ($subline[0] =~ /^$nodomid$/){
			if ($subline[1] > 89){ ## ID filter

				my $klength = $subline[4];
				my $slength = $subline[5];
				my $alnlength = $subline[6];
				my $expected = $slength * 0.9;

				if ($alnlength > $expected){ ### Length filter
					$good = 1;
				}
			}
		}

	}
	if ($good == 1){
		print Results2 "$key\n";
	} elsif ($good == 0){
		print Results "$key\n";
	} 
}
close Results;
close Results2;




