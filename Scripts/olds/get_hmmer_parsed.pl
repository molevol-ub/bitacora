#!/usr/bin/perl
use strict;
use warnings;

#usage:

#  perl get_hmmer_frames_parsed.pl nombreficherosdeframes_tblout_hastaantes_del_frame outname

# 17 y 18 posiciones inicio y fin


my $hits = "";
my ($line, $name, $id);
my (%frame, %ini, %fin);

my $file= "$ARGV[0]";
open (Blastfile , "<", $file);
while (<Blastfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	$line =~ s/\s+/\t/g;$line =~ s/\t+/\t/g;
	my @subline = split (/\t/, $line);
	if ($subline[6] <= "1e-5"){ # E-value para filtrar
		next if ($subline[12] > "1e-5");
		if ($hits !~ /$subline[0]\_/) {
			$ini{$subline[0]} = $subline[17];
			$fin{$subline[0]} = $subline[18];
			$hits .= $subline[0];
			$hits .= "_";
		}
		else {
			if ($subline[17] < $ini{$subline[0]} ) {
				$ini{$subline[0]} = $subline[17];
			}
			if ($subline[18] > $fin{$subline[0]} ) {
				$fin{$subline[0]} = $subline[18];
			}
		}
	}
}
close Blastfile;



open (Results , ">", "$ARGV[1]hmmer_parsed_list.txt");
foreach my $key (sort keys %ini) {
	print Results "$key - $ini{$key} $fin{$key} hmmer\n";	
}







