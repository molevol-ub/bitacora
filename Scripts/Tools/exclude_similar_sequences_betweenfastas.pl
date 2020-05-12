#!/usr/bin/perl
use strict;
use warnings;

# Script to exlude similar sequences between two fasta files

die "\nUsage: insert the protein fasta files to be compared. BLAST needs to be installed in the system.\nperl exclude_similar_sequences_bewtweenfastas.pl input1.fasta input2.fasta\n\nAdditional parameters:\nThe % identity used to filter similar sequences can be modified in line 10 of the script (default is 98)\nOutput name to be added in each sequence can be included in line 11\nThe number of threads to be used in blast search can be input in line 12 (default is 1)\n\n" unless @ARGV == 2;


my $ident = "98"; # Percent of identity to filter sequences
my $nametag = ""; # Name to add in out sequences
my $threads = "1"; # Threads to use in blastp search

=h
my %length = ( # Longitud para filtrar cada familia por completas (50% de esa length para filtro de este script)
	'CCP' => '50',
	'CSP' => '65',
	'GR' => '235',
	'IR' => '180',
	'IRlong' => '200',
	'NPC2' => '100',
	'OBP' => '75',
	'SNMP' => '280',
	'OR' => '235',
);
=cut

my $in = $ARGV[0];
my $sin = $ARGV[1];

my %fasta; my $name;

open (Fasta , "<", "$in"); 
while (<Fasta>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ />(\S+)/){
		$name = $1;
	} else {
		$fasta{$name} .= $line;
	}

}
close Fasta;

my %blast;

system ("makeblastdb -in $sin -dbtype prot >/dev/null");
system ("blastp -query $in -db $sin -out $in\.outfmt6 -evalue 1e-3 -num_threads $threads -outfmt \"6 std qlen slen\" ");

open (Blastfile , "<", "$in\.outfmt6"); 
while (<Blastfile>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > 1e-5); # Evalue filter
	push (@{$blast{$subline[0]}}, join("\t",$subline[1], $subline[10], $subline[3], $subline[12], $subline[13], $subline[8], $subline[9], $subline[2]));
}
close Blastfile;

my $exclude = "";
open (Results2, ">", "$in\_blast_table.txt");
foreach my $key (sort keys %blast) {
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[3]*0.8);
		my $filtro2 = ($subline[4]*0.8);

		next if ($subline[2] < $filtro1 && $subline[2] < $filtro2); # Alignment length filter

		### New script

		#next if ($key =~ /^$subline[0]$/);

		if ($subline[7] >= $ident){ # Identity filter
			my $length1 = length($fasta{$key});
			my $length2 = $subline[4];	

			print Results2 "$key\t$subline[0]\t$subline[7]\t$subline[2]\t$length1\t$length2\n";

=h
			if ($length1 > $length2){
				$exclude .= "$subline[0] ";
			} elsif ($length1 < $length2) {
				$exclude .= "$key ";
			} else { # same length proteins
				next if ($exclude =~ /$key\s/);
				next if ($exclude =~ /$subline[0]\s/);
				$exclude .= "$subline[0] ";
				#print "$key\n"; # debugging
			}
=cut			
			$exclude .= "$key ";
		}
	}

}
close Results2;

open (Results, ">", "$in\_seqsinref.fasta");
open (Results2, ">", "$in\_seqsnotinref.fasta");
foreach my $key (sort keys %fasta) {
	my $length3 = length($fasta{$key});

	#next if ($length3 < $filtlen); # Length filter

	my $escapedk = $key;
	$escapedk =~ s/\|/\\\|/g; # Escape \|		
	
	if ($exclude =~ /$escapedk\s/) {
		print Results ">$nametag"."$key\n$fasta{$key}\n";
	} else {		
		print Results2 ">$nametag"."$key\n$fasta{$key}\n";
	}
}
close Results;
close Results2;

print "\nDone, output files have been saved in $in\_seqsinref.fasta (sequences found in both fasta files), $in\_seqsnotinref.fasta (sequences only found in $in) and $in\_blast_table.txt\n\n";




