#!/usr/bin/perl
use strict;
use warnings;

# Script to exlude similar sequences in a fasta file, retaining the longest sequence as representative (useful to cluster isoforms or proteins resulting from putative assembly artifacts)

die "\nUsage: insert the protein fasta to be filtered and the minimum length (aa) required to retain sequences. BLAST needs to be installed in the system.\nperl exclude_similar_sequences_infasta.pl input.fasta 100\n\nAdditional parameters:\nThe % identity used to filter similar sequences can be modified in line 11 of the script (default is 98)\nOutput name to be added in each sequence can be included in line 12\nThe number of threads to be used in blast search can be input in line 13 (default is 1)\n\n" unless @ARGV == 2;


my $filtlen = "$ARGV[1]"; # Length to filter sequences
my $ident = "98"; # Percent of identity to filter sequences
my $nametag = ""; # Name to add in out sequences
my $threads = "1"; # Threads to use in blastp search


my $in = $ARGV[0];

my %fasta; my $name;

my $outname="";
if ($in =~ /(\S+)\.fa/){
	$outname = $1;
} else {
	$outname = $in;
}

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

system ("makeblastdb -in $in -dbtype prot >/dev/null");
system ("blastp -query $in -db $in -out $in\.outfmt6 -evalue 1e-3 -num_threads $threads -outfmt \"6 std qlen slen\" ");

open (Blastfile , "<", "$in\.outfmt6"); 
while (<Blastfile>) {
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\t/, $line);
	next if ($subline[10] > 1e-5); # Evalue filter
	push (@{$blast{$subline[1]}}, join("\t",$subline[0], $subline[10], $subline[3], $subline[13], $subline[12], $subline[8], $subline[9], $subline[2]));
}
close Blastfile;

my $exclude = "";
open (Results2, ">", "$outname\_idseqsclustered_table.txt");
foreach my $key (sort keys %blast) {
	foreach my $blastresult (@{$blast{$key}}) {
		my @subline = ();
		@subline = split (/\t/, $blastresult);
		my $hitblast = $key;
		my $filtro1 = ($subline[4]*0.8);
		my $filtro2 = ($subline[3]*0.8);

		next if ($subline[2] < $filtro1 && $subline[2] < $filtro2); # Alignment length filter

		### New script

		my $escaped = $subline[0];
		$escaped =~ s/\|/\\\|/g; # Escape \|

		my $escapedk = $key;
		$escapedk =~ s/\|/\\\|/g; # Escape \|			

		next if ($key =~ /^$escaped$/);

		if ($subline[7] >= $ident){ # Identity filter
			my $length1 = length($fasta{$key});
			my $length2 = length($fasta{$subline[0]});	

			print Results2 "$key\t$subline[0]\t$subline[7]\t$subline[2]\t$length1\t$length2\n";

			if ($length1 > $length2){
				$exclude .= "$subline[0] ";
			} elsif ($length1 < $length2) {
				$exclude .= "$key ";
			} else { # same length proteins
				next if ($exclude =~ /$escapedk\s/);
				next if ($exclude =~ /$escaped\s/);
				$exclude .= "$subline[0] ";
				#print "$key\n"; # debugging
			}

		}
	}

}
close Results2;

open (Results, ">", "$outname\_idseqsclustered.fasta");
foreach my $key (sort keys %fasta) {
	my $length3 = length($fasta{$key});

	my $escapedk = $key;
	$escapedk =~ s/\|/\\\|/g; # Escape \|		
	next if ($exclude =~ /$escapedk\s/); # Exclude similar seqs

	next if ($length3 < $filtlen); # Length filter
	
	print Results ">$nametag"."$key\n$fasta{$key}\n";
}
close Results;


print "\nDone, output files have been saved in $outname\_idseqsclustered.fasta and $outname\_idseqsclustered_table.txt\n\n";




