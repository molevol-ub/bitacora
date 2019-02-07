#!/usr/bin/perl
use strict;
use warnings;

# usage: perl get_genomic_genes_from_exons.pl tblastn_parsed_list.txt genomic_exons name_out max_intron_size

my ($line, $name, $nameout);
my (%blast, %fasta, %hits);

my $maxintron = $ARGV[3];
my $maxintronprot = $maxintron/3;

# Reading exon sequences
open (Fasta , "<", $ARGV[1]);
while (<Fasta>) {
chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
	}
	else {
		$fasta{$name} .= $line;
	}
}
close Fasta;

# Reading exon positions

open (File , "<", $ARGV[0]); 
while (<File>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	my @subline = split (/\s/, $line);
	my $key = "";
	if ($subline[0] =~ /(\S+)\_exon\d+$/){
		$key = $1;
	}
	else {
		print "ERROR in get_genomic_genes_from_exons.pl: $subline[0] is not identified\n";
	}

	my $frame = "";
	if ($subline[1] =~ /Frame1/ || $subline[1] =~ /Frame2/ || $subline[1] =~ /Frame3/){
		$frame = "Forward";
	} elsif ($subline[1] =~ /Frame4/ || $subline[1] =~ /Frame5/ || $subline[1] =~ /Frame6/) {
		$frame = "Reverse";
	} else {
		die "ERROR in get_genomic_genes_from_exons.pl: $ARGV[0] in line $line does not contain frame\n";
	}

	$hits{$key}{$frame}{$subline[3]} = join("\t",$subline[0], $subline[3], $subline[4]);
}
close File;


# Joining exons into genes

open (Results, ">", "$ARGV[2]\_genomic_genes_proteins.fasta");
foreach my $key (sort keys %hits) {
	foreach my $frame (keys %{$hits{$key}}){
		my $prev_ini = 0;
		my $num_exons = 0;
		my $key_fasta = "";
		my $seq_exon = "";
		foreach my $ini (sort {$a <=> $b} keys %{$hits{$key}{$frame}}) {
			my @subhit = split (/\t/, $hits{$key}{$frame}{$ini});
			my $fin = $subhit[2];
			my $exon = $subhit[0];

			if (scalar(keys %{$hits{$key}{$frame}}) == 1 ){
				print Results ">$exon\n$fasta{$exon}\n";
				next;
			}
			else {
				if ($prev_ini == 0){
					$num_exons++;
					$key_fasta = $exon;
					$prev_ini = $ini;
					$seq_exon = $fasta{$exon};
				}
				elsif ($ini > $prev_ini && $ini <= ($prev_ini + $maxintronprot)){
					$num_exons++;
					$prev_ini = $ini;
					$seq_exon .= $fasta{$exon};
					if ($exon =~ /\S+(_exon\d+)/){
						$key_fasta .= $1;
					}
					else {
						die "ERROR in get_genomic_genes_from_exons.pl: It does not recognized exon id in $exon\n";
					}
				}			
				elsif ($ini > ($prev_ini + $maxintronprot)){
					print Results ">$key_fasta\n$seq_exon\n";
					$num_exons = 1;
					$prev_ini = $ini;
					$seq_exon = $fasta{$exon};
					$key_fasta = $exon;
				}
				else {
					die "ERROR in get_genomic_genes_from_exons.pl: A initial position is lower than previous one, bad order in $ARGV[2]: $hits{$key}{$frame}{$ini}\n";
				}
			}
		
		}
		if ($num_exons >= 1){
			print Results ">$key_fasta\n$seq_exon\n";
		}
	}
}

close Results;




