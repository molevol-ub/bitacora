#!/usr/bin/perl
use strict;
use warnings;

# Remove from fasta duplicated sequences

#usage: perl get_fasta_uniqueseqs.pl infasta outfasta

my ($name, $line, $nameout, $contig, $frame, $name2);
my (%nrfa, %inicio, %fin, %fastacutted, %fastanames);
my $skip = 0;
my $ids = "";
my $numseqsinput = 0;

# read and print fasta with one seq per line
open (Res, ">", "$ARGV[1]\.tmp");
my $file= "$ARGV[0]";
my %fasta;
open (Fasta , "<", $file);
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		$numseqsinput++;
		if (exists $fasta{$name}){
			$skip = 1;
		} else {
			$skip = 0;
		}
	} else {
		next if ($skip == 1);
		$fasta{$name} .= $line;
	}
}
close Fasta;

foreach my $key (sort keys %fasta){
	print Res ">$key\n$fasta{$key}\n";
}
close Res;

if ($numseqsinput > 20000){
	print "Warning, you input database contains more that 20,000 protein sequences. This may take a while..\n";
}

# delete repeated seqs

$skip = 0;

open (Results, ">", "$ARGV[1]");

open (Fasta , "<", "$ARGV[1]\.tmp");
while (<Fasta>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>(\S+)/) {
		$name = $1;
		if (exists $nrfa{$name}){
			$skip = 1;
		} else {
			$skip = 0;
		}
	}
	else {
		next if ($skip == 1);
		if ($numseqsinput > 20000){
			$nrfa{$name} = $line;
			next;
		}

			my $match = "0";
			my $hit = "";
			my $hitl = "";
			foreach my $kseq (keys %nrfa){
				my $individualseq = $nrfa{$kseq};
				if ($individualseq =~ /$line/){
					$match++;
					$hit = $kseq;
					$hitl = length($individualseq);
				}
				elsif ($line =~ /$individualseq/){
					$match++;
					$hit = $kseq;
					$hitl = length($individualseq);
				}
			}

			if ($match > 0){
				my $seql = length($line);
				if ($seql <= $hitl){
					next;
				} elsif ($seql > $hitl){
					$nrfa{$name} = $line;
					delete $nrfa{$hit};
				}
			} else {
				$nrfa{$name} = $line;
			}

	}
}
close Fasta;

foreach my $key (sort keys %nrfa){
	my $nkey = $key;
	$nkey =~ s/\_//g;
	#$nkey =~ s/\.//g;
	#$nkey =~ s/\-//g;
	$nkey =~ s/\,//g;
	$nkey =~ s/\;//g;
	$nkey =~ s/\://g;
	$nkey =~ s/\|//g;
	print Results ">$nkey\n$fasta{$key}\n";
}

close Results;

