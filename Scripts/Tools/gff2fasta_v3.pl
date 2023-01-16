#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
my $dirname = dirname(__FILE__);

# Creates a CDS and protein fasta from a GFF and the genome file
# 

die "Usage: gff2fasta_v3.pl genome.fasta anotation.gff OutPrefix\n" unless @ARGV == 3;

my $genome = $ARGV[0];
my $gff = $ARGV[1];
my $gffcds = "$ARGV[2].cds.fasta";
my $gffprot = "$ARGV[2].pep.fasta";


system ("cat $gff | perl $dirname/gffGetmRNA.pl --genome=$genome --mrna=$gffcds");
system ("perl $dirname/translating_seqs_v3_onlyframe1.pl $gffcds $gffprot");
