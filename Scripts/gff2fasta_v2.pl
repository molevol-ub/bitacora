#!/usr/bin/perl
use strict;
use warnings;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Fasta;

die "Usage: gff2fasta_v2.pl genome.fasta anotation.gff OutPrefix\n" unless @ARGV == 3;

#add a help message here
#my $num_args=$#ARGV + 1;
#if ($num_args != 4) {
#	print "\nUsage: gff2perl Genome.fasta Annotation.gff OutputPrefix \n\n";
#	exit;
#}

$| = 1;    # Flush output
my $outfile_cds = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].cds.fasta" );
my $outfile_pep = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].pep.fasta" );
#my $outfile_cdna = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].cdna.fasta" );
#my $outfile_gene = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].gene.fasta" );
#my $outfile_upstream3000 = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].upstream3000.fasta" );
#my $outfile_exon = Bio::SeqIO->new( -format => 'fasta', -file => ">$ARGV[2].exon.fasta");

###### Output type description ######
# cds - translated sequence (starting with ATG and ending with a stop codon included)
# cdna - transcribed sequence (devoid of introns, but containing untranslated exons)
# protein - cds translated (includes a * as the stop codon)
# gene - the entire gene sequence (including UTRs and introns)
# upstream3000 - the 3000 upstream region of the gene (likely including the promoter)

### First, index the genome
my $file_fasta = $ARGV[0];
my $db = Bio::DB::Fasta->new($file_fasta);
print ("Genome fasta parsed\n");



### Second, parse the GFF3
my %CDS;
my %CDNA;
my %EXON;
my $mRNA_name;
my $frame;
open GFF, "<$ARGV[1]" or die $!;
while ( my $line = <GFF> ) {
    chomp $line;
    next if ($line !~ /\S+/);
    next if ($line =~ /^#/);

    my @array = split( "\t", $line );
    my $type = $array[2];

#CDS
    if ( ( $type eq 'mRNA' || $type eq 'transcript' ) and ( $. > 2 ) ) {
        # CDS: Collect CDSs and extract sequence of the previous mRNA
        my $mergedCDS_seq;
	   # WARNING we must sort by $cds_coord[1]


        foreach my $key (sort {$a <=> $b} keys %CDS) { # Ascending numeric sort of the starting coordinate
            my $coord = $CDS{$key};
            my @cds_coord = split( " ", $coord );
            my $cds_seq = $db->seq( $cds_coord[0], $cds_coord[1], $cds_coord[2] );
            $mergedCDS_seq .= $cds_seq;
        }
        

        my $output_cds = Bio::Seq->new(
            -seq        => $mergedCDS_seq,
            -id         => $mRNA_name,
            -display_id => $mRNA_name,
            -alphabet   => 'dna',
        );
        if ($frame eq '-') {
            $output_cds = $output_cds->revcom();
        }
	    #translate CDS to peptide for protein sequence 
        my $output_pep = $output_cds->translate();
    	#write to file
    	if (length($mergedCDS_seq) != 0) {
            $outfile_cds->write_seq($output_cds);
    	}
    	if (length($mergedCDS_seq) != 0) {
            $outfile_pep->write_seq($output_pep);
    	}

#exons
#should be able to add exon output here since exons will be useful in gene models for other organisms can be added in the EVM program

        # Now initialize the next mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $frame=$array[6];
        %CDS = (); %CDNA = (); # Empty the chunk arrays
	%EXON = (); %EXON = (); #Empty the EXON chunk arrays
    }
    elsif ( $type eq 'mRNA' || $type eq 'transcript' ) {    # First mRNA
        my @attrs = split( ";", $array[8] );
        $attrs[0] =~ s/ID=//;
        $mRNA_name = $attrs[0];
        $frame=$array[6];
    }
    elsif ( $type eq 'CDS' ) {
        my $cds_coord = $array[0] . " " . $array[3] . " " . $array[4];
        $CDS{$array[3]}=$cds_coord;
        $CDNA{$array[3]}=$cds_coord;
    }
    elsif (($type eq 'five_prime_UTR')||($type eq 'three_prime_UTR')) {
        my $utr_coord = $array[0] . " " . $array[3] . " " . $array[4];
        $CDNA{$array[3]}=$utr_coord;
    }
    elsif ($type eq 'exon' ) {
	my $exon_coord = $array[0] . " " . $array[3] . " " . $array[4];
	$EXON{$array[3]}=$exon_coord;
    }
}

close GFF;