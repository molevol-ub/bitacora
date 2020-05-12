#!/usr/bin/perl
use strict;
use warnings;

if(!$ARGV[4]) {  
	die "Split the input fasta to run parallel blast\nperl run_parallel_blast.pl query database num_threads evalue output_file\n";
}

my $query = $ARGV[0];
my $database = $ARGV[1];
my $nthreads = $ARGV[2];
my $eval = $ARGV[3];
my $outfile= $ARGV[4];

my $seqsperfile = ""; 

my $numseq = 0; 
my $n = 1;

my %fasta;
my $name;
my @numfiles;

# Read fasta
open (File, "<", $query);
while (<File>){
	chomp;
	my $line = $_;
	next if ($line !~ /\S+/);
	if ($line =~ /^>/){
		$name = $line;
		$numseq ++;
	}
	else {
		$fasta{$numseq}{$name} .= $line;
	}
}
close File;

$seqsperfile = int($numseq/$nthreads + 0.5);

# Split fasta
for (my $i = 1; $i <= $numseq; $i+=  $seqsperfile){

	open (Results, ">", "$query\_$n.tmpfa");
	my $endfile = ($i + $seqsperfile);
	for ( my $ii = $i; $ii < $endfile ; $ii ++){	
		foreach my $seq (keys %{$fasta{$ii}}) {
			print Results "$seq\n$fasta{$ii}{$seq}\n";
		}		
	}
	close Results;

	push (@numfiles, $n);
	$n++;

}

# Run parallel blast

foreach my $i (@numfiles) {
    my $pid = fork();
    if ($pid==0) { # child
        exec("tblastn -query $query\_$i\.tmpfa -db $database -out $query\_$i.blastout -outfmt \"6 std sallseqid score nident positive gaps ppos qframe sframe qseq sseq qlen slen salltitles\" -evalue $eval");
        die "ERROR in run_parallel_tblastn_gemoma.pl Exec $i failed: $!\n";
    } elsif (!defined $pid) {
        warn "ERROR in run_parallel_tblastn_gemoma.pl Fork $i failed: $!\n";
    }
}

1 while wait() >= 0;


# Join blast results and delete intermediate files

system ("cat $query\_*.blastout > $outfile");
system ("rm $query\_*.tmpfa $query\_*.blastout")



