#!/usr/bin/perl
use strict;
use warnings;

## Script to generate statistics about intron length in a specific GFF

# usage: perl get_intron_size_fromgff.pl gff3_file


die "Usage: insert the GFF to estimate intron length statistics. You can modify inside the script the intron length to calculate the number of introns longer than the input number. Default is 15,000 bp\n" unless @ARGV == 1;


#######################################################################################################
# Edit this variable to test number of introns longer than the specified number. Default is 16,000 bp #
#######################################################################################################

my $maxintronlength = "15000";

###


my $totalintrons = 0;
my $totalintronsize = 0;
my $longerintrons = 0;
my $maxintron = "0";
my $minintron = "99999999999999999999999999999999";
my @intronslength;

my $input=$ARGV[0];
my ($eachline,@exons);
my $first=0;
open (IN, "<$input") or die ("no such file!");
while(defined($eachline=<IN>)){
	if($eachline=~/\tmRNA\t/ || $eachline=~/\ttranscript\t/ ){
		$first++;
		if($first != 1){
			if(scalar(@exons)>2){
				my @ordered_exons=sort {$a<=>$b} @exons;
				for (my $i=1;$i<=scalar(@ordered_exons)-3;$i=$i+2){
					my $each_intron_size=$ordered_exons[$i+1]-$ordered_exons[$i]-1;
					#print "$each_intron_size\t";
					$totalintrons++;
					$totalintronsize += $each_intron_size;
					push(@intronslength, $each_intron_size);
					if ($each_intron_size > $maxintron){
						$maxintron = $each_intron_size;
					}
					if ($each_intron_size < $minintron){
						$minintron = $each_intron_size;
					}
					if ($each_intron_size > $maxintronlength){
						$longerintrons++;
					}
				}
  			}#else{print "0";}
  		#print "\n";
		}
		@exons=();
	}elsif($eachline=~/\tCDS\t/){ 
		my @eachline=split(/\t/,$eachline);
		if ($eachline[3] <= $eachline[4]){
			push (@exons, $eachline[3],$eachline[4]);
		} else {
			push (@exons, $eachline[4],$eachline[3]);
		}
	}
}

if(scalar(@exons)>2){
	my @ordered_exons=sort {$a<=>$b} @exons;
	for (my $i=1;$i<=scalar(@ordered_exons)-3;$i=$i+2){
		my $each_intron_size=$ordered_exons[$i+1]-$ordered_exons[$i]-1;
		#print "$each_intron_size\t";
		$totalintrons++;
		$totalintronsize += $each_intron_size;
		push(@intronslength, $each_intron_size);
		if ($each_intron_size > $maxintron){
			$maxintron = $each_intron_size;
		}
		if ($each_intron_size < $minintron){
			$minintron = $each_intron_size;
		}		
		if ($each_intron_size > $maxintronlength){
			$longerintrons++;
		}
	}
}#else{print "0";}
#print "\n";

my $avelength = sprintf("%.4f", $totalintronsize/$totalintrons);
my $longerperc =  sprintf("%.4f", (($longerintrons/$totalintrons) * 100));
my $medianlength = median(@intronslength);
my $perc95 = percentile95(@intronslength);
my $perc99 = percentile99(@intronslength);
print "Totalintrons: $totalintrons\nAverageIntronLength: $avelength\nMedianIntronLength: $medianlength\n";
print "Percentile95 IntronLength: $perc95\nPercentile99 IntronLength: $perc99\n";
print "MaxIntronLength: $maxintron\nMinIntronLength: $minintron\nIntronLongerThan$maxintronlength"."bp: $longerintrons\n";
print "PercentLongerIntrons:$longerperc%\n";

sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}

sub percentile95
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    my $perc = $len*0.95;
    return ($vals[int($perc)]);
}

sub percentile99
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    my $perc = $len*0.99;
    return ($vals[int($perc)]);
}
