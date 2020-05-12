#!/usr/bin/perl
use strict;
use warnings;

# Script to detect genes in overlapping positions (being putative isoforms or bad-annotations)
# It uses the positions in CDS fields for a gene, it is more precise than using the mRNA field, but it takes long to run with big gffs

# usage: perl get_overlapping_genes_fromgff.pl file.gff3
die "Script to detect genes in overlapping positions \n\nperl get_overlapping_genes_fromgff.pl file.gff3 \n\n" unless @ARGV == 1;


my ($line, $name, $nameout);
my %overlap; my %mrna;

open (Results, ">", "$ARGV[0]\_overlapping_genes.txt");

open (GFFfile , "<", $ARGV[0]); 
while (<GFFfile>) {
	chomp;
	$line = $_;
	next if ($line !~ /\S+/);
	next if ($line =~ /^#/);
	my @subline = split (/\t/, $line);

	if ($subline[2] =~ /CDS/){
	
		my $genename = "";
		if ($subline[8] =~ /Parent=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {die "ERROR in get_overlapping_genes_fromgff.pl: It fails detecting Parent ID in $line\n";}



#		if ($genename =~ /(\S+)\_\d+dom/){
#			$genename = $1;
#		} 



		#Overlap checking
		my $start = "";
		my $end = "";
		my $frame = "$subline[6]";
		$frame =~ s/\+/\\\+/g;
		$frame =~ s/\-/\\\-/g;
		if ($subline[4] > $subline[3]){
			$start = $subline[3];
			$end = $subline[4];
		} else {
			$start = $subline[4];
			$end = $subline[3];			
		}

		my $escapedgenename = $genename;
		$escapedgenename =~ s/\./\\\./g;
		$escapedgenename =~ s/\|/\\\|/g;


		foreach my $key (keys %{$mrna{$subline[0]}}){			
			my @subl2 = split (/\s/, $mrna{$subline[0]}{$key});
			next if ($subl2[2] !~ /$frame/); # First check chain
			my $start2 = $subl2[0];
			my $end2 = $subl2[1];
			my $over = 0;
			if ($start <= $start2 && $end >= $end2){
				$over++
			} elsif ($start2 <= $start && $end2 >= $end){
				$over++
			} elsif ($start <= $start2 && $end >= $start2){
				$over++
			} elsif ($start2 <= $start && $end2 >= $start){
				$over++
			} 

			if ($over > 0){
				#if (exists $overlap{$genename}){
				#	$overlap{$genename} .= "$key ";
				#} elsif (exists $overlap{$key}){
				#	$overlap{$key} .= "$genename ";
				#} else {

					if ($key =~ /^$escapedgenename$/){
						next;
					} elsif (exists $overlap{$subline[0]}{$key}){
						if ($overlap{$subline[0]}{$key} =~ /$escapedgenename\s/){
							next;
						}
					}

					$overlap{$subline[0]}{$genename} .= "$key ";
					$overlap{$subline[0]}{$key} .= "$genename ";
				#}
			}

		}

		if (exists ($mrna{$subline[0]}{$genename})){
			my @spl = split(/\s/, $mrna{$subline[0]}{$genename});
			my $nstart = $spl[0]; 
			my $nend = $spl[1];
			if ($start < $spl[0]){
				$nstart = $start;
			}
			if ($end > $spl[1]){
				$nend = $end;
			}
			$mrna{$subline[0]}{$genename} = "$nstart $nend $frame";
		} else {
			$mrna{$subline[0]}{$genename} = "$start $end $frame";

		}


	}
	elsif ($subline[2] =~ /mRNA/){

###Using CDS
next;

		next if ($line =~ /^END/);

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "ERROR in get_overlapping_genes_fromgff.pl: It fails detecting ID in $line\n";}

#		if ($genename =~ /(\S+)\_\d+dom/){
#			$genename = $1;
#		} 


		#Overlap checking
		my $start = "";
		my $end = "";
		my $frame = "$subline[6]";
		$frame =~ s/\+/\\\+/g;
		$frame =~ s/\-/\\\-/g;
		if ($subline[4] > $subline[3]){
			$start = $subline[3];
			$end = $subline[4];
		} else {
			$start = $subline[4];
			$end = $subline[3];			
		}

		foreach my $key (keys %{$mrna{$subline[0]}}){			
			my @subl2 = split (/\s/, $mrna{$subline[0]}{$key});
			next if ($subl2[2] !~ /$frame/); # First check chain
			my $start2 = $subl2[0];
			my $end2 = $subl2[1];
			my $over = 0;
			if ($start <= $start2 && $end >= $end2){
				$over++
			} elsif ($start2 <= $start && $end2 >= $end){
				$over++
			} elsif ($start <= $start2 && $end >= $start2){
				$over++
			} elsif ($start2 <= $start && $end2 >= $start){
				$over++
			} 

			if ($over > 0){
				#if (exists $overlap{$genename}){
				#	$overlap{$genename} .= "$key ";
				#} elsif (exists $overlap{$key}){
				#	$overlap{$key} .= "$genename ";
				#} else {
					$overlap{$subline[0]}{$genename} .= "$key ";
					$overlap{$subline[0]}{$key} .= "$genename ";
				#}
			}

		}

		$mrna{$subline[0]}{$genename} = "$start $end $frame";


	}

	elsif ($subline[2] =~ /gene/){

####### No needed in this script
next;
#######
		next if ($line =~ /^END/);

		my $genename = "";
		if ($subline[8] =~ /ID=([^;]+)/){
		#if ($subline[8] =~ /transcript_id \"(\S+)\"/){
			$genename = $1;
		}
		else {print "ERROR in get_overlapping_genes_fromgff.pl: It fails detecting ID in $line\n";}

		if ($genename =~ /(\S+)\_\d+dom/){
			$genename = $1;
		} 


	}
}
close GFFfile;

my %overlapgood;
foreach my $scaf (sort keys %overlap){
	my $exclude;
	foreach my $key (sort keys %{$overlap{$scaf}}){

		#print "$key $overlap{$scaf}{$key}\n"; #### Debug

		my @subl = split (/\s/, $overlap{$scaf}{$key});

		if (exists $overlapgood{$scaf}{$key}){
			foreach my $sgen (@subl){
				my $list = $overlapgood{$scaf}{$key};
				#$list =~ s/\./\\\./g;
				$list =~ s/\|/\\\|/g;
				if ($list =~ /$sgen /){
					next;
				} else {
					$overlapgood{$scaf}{$key} .= "$sgen ";
				}
			}
		} else {
			my $cont = 0;
			foreach my $sgen (@subl){
				if (exists $overlapgood{$scaf}{$sgen}){
					$cont++;
					my $list = "$sgen $overlapgood{$scaf}{$sgen}";
					#$list =~ s/\./\\\./g;
					$list =~ s/\|/\\\|/g;

					if ($list !~ /$key /){
						$overlapgood{$scaf}{$sgen} .= "$key ";
					}

					foreach my $ssgen (@subl){

						if ($list =~ /$ssgen /){
							next;
						} else {
							$overlapgood{$scaf}{$sgen} .= "$ssgen ";
						}
					}					

				}

			}

			if ($cont == 0){
				$overlapgood{$scaf}{$key} .= "$overlap{$scaf}{$key}";
			} elsif ($cont > 1){
				#print "Warning, two different fields have been saved for the same gene in $key $overlap{$scaf}{$key}\n";
			}
		}

	}
}

my %overlapgood2;
foreach my $scaf (sort keys %overlapgood){
	my $exclude;
	foreach my $key (sort keys %{$overlapgood{$scaf}}){

		#print "$key $overlapgood{$scaf}{$key}\n"; #### Debug

		my @subl = split (/\s/, $overlapgood{$scaf}{$key});

		if (exists $overlapgood2{$key}){
			foreach my $sgen (@subl){
				my $list = $overlapgood2{$key};
				#$list =~ s/\./\\\./g;
				$list =~ s/\|/\\\|/g;
				if ($list =~ /$sgen /){
					next;
				} else {
					$overlapgood2{$key} .= "$sgen ";
				}
			}
		} else {
			my $cont = 0;
			foreach my $sgen (@subl){
				if (exists $overlapgood2{$sgen}){
					$cont++;
					my $list = "$sgen $overlapgood2{$sgen}";
					#$list =~ s/\./\\\./g;
					$list =~ s/\|/\\\|/g;

					if ($list !~ /$key /){
						$overlapgood2{$sgen} .= "$key ";
					}

					foreach my $ssgen (@subl){

						if ($list =~ /$ssgen /){
							next;
						} else {
							$overlapgood2{$sgen} .= "$ssgen ";
						}
					}					

				}

			}

			if ($cont == 0){
				$overlapgood2{$key} .= "$overlapgood{$scaf}{$key}";
			} elsif ($cont > 1){
				#print "Warning, two different fields have been saved for the same gene in $key $overlapgood{$scaf}{$key}\n";
			}
		}

	}
}

foreach my $key (sort keys %overlapgood2){
	print Results "$key $overlapgood2{$key}\n";

}

close Results;

