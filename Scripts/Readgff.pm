package Readgff;
use strict;
use warnings;
use Exporter qw(import);
our @EXPORT_OK = qw(readgff);

# Perl module to read a GFF file


sub readgff {
	my ($inputgff) = @_;
	my %gffcds; my %gffgene; my %gffscafcds;
	open (GFFfile , "<", $inputgff) 
		or die "ERROR in Scripts/Readgff.pm:Could not open file '$inputgff' $!";
	while (<GFFfile>) {
		chomp;
		my $line = $_;
		next if ($line !~ /\S+/);
		next if ($line =~ /^#/);
		next if ($line !~ /\t/); # Avoid extra lines, i.e. GFF containing Fasta sequences at the end
		my @subline = split (/\t/, $line);

		if ($subline[2] =~ /CDS/){
			my $genename = "";
			if ($subline[8] =~ /Parent=([^;]+)/){
			#if ($subline[8] =~ /(transcript_id|Transcript) \"([^"]+)\"/){  # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
				$genename = $1;

				## Extra for ensembl GFF3
				if ($genename =~ /ranscript\:(\S+)/){
					$genename = $1;
				}
				#if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
				#	$genename = "$1"."-P$2"; ## Change P for R if protein is named as transcript
				#}

			} elsif ($subline[8] =~ /(transcript_id|Transcript)."([^"]+)"/){  # Reading GTF formats.
				$genename = $2;

			}
			else {die "ERROR in Scripts/Readgff.pm: It does not recognize Parent ID in the GFF3 in: $line\n";}


			if ($subline[6] =~ /\+/){ # Forward strand
				#OK
			} elsif ($subline[6] =~ /\-/){
				#OK
			} else {die "ERROR in Scripts/Readgff.pm: No forward/reverse in: $line\n";}


			if ($subline[7] =~ /\d+/){ # Forward strand
				#OK
			} else {die "ERROR in Scripts/Readgff.pm: No frame in CDS in: $line\n";}


			if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
				push ( @{$gffcds{$genename}{$subline[3]}}, $line);
				push ( @{$gffscafcds{$subline[0]}{$genename}}, join ("\t" , $subline[3], $subline[4])) ;		
			}
			 else {
			 	push ( @{$gffcds{$genename}{$subline[4]}}, join ("\t", $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]));
				push ( @{$gffscafcds{$subline[0]}{$genename}}, join ("\t" , $subline[4], $subline[3])) ;
			}
			
		}
		elsif ($subline[2] =~ /mRNA/ || $subline[2] =~ /transcript/){
			next if ($subline[2] =~ /transcriptional/);
			next if ($subline[2] =~ /transcription_end_site/);
			next if ($subline[2] =~ /transcription_start_site/);
			my $genename = "";
			if ($subline[8] =~ /ID=([^;]+)/){
			#if ($subline[8] =~ /(transcript_id|Transcript) \"([^"]+)\"/){  # Example to use in GTF formats. Uncomment this line (delete #) and comment previous line
				$genename = $1;

				## Extra for ensembl GFF3
				if ($genename =~ /ranscript\:(\S+)/){
					$genename = $1;
				}
				#if ($genename =~ /(\S+)\-R(\S+)/){ # Saving -RA as -PA to easily identify the annotated protein in next step
				#	$genename = "$1"."-P$2"; ## Change P for R if protein is named as transcript
				#}
			} elsif ($subline[8] =~ /(transcript_id|Transcript)."([^"]+)"/){  # Reading GTF formats.
				$genename = $2;

			}
			else {die "ERROR in Scripts/Readgff.pm: It does not recognize ID in the GFF3 in: $line\n";}

			if (exists $gffgene{$genename}){ # Control for duplicated genes in GFF3
				die "ERROR in Scripts/Readgff.pm: Gene $genename is duplicated in the GFF3, found duplicate in $line\nPlease, take a look into your GFF3 and delete duplicated genes\n";
			}


			if ($subline[6] =~ /\+/){ # Forward strand
				#OK
			} elsif ($subline[6] =~ /\-/){
				#OK
			} else {die "ERROR in Scripts/Readgff.pm: No forward/reverse in: $line\n";}


			if ($subline[4] > $subline[3]){ # Make sure that positions are ordered
				$gffgene{$genename} = $line;
			}
			 else {
			 	$gffgene{$genename} = join ("\t", $subline[0],$subline[1],$subline[2],$subline[4],$subline[3],$subline[5],$subline[6],$subline[7],$subline[8]);
			}

		}
	}
	close GFFfile;

	#print "GFF parsed correctly\n";

	return ( \%gffgene, \%gffcds, \%gffscafcds);

}

1;