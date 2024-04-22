#!/usr/bin/perl

use strict;

die "perl $0 <gff> <gff> > outfile \n" if @ARGV < 2;

my $Ref=shift;
my $File=shift;

my %gff;
read_gff($File,\%gff);

my %ref;
read_gff($Ref,\%ref);

print "#id1\tid2\tm_len1\tm_len2\tm_over\tper1\tper2\tc_len1\tc_len2\tc_over\tper1\tper2\n";
for my $name(keys %ref){
	for my $id(sort {$ref{$name}{$a}{mRNA}->[3] <=> $ref{$name}{$b}{mRNA}->[3]} keys %{$ref{$name}}){
		my $q=$ref{$name}{$id};
		my @aa=sort {$gff{$name}{$a}{mRNA}->[3] <=> $gff{$name}{$b}{mRNA}->[3]} keys %{$gff{$name}};
		for(@aa){
			my $p=$gff{$name}{$_};
			next if($p->{mRNA}->[4] < $q->{mRNA}->[3]);
			last if($p->{mRNA}->[3] > $q->{mRNA}->[4]);
			my $mRNA_over=overlap($q->{mRNA}->[3],$q->{mRNA}->[4],$p->{mRNA}->[3],$p->{mRNA}->[4]);
			my @q_cds=sort {$a->[3] <=> $b->[3]} @{$ref{$name}{$id}{CDS}};
			my @p_cds=sort {$a->[3] <=> $b->[3]} @{$gff{$name}{$_}{CDS}};
			my ($len_q_cds,$len_p_cds,$cds_over)=(0,0,0);
			for(@q_cds){
				$len_q_cds+=($_->[4]-$_->[3]+1);
			}
			for(@p_cds){
				$len_p_cds+=($_->[4]-$_->[3]+1);
			}
			for my $i(@q_cds){
				for my $j(@p_cds){
					next if($j->[4] < $i->[3]);
					last if($j->[3] > $i->[4]);
					$cds_over+=overlap($i->[3],$i->[4],$j->[3],$j->[4]);
				}
			}
			print "$id\t$_\t".($q->{mRNA}->[4]-$q->{mRNA}->[3]+1)."\t".($p->{mRNA}->[4]-$p->{mRNA}->[3]+1);
			print "\t$mRNA_over\t".($mRNA_over/($q->{mRNA}->[4]-$q->{mRNA}->[3]+1))."\t".($mRNA_over/($p->{mRNA}->[4]-$p->{mRNA}->[3]+1));
			print "\t$len_q_cds\t$len_p_cds\t$cds_over\t".($cds_over/$len_q_cds)."\t".($cds_over/$len_p_cds)."\n";
		}
	}
}


############
sub overlap{
	my ($aa,$ab,$ba,$bb)=@_;
	my $max=$ab>$bb?$ab:$bb;
	my $min=$aa<$ba?$aa:$ba;
	my $ov=0;
	$ov=($bb-$ba+$ab-$aa+2)-($max-$min+1);
	$ov=0 if ($ov < 0);
	return $ov;
}
sub read_gff{
	my ($file,$hash)=@_;
	open IN,$file or die $!;
	while(<IN>){
		next if /^#/;
		chomp;
		my @a=split /\t+/;
		next if($a[2] ne "mRNA" and $a[2] ne "CDS");
		my $name="$a[0]$a[6]";
		@a[3,4]=@a[4,3] if ($a[3] > $a[4]);
		my $id;
		if($a[2] eq 'mRNA' && $a[8] =~ /ID=([^;\s]+)/){
			$id=$1;
			@{$$hash{$name}{$id}{mRNA}}=@a;
		}elsif($a[2] eq 'CDS' && $a[8] =~ /Parent=([^;\s]+)/){
			$id=$1;
			push @{$$hash{$name}{$id}{CDS}},[@a];
		}
	}
	close IN;
}
