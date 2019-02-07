#!/usr/bin/perl
use strict;
use warnings;

my $name;
my %fasta;

open (Fasta, "<", $ARGV[0]);
while (<Fasta>){
    chomp;
    my $line = $_;
    if ($line =~ />(\S+)/){
        $name = $1;
    } else {
        $fasta{$name} .= $line;
    }
}
close Fasta;

open (Res, ">", "$ARGV[0].nostop.fasta");
foreach my $key (sort keys %fasta){
    my $seq = $fasta{$key};
    if ($seq =~ /(\S+)\*$/){
        print Res ">$key\n$1\n";
    } else {
        print Res ">$key\n$seq\n";
    }

}
close Res;