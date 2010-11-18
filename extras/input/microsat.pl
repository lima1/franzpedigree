#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Fatal qw(open close);

open my $IN, '<', $ARGV[0];

my %freqs;
while (my $line = <$IN>) {
    chomp $line;
    next if $line =~ m{\A \s+ \d}xms;
    my ( $pop_id, $locus, $allele, $abs ) = split /\s+/, $line;
    $freqs{$locus}{$allele} += $abs;
}

warn Dumper \%freqs;
print scalar(keys %freqs) . "\n";
for my $locus (sort keys %freqs) {
    my @sorted = sort { $a <=> $b } keys %{$freqs{$locus}};
    print scalar(keys %{$freqs{$locus}}) . " $sorted[0] $sorted[-1] $locus\n";
    my $sum = 0;
    for my $allele ( @sorted ) {
        $sum += $freqs{$locus}{$allele};
    }    
    for my $allele ( @sorted ){
        print "$allele " . ($freqs{$locus}{$allele}/$sum) . "\n";
    }    
}    
