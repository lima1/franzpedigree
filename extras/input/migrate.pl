#!/usr/bin/perl

use strict;
use warnings;

use Fatal qw(open close);
use Data::Dumper;

my $ignore_missing = $ARGV[1];

my $i = -1;
my @lines;
my %genotypes;

open my $IN, '<', $ARGV[0];
LINE:
while (my $line = <$IN>) {
    chomp $line;
    $i++;
    if ($i < 2) {
        push @lines, $line;
        next LINE;
    }
    my $id = substr $line, 0, 10, q{};
    next LINE if ($ignore_missing && $line =~ m{\?}xms);

    if (!exists $genotypes{$line}) {
        $genotypes{$line} = [];
    }
    push @{$genotypes{$line}}, $id;
}    
close $IN;

$lines[1] = scalar keys %genotypes;

push @lines, map { $genotypes{$_}->[0] . q{ } . scalar(@{$genotypes{$_}}) . q{ } .
$_} keys %genotypes;

print join "\n", @lines;    
