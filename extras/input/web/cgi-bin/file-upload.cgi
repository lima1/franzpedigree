#!/usr/bin/perl
#
use strict;
use warnings;

use CGI;
use JSON;

my $Query = new CGI;
my $Json  = new JSON;

print $Query->header();

my $fh = $Query->upload('csv-path');

if (!$fh && $Query->cgi_error()) {
    JSONerror('Upload failed!');
}

my $sep_char         = defined $Query->param('sep_char') ? $Query->param('sep_char') : q{,};
my $has_header       = defined $Query->param('has_header') ? 1 : 0;
my $read_locus_names = defined $Query->param('read_locus_names') ? 1 : 0;
my $data_col         = defined $Query->param('data_col') ? ($Query->param('data_col')-1) : 1;
my $alleles_per_col  = 1;
my $missing_allele   = '0';
my $id_col           = defined $Query->param('id_col') ? ($Query->param('id_col')-1) : 0;
my $birth_col        = defined $Query->param('birth_col') ? ($Query->param('birth_col')-1) : -1;
my $death_col        = defined $Query->param('death_col') ? ($Query->param('death_col')-1) : -1;
my $sex_col          = defined $Query->param('sex_col') ? ($Query->param('sex_col')-1) : -1;
my $mother_col       = defined $Query->param('mother_col') ? ($Query->param('mother_col')-1) : -1;
my $location_col     = defined $Query->param('location_col') ? ($Query->param('location_col')-1) : -1;
my $sex_chars        = "F:M:?";
my $collapse         = 0;
my $dataset_title    = defined $Query->param('dataset_title') ? $Query->param('dataset_title') : 'FRANz';
my $line_number      = 0;
my @data;
my @locus_names;

LINE:
while ( my $line = <$fh> ) {
    $line_number++;
    if ($line_number == 1 && $has_header && $read_locus_names) {
        my @fields = split $sep_char, $line;
        COL:
        for (my $i = $data_col; $i <= $#fields; $i+=2) {
            my $locus_name = trim($fields[$i]);
            $locus_name =~ s/ a \z//xms;
            next COL if $locus_name =~ m{\A \s* \z}xms;
            push @locus_names, $locus_name;
        }
    }
    next LINE if ($line_number == 1 && $has_header);
    chomp $line;
    next LINE if $line =~ m{\A \s* \z}xms;

    my @fields = map { trim($_) } split $sep_char, $line;

    if (scalar @fields < 2) {
        JSONerror('File does not look like a CSV file (line ' . $line_number
        . ' has only one column)' );
    }    
    push @data, [ @fields ];
}


close $fh;

my ( $fc, $mc, $uc ) = split q{:}, $sex_chars;
my @pedigree_ids;
my @pedigree_arcs;

my %rows;
ROW:
for my $data (@data) {
    my $row =  create_data_row($data);
    next ROW if !$row;
    my $location = 0;
    if ($location_col > -1) {
        $location = $data->[$location_col];
    }
    if (!defined $rows{$location}) {
        $rows{$location} = [ $row ];
    }
    else {
        push @{$rows{$location}}, $row;
    }    
}

my $num_loci = calc_number_loci();

my %json_msg = ( success => JSON::true );

my $output = scalar(keys %rows) . " $num_loci / $dataset_title";

if (@locus_names) {
    if ($num_loci != scalar @locus_names) {
        JSONerror(
            'Could not read locus names in header. Uncheck the option or fix the header (' . join(q{,}, @locus_names) . ").\n" .
            'Make sure the locus names are in the same columns as the ' .
            'corresponding allele data.'
        );
    }
    if ($locus_names[0] =~ m{\A \d+ \z}xms) {
        JSONerror(
            'Locus names numeric, change at least the first ID it to non-numeric, for example L'. $locus_names[0] .' (' . join(q{,}, @locus_names) . ').' );

    }    
    $output .=  "\n" . join("\n", @locus_names);
}

for my $key (sort keys %rows) {
    my $lid = q{};
    if ($location_col > -1) {
        $lid = q{ } . $key;
    }
    $output .= "\n" . scalar(@{$rows{$key}}) . "$lid\n";
    $output .=  join "\n", @{$rows{$key}};
}

$json_msg{data} = $output;

my $pedigree = q{};
if ($mother_col> -1) {
    $pedigree = scalar(@pedigree_ids) . "\n" . join("\n", @pedigree_ids) .
    "\n" . join("\n", @pedigree_arcs);
    $json_msg{ped} = $pedigree;
}

print $Json->encode(\%json_msg);

sub calc_number_loci {
    my $max = -1;
    for my $data (@data) {
        if (scalar @{$data} > $max) {
            $max = scalar @{$data};
        }
    }
    if (($max - $data_col) % 2 != 0) {
        JSONerror('Uneven number of allele columns');
    }
    return ($max - $data_col) / (3-$alleles_per_col);
}

sub create_data_row {
    my ($data) = @_;
    my $line = sprintf "%10s", trim($data->[$id_col]);
    my @loci;
    
    if ($alleles_per_col == 1) {
        for my $i ( $data_col .. (scalar(@{$data})-1)) {
            if ($data->[$i] eq $missing_allele || $data->[$i] =~ m{\A \s* \z}xms) {
                $data->[$i] = '?';
            }    
        }    
        for (my $i=$data_col; $i < scalar(@{$data})-1; $i += 2) {
            push @loci, "$data->[$i]/$data->[$i+1]";
        }
        my $sex = '?';

        if ($sex_col > -1) {
            if (lc($data->[$sex_col]) eq lc($fc)) {
                $sex = 'F';
            }  
            elsif (lc($data->[$sex_col]) eq lc($mc)) {
                $sex = 'M';
            }
        }
        
        if ($mother_col > -1) {
            push @pedigree_ids, $line;
            my $mother_id = trim($data->[$mother_col]);

            if ($mother_id !~ m{\A \s* \z}xms) {
                push @pedigree_arcs, sprintf("%10s%10s", $mother_id,$line);
            }
        }
        my $birth = '?';
        my $death = '?';
        if ($birth_col > -1 && $data->[$birth_col] =~ /\d\d\d\d/) {
            $birth = $data->[$birth_col];
        }    
        if ($death_col > -1 && $data->[$death_col] =~ /\d\d\d\d/) {
            $death = $data->[$death_col];
        }    
        $line .= " 1 $birth $death $sex " . join(q{ }, @loci);
    }
    return $line;
}

sub JSONerror {
    my ( $msg ) = @_;
    my %json = ( success => JSON::false, msg => $msg );
    print $Json->encode(\%json);
    exit;
}


sub trim {
    my ($s) = @_;
    $s =~ s/^\s+|\s$//;
    return $s;
}

