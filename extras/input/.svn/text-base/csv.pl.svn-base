#!/usr/bin/perl

use strict;
use warnings;

use Fatal qw(open close);
use Getopt::Long;
use Pod::Usage;
use Text::CSV;
use Data::Dumper;

our $VERSION = '0.1';

my ( $infile, $help, $man, $version );

my $sep_char        = q{,};
my $has_header      = 0;
my $data_col        = 1;
my $alleles_per_col = 1;
my $missing_allele  = '0';
my $sex_col         = -1;
my $id_col          = 0;
my $birth_col       = -1;
my $death_col       = -1;
my $sex_chars       = "F:M:?";
my $collapse        = 0;

my $options_ok = GetOptions(
    'in=s'         => \$infile,
    'sep_char=s'   => \$sep_char,
    'data_col=i'   => \$data_col,
    'alleles_per_col=i'   => \$alleles_per_col,
    'has_header'   => \$has_header,
    'collapse'     => \$collapse,
    'sex_col=i'    => \$sex_col,
    'id_col=i'     => \$id_col,
    'birth_col=i'  => \$birth_col,
    'death_col=i'  => \$death_col,
    'sex_chars=s'  => \$sex_chars,
    'missing_allele=s' => \$missing_allele,
    'help|?'       => \$help,
    'version|v'    => \$version,
    'man'          => \$man,
) or pod2usage(2);

if ($version) {
    print "$0 $VERSION\n";
    exit;
}

if ($man) {
    pod2usage( -exitstatus => 0, -verbose => 2 );
}
if ( $help || !defined $infile || ( $alleles_per_col != 1 && $alleles_per_col
    != 2) ) {
    pod2usage(1);
}

open my $IN, '<', $infile;
my $csv = Text::CSV->new(
    {   binary           => 1,
        sep_char         => $sep_char,
        allow_whitespace => 1,
    }
);

my @data;

my $line_number = 0;

LINE:
while ( my $line = <$IN> ) {
    $line_number++;
    next LINE if ($line_number == 1 && $has_header);
    chomp $line;
    next LINE if $line =~ m{\A \s* \z}xms;

    my $status = $csv->parse($line);
    push @data, [ $csv->fields() ];
}


close $IN;

my ( $fc, $mc, $uc ) = split q{:}, $sex_chars;

print "1 " . calc_number_loci() . " . csv.pl\n";
#print Dumper \@data;
my @rows;
ROW:
for my $data (@data) {
    my $row =  create_data_row($data);
    next ROW if !$row;
    push @rows, $row . "\n";
}

if ($collapse) {
    my $max = 0;
    my %ids;
    for my $row ( @rows ) {
        my $id = substr $row, 0, 10;
        $ids{$id}++;
        if ($ids{$id} > $max) {
            $max = $ids{$id};
        }    
    }
    my %seen;
    my @rows_collapsed;
    for my $row ( @rows ) {
        my $id = substr $row, 0, 10;
        my $rest = substr $row, 13;
        if (!defined $seen{$id}) {
            $seen{$id} = 1;
            push @rows_collapsed, "$id " . sprintf("%*s", length($max),
                $ids{$id}). " $rest";
        }
    }
    @rows = @rows_collapsed;
}

print scalar(@rows) . "\n";
print join q{}, sort @rows;

sub calc_number_loci {
    my $max = -1;
    for my $data (@data) {
        if (scalar @{$data} > $max) {
            $max = scalar @{$data};
        }
    }
    return ($max - $data_col) / (3-$alleles_per_col);
}

sub create_data_row {
    my ($data) = @_;
    my $line = sprintf "%10s", $data->[$id_col];
    my @loci;
    
    if ($alleles_per_col == 1) {
        for my $i ( $data_col .. (scalar(@{$data})-1)) {
            if ($data->[$i] eq $missing_allele || $data->[$i] =~ m{\A \s* \z}xms) {
                $data->[$i] = '?';
            }    
        }    
        for (my $i=$data_col; $i < scalar(@{$data})-1; $i += 2) {
            push @loci, "$data->[$i].$data->[$i+1]";
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
        my $birth = '?';
        my $death = '?';
        if ($birth_col > -1 && $data->[$birth_col] =~ /\d\d\d\d/) {
            $birth = $data->[$birth_col];
        }    
        if ($death_col > -1 && $data->[$death_col] =~ /\d\d\d\d/) {
            $death_col = $data->[$death_col];
        }    
        $line .= " 1 $birth $death $sex " . join(q{ }, @loci);
    }
    return $line;
}    
__END__

=head1 NAME

csv.pl

=head1 SYNOPSIS

  csv.pl [OPTIONS] --in in.csv

=head1 OPTIONS

=over

=item C<--sep_char c>

The separator character. Default comma ','.

=item C<--has_header>

Skips first line.

=item C<--id_col i>

The id of the genotype is at this column. First column is column 0. Default is
0.

=item C<--data_col i>

The marker data starts at this column. First column is column 0. Default is 1.

=item C<--birth_col i>

The column with the year of birth.

=item C<--death_col i>

The column with the year of death.

=item C<--sex_col i>

The column with sex_data.

=item C<--sex_chars f:m:u >

The chars for the sex in format f:m:u:

  --sex_chars F:M:?

=item C<--collapse>

Combines non-unique genotype ids to one genotype. The number of observed
genotypes is stored in the corresponding field. Ment for clonal organisms.

=item C<--man>

Display manpage.

=item C<--version>

Print version number of this software.

=back

=head1 DESCRIPTION

Converts a CSV file to the FRANz input format.  

=head1 CONFIGURATION AND ENVIRONMENT

C<csv.pl> does not support configuration files or environment variables.

=head1 DEPENDENCIES

L<Getopt::Long>, L<Text::CSV>

=head1 BUGS AND LIMITATIONS

No bugs have been reported. 

Please report any bugs or feature requests to
C<bug-latex-table@rt.cpan.org>, or through the web interface at
L<http://rt.cpan.org>. 

=head1 AUTHOR

Markus Riester  C<< <mriester@gmx.de> >>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2008, Markus Riester C<< <markus@bioinf.uni-leipzig.de> >>. 

This module is free software; you can redistribute it and/or
modify it under the same terms as Perl itself. See L<perlartistic>.

=head1 DISCLAIMER OF WARRANTY

BECAUSE THIS SOFTWARE IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY
FOR THE SOFTWARE, TO THE EXTENT PERMITTED BY APPLICABLE LAW. EXCEPT WHEN
OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES
PROVIDE THE SOFTWARE "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER
EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE
ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE SOFTWARE IS WITH
YOU. SHOULD THE SOFTWARE PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL
NECESSARY SERVICING, REPAIR, OR CORRECTION.

IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED TO IN WRITING
WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY AND/OR
REDISTRIBUTE THE SOFTWARE AS PERMITTED BY THE ABOVE LICENCE, BE
LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE
THE SOFTWARE (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A
FAILURE OF THE SOFTWARE TO OPERATE WITH ANY OTHER SOFTWARE), EVEN IF
SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGES.

=cut

# vim: ft=perl sw=4 ts=4 expandtab
