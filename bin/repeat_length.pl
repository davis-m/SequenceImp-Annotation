#!/usr/bin/env perl

###########################################################################
#
##This annotation pipeline is part of the Kraken framework which aims to
##facilitate RNA sequence analysis in a streamlined and efficient manner.
##Copyright (C) 2011 2012 2013 2014 EMBL - European Bioinformatics Institute 
#
##This program is free software: you can redistribute it and/or modify
##it under the terms of the GNU General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
#
##This program is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU General Public License for more details.
#
##You should have received a copy of the GNU General Public License
##along with this program (Please see the COPYING file for details).
##If not, see <http://www.gnu.org/licenses/>.
#
##Please send bug reports to:
##kraken@ebi.ac.uk
#
############################################################################

use warnings;
use strict;

my $inFile = "";
my $outFile = "";
my $baseNo = 0 ;
my $repeatID = "";

die "This script needs 2 argument" unless @ARGV == 2;

$inFile = shift;
$outFile = shift;

die "Need a specified input file: --inFile\n" if length($inFile)==0;
die "There is a problem with in input file" if (! -s $inFile);
die "Need a specified output file: --outFile\n" if length($outFile)==0;

open (REPEAT, "< $inFile");

my $Firstline = <REPEAT>;
chomp $Firstline;
#print "$Firstline\n";

die "The Fasta file seems to be incorrectly formatted, please check.\n" if !($Firstline =~ /^>/);

if ($Firstline =~ /^>(\S+)/){
   $repeatID = $1;
}

die "No repeat ID found at start of Fasta file\n" if length($repeatID) == 0;

while(<REPEAT>){
   chomp;
   die "There are multiple sequences in the Fasta file\n" if ($_ =~ /^>/);
   $baseNo += length($_);
}

open(OUTPUT, "> $outFile");

print OUTPUT "Chr\tLength\n";
print OUTPUT "$repeatID\t$baseNo\n";

close REPEAT;
close OUTPUT;

#print "$baseNo\n";
