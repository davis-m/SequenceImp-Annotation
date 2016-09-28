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

### This script takes the standard output from RepeatMasker and reformats it to
# a format suitable for pipline R scripts (tab delineated, same column numbers, remove
# leading space).

my $line_number;

while(<>){
   $line_number++;
   if($line_number > 3){
      chomp;
      s/\s+/\t/g;
      s/^\t//;
      s/\t$/\tNA/;
      print "$_\n";
   }
}

