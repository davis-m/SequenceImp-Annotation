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

### This program accepts downloaded miRBase gff files and parses them to retrieve information required.
### For Genomic Ranges objects. This info includes hairpin name, miRBase ID  and genome coordinates.
### 040714 - The system was updated to require miRBase GFF3 files as GFF2 annotation has been discontinued in database.

use strict;
use warnings;
use Getopt::Long;

my $build = "";
my $format = 0;
my $help = 0;
my $miRNA_gff = "";
my $miRNA_fasta = "";

my %precursor_hash;

if(! GetOptions(
      "gff3=s" => \$miRNA_gff,
      "fasta=s" => \$miRNA_fasta,
      "help" => \$help
   )
){
   die "Could not parse the command line options\n";
}

if($help){
   print <<EOH;
   INPUT
   --gff3   The miRBase GFF3 file for the species
   --fasta  The corresponding file of miRNA FASTA sequences
   
   OUTPUT
   Output is printed to standard out
EOH
   exit(0);
}

# Check options
die "Need a miRNA gff3 file (See --help)\n" if length($miRNA_gff)==0;
die "Need a miRNA fasta file (See --help)\n" if length($miRNA_fasta)==0;

die "miRNA gff3 file ($miRNA_gff) does not exist\n" if (! -s $miRNA_gff);
die "miRNA fasta file ($miRNA_fasta) does not exist\n" if (! -s $miRNA_fasta);

# Read lines from files

open(GFF,"< $miRNA_gff") || die "Could not open the GFF3 file\n";
my @lines = <GFF>;
close(GFF) || die "Could not close the GFF3 file\n";

open(FASTA, "gzip -cd $miRNA_fasta |") || die "Could not open the miRBase fasta file\n";
my $line_count = 0;
my %fasta_seqs;
my $miRNA_id = "";
while(<FASTA>){
   $line_count++;
   chomp;
   my $this_line = $_;
   if(($line_count%2)==1){
      if($this_line =~ /^>[\w\-\.\?]+\s+(MIMAT\d+)/){
         $miRNA_id = $1;
      }else{
         die "Could not find the miRBase unique ID: $this_line\n";
      }
   }else{
      if($this_line =~ /^[UGCA]+$/){
         die "Same miRBase ID used twice: $miRNA_id\n" if (exists($fasta_seqs{$miRNA_id}));
         $fasta_seqs{$miRNA_id} = $this_line;
      }else{
         die "Unexpected symbol in sequence: \n";
      }
   }
}
close(FASTA) || die "Could not close the fasta file\n";

#use Data::Dumper;
#print Dumper \%fasta_seqs;
#exit (0);
###

my @comments = grep { /^#/ } @lines;
my @other    = grep { !/#/ } @lines;

print STDERR "---\n";
print STDERR @comments;
print STDERR "---\n";
print STDERR @other[0..5];
print STDERR "---\n";

my $buildCount = 0;

for my $caption ( @comments ){ 
    if ($caption =~ /gff-version\s+(\d+)/){
       $format = $1;
    }
    if ($caption =~ /^\#\s+genome-build-id[:]*\s+(\S+)/ || $caption =~ /^\#\s+Genome assembly:\s+(\S+)/ ){
       $build = $1;
       $buildCount++;
    }
    next;
}

die "ERROR: miRBase annotation does not appear to be gff3 format\n" if ($format !=3);
die "ERROR: Unclear definition of genome build in gff file\n" if $buildCount != 1;

print STDERR "miRBase genome build: $build\n";

print "ID\tName\tFeature\tChr\tStrand\tStart\tEnd\tPrecursor\tSequence\n";

my %recordcounts = map { my $x =()= /\t/g; ($x, 1) } @other;
die "record count conflict" unless keys (%recordcounts) == 1;

#use Data::Dumper;
#print Dumper \%recordcounts;
#exit(0);

for my $precursor_search(@other){
   
   my @splitSearch = split /\t/, $precursor_search;
   die "Expecting 9 fields per line in GFF3\n" if (scalar(@splitSearch) != 9);
   next if ($splitSearch[2] ne "miRNA_primary_transcript");

   if($splitSearch[8] =~ /Name=([\w\-\d]+)/){
      my $precursor_id = $1;
      if($splitSearch[8] =~ /ID=(MI\d+);/){
         $precursor_hash{$1} = $precursor_id;
      }else{
         die "Could not find the miRBase ID for $precursor_id\n";
      }
   }else{
      die "Could not find the precursor ID\n";
   };
}

for my $current (@other){
   my $ID="";
   my $chr="";
   my $strand="";
   my $start="";
   my $end="";
   my $feature="";
   my $hairpin = "";
   my $sequence = "";
   my $name = "";

   my @splitLine = split /\t/, $current;
   next if ($splitLine[2] ne "miRNA");
   my @IDlong = grep { /Name=/ } @splitLine;
   die "Too many IDs: @IDlong\n" if scalar @IDlong !=1;
   #print "$IDlong[0]\n";
   if ($IDlong[0] =~ /Name=([A-Za-z0-9\-\.]+);/){
      $name = $1;
   }else{
      die "Can't find an appropriate miRNA name within: @splitLine\n";
   }
  
   if($IDlong[0] =~ /Derives_from=(MI\d+)/){
      die "Precursor not found in gff3 file for $1\n" if (! exists($precursor_hash{$1}));
      $hairpin = $precursor_hash{$1};
   }else{
      die "Can't find the precursor ID link in:".$IDlong[0];
   }

   if($IDlong[0] =~ /ID=(MIMAT\d+)/){
      die "miRBase ID not found in fasta for $1\n" if (! exists($fasta_seqs{$1}));
      $sequence = $fasta_seqs{$1};
      $ID = $1;
   }else{
      die "Can't find the precursor ID link in:".$IDlong[0];
   }
   
   $chr = $splitLine[0];
   $chr=~s/^Chr(\w+)/$1/i;
   $feature = $splitLine[2];
   $start = $splitLine[3];
   $end = $splitLine[4];
   $strand = $splitLine[6];

   print "$ID\t$name\t$feature\t$chr\t$strand\t$start\t$end\t$hairpin\t$sequence\n";
}
