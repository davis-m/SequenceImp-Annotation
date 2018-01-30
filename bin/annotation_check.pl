#!/usr/local/bin/perl

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

# This script will take the annotation files past by the user
# and perform simple checks to ensure that the files correspond
# to each other where appropriate.

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $chrClassFile = "";
my $chrLenFile = "";
my $tRNAFile = "";
my $miRFastaFile = "";
my $miRGFFFile = "";
my $repFile = "";
my $genesFile = "";
my $genomeFasta = "";
my $ncbi_seqs = "";
my $skip_gtf_chr = 0;
my %failed;
my $debug = "";

print STDERR "\n### Checking annotation file formats ###\n";

if(! GetOptions(
   "chrlen=s" => \$chrLenFile,
   "chrclass=s" => \$chrClassFile,
   "fasta=s"   => \$genomeFasta,
   "tRNAs=s" => \$tRNAFile,
   "miR_fasta=s" => \$miRFastaFile,
   "miR_gff=s" => \$miRGFFFile,
   "reps=s" => \$repFile,
   "genes=s" => \$genesFile,
   "ncbi=s" => \$ncbi_seqs,
   "no_gtf_chr_check" => \$skip_gtf_chr,
   "debug"  => \$debug
)){
   die "Could not parse the command line options in annotation_check.pl\n";
}

die "Need a chromosome length file\n" if length($chrLenFile) == 0;
die "Need a chromosome class file\n" if length($chrClassFile) == 0;
#die "Need a tRNA file\n" if $tRNAFile ==0
#die "Need a repeat file\n" if $repFile ==0
die "Need a miRBase Fasta file\n" if length($miRFastaFile) ==0;
die "Need a miRBase GFF 3 file\n" if length($miRGFFFile) ==0;
die "Need a gene GTF file\n" if length($genesFile) ==0;
die "Need a genome FASTA file\n" if length($genomeFasta) ==0;

die "$genomeFasta not found or empty\n" if (! (-s $genomeFasta));
die "$chrLenFile not found or empty\n" if (! (-s $chrLenFile));
die "$chrClassFile not found or empty\n" if (! (-s $chrClassFile));
die "$miRFastaFile not found or empty\n" if (! (-s $miRFastaFile));
die "$miRGFFFile not found or empty\n" if (! (-s $miRGFFFile));
die "$genesFile not found or empty\n" if (! (-s $genesFile));

die "$repFile not found or empty\n" if ((length($repFile)!=0) && (! (-s $repFile)));
die "$tRNAFile not found or empty\n" if ((length($tRNAFile)!=0) && (! (-s $tRNAFile)));

my @all_ncbi;
if(length($ncbi_seqs)!=0){
   @all_ncbi = split ",", $ncbi_seqs;
   foreach(@all_ncbi){
      die "$_ not found or empty\n" if (! (-s $_))
   }
}

# OK

# Read in the chromosome length and chromosome class files into a hash for comparisons

# Collecting chromosome classes
my %chromosome_class;
my %class_count;

open (CLASS, "< $chrClassFile") || die "Could not open $chrClassFile\n";
my $class_header = <CLASS>;
chomp $class_header;
die "Class file header not as expected: $class_header\n" if(!($class_header=~/Chr\tClass/));

while(<CLASS>){
   chomp;
   my($key,$value) = split "\t";
   $chromosome_class{$key} = $value;
   $class_count{$key}++;
}

close(CLASS) || die "Could not close the chromosome class file\n";

print STDERR (Dumper(\%chromosome_class)) if $debug;

##### Check chromosome classes - avoid patches
my @class_crunch;
for (keys %class_count){
   if ($chromosome_class{$_} ne "chromosome" && $chromosome_class{$_} ne "scaffold"){
      push @class_crunch, $chromosome_class{$_};
   }
}

print STDERR "\nCHECKING Genome FASTA classifications\n";
scalar(@class_crunch) != 0 ? $failed{"genome_class"} = "FAIL" : print STDERR "Genome FASTA: Chromosome class - PASS\n";

# Collecting chromosome lengths
my %chromosome_length;

open (LEN, "< $chrLenFile") || die "Could not open $chrLenFile\n";
my $length_header = <LEN>;
chomp $length_header;
die "Length file header not as expected: $length_header\n" if(!($length_header=~/Chr\tLength/));

while(<LEN>){
   chomp;
   my($key,$value) = split "\t";
   $chromosome_length{$key} = $value;
}
close(LEN) || die "Could not close the chromosome length file\n";

print STDERR (Dumper(\%chromosome_length)) if $debug;

# Collecting genome version - eg. data/drosophila_short/Drosophila_melanogaster.BDGP5.dna_chr4.toplevel.fa.gz
my $Ensembl_version = "";
my $clipped_genomeFasta = $genomeFasta;
$clipped_genomeFasta =~ s/^.*\/([^\/]+$)/$1/;
if($clipped_genomeFasta =~ /^\w+\.(\w+)\.\w+\.\w+\.\w+(\.gz){0,1}+$/){
   $Ensembl_version = $1;
}

# Most significant tests:
# Genome version where possible
# Chromosome names in files
# Chromosome classes in files
# Coordinates within expected chromosome ranges

# Check the genes Ensembl gtf - eg. data/drosophila_short/Drosophila_melanogaster.BDGP5.78_chr4.gtf.gz
print STDERR "\nCHECKING Ensembl GTF: $genesFile\n";

sub check_a_gtf { 
  
   my %temp_gene_chr_crunch;
   my @temp_gene_coord_crunch;
   my $found_vers = "";
   my $origin = "NOPE";

   die "&check_a_gtf requires 3 options\n" if scalar(@_) != 3;

   my $aGTF = shift;
   $origin = shift;
   my $length_ref = shift;
   my %chromosome_length_temp = %{$length_ref};

   die "GTF/GFF origin must be either 'miRBase' OR 'Ensembl'\n" if (($origin ne 'Ens') && ($origin ne 'miR'));

   if($origin eq 'Ens'){

      my $slicedGTFFile = $aGTF;
      
      $slicedGTFFile =~ s/^.*\/([^\/]+$)/$1/;
      if($slicedGTFFile =~ /^\w+\.(\w+)\.\w+\.gtf(\.gz){0,1}+$/){
         $found_vers = $1;
      }
   }

   my $geneLineNo = 0;
   
   open (GENES, "gzip -cdf $aGTF |") || die "Could not open the genes GTF file\n";
   while(<GENES>){
      chomp;
      $geneLineNo++; 
      
      if(/^#/){
         if($origin eq 'miR'){
            if(/^#\s+genome-build-id:\s+([^\.]+)(\.){0,1}\d*$/){
               $found_vers = $1;
               next;
            }else{
               next; 
            }
         }else{
            next;
         }
      }

      my @gene_line = split "\t";
      die "ERROR: $aGTF line $geneLineNo - Requires 9 fields\n" if scalar (@gene_line) != 9;
      die "ERROR: $aGTF line $geneLineNo - 4th field should be a coordinate\n" if (! $gene_line[3] =~ /^\d+$/);
      die "ERROR: $aGTF line $geneLineNo - 5th field should be a coordinate\n" if (! $gene_line[4] =~ /^\d+$/);
     

      if (! exists($chromosome_length_temp{$gene_line[0]})){
         $temp_gene_chr_crunch{$gene_line[0]}++;
      }else{
         if(($gene_line[3] < 0) | ($gene_line[4] > $chromosome_length_temp{$gene_line[0]})){
            push @temp_gene_coord_crunch, $geneLineNo;
         }
      }
   }
   close (GENES) || die "Could not close $aGTF\n";
  
   die "Haven't found a genome version for $aGTF\n" if length($found_vers) == 0;
   my @gtf_reply = (\%temp_gene_chr_crunch,\@temp_gene_coord_crunch, $found_vers);

   return (\@gtf_reply);
}

my @gtf_return_refs = @{&check_a_gtf($genesFile,"Ens",\%chromosome_length)};

my %gene_chr_crunch = %{$gtf_return_refs[0]};
my @gene_coord_crunch = @{$gtf_return_refs[1]};
my $Genes_GTF_version = $gtf_return_refs[2];

($Ensembl_version ne $Genes_GTF_version) ? $failed{"gtf_genome"} = "FAIL" : print STDERR "Ensembl GTF: Genome version - PASS\n";
if(! $skip_gtf_chr){
   (scalar(keys(%gene_chr_crunch))>0) ? $failed{"gtf_chr"} = "FAIL" : print STDERR "Ensembl GTF: Chromosome check - PASS\n";
}else{
   print STDERR "Ensembl GTF: Chromosome check - SKIPPED\n";
}
(scalar(@gene_coord_crunch)>0) ? $failed{"gtf_coords"} = "FAIL" : print STDERR "Ensembl GTF: Coordinate bounds - PASS\n";

# Check the miRNA miRBase gtf - eg. data/drosophila_short/dme_short_chr4.gff3

print STDERR "\nCHECKING miRBase GFF: $miRGFFFile\n";

my @miR_return_refs = @{&check_a_gtf($miRGFFFile,"miR",\%chromosome_length)}; 

my %miR_chr_crunch = %{$miR_return_refs[0]};
my @miR_coord_crunch = @{$miR_return_refs[1]};
my $miR_GFF_version = $miR_return_refs[2];

($Ensembl_version ne $miR_GFF_version) ? $failed{"miR_genome"} = "FAIL" : print STDERR "miRBase GFF: Genome version - PASS\n";
(scalar(keys(%miR_chr_crunch))>0) ?  $failed{"miR_chr"} = "FAIL": print STDERR "miRBase GFF: Chromosome check - PASS\n";
(scalar(@miR_coord_crunch)>0) ? $failed{"miR_coords"} = "FAIL" : print STDERR "miRBase GFF: Coordinate bounds - PASS\n";

# Check the miRNA miRBase FASTA file - eg. data/drosophila/mature.fa.gz

print STDERR "\nCHECKING miRBase FASTA: $miRFastaFile\n";

open (MIRFAST, "gzip -cdf $miRFastaFile |") || die "Could not open $miRFastaFile\n";
my $miRline = 1;

my @miR_wrong_format;
my @miR_wrong_seq;
my @miR_wrong_length;

while(<MIRFAST>){
   chomp;
   if (($miRline % 2)==1){
      if (! /^>/){
         push @miR_wrong_format, $miRline;
      }
   }else{
      if (! /^[AUGC]+$/){
         push @miR_wrong_seq, $miRline;
      }

      if (length($_) > 40){
         push @miR_wrong_length, $miRline;
      }
   }
   $miRline++;
}
close (MIRFAST) || die "Could not close $miRFastaFile\n";
      
(scalar(@miR_wrong_format) != 0) ? $failed{"miR_fasta_format"} = "FAIL" : print STDERR "miRBase FASTA: Expected IDs - PASS\n";   
(scalar(@miR_wrong_seq) != 0) ? $failed{"miR_seqeunce"} = "FAIL" : print STDERR "miRBase FASTA: Expected nucleotides - PASS\n";
(scalar(@miR_wrong_length) != 0) ? $failed{"miR_length"} = "FAIL" : print STDERR "miRBase FASTA: Mature miRNA length - PASS\n";


# Check the repeat annotation file - eg. data/drosophila_short/dm3_chr4.fa.out

my @rep_field_crunch;
my %rep_chr_crunch;
my @rep_coord_crunch;

if (length($repFile) >0){
   print STDERR "\nCHECKING RepeatMasker File: $repFile\n";
   open (REPS, "gzip -cdf $repFile |") || die "Could not open $repFile\n";
   my $header_one = <REPS>;
   my $header_two = <REPS>;
   my $header_three = <REPS>;
  
   if (
      (! ($header_one   =~ /^\s*SW\s+perc\s+perc\s+perc\s+query\s+position\s+in\s+query\s+matching\s+repeat\s+position\s+in\s+repeat$/)) ||
      (! ($header_two   =~ /^\s*score\s+div\.\s+del\.\s+ins\.\s+sequence\s+begin\s+end\s+\(left\)\s+repeat\s+class\/family\s+begin\s+end\s+\(left\)\s+ID$/)) ||
      (! ($header_three =~ /^\s+$/))
   ){
      $failed{"rep_header"} = "FAIL";
   }else{
      print STDERR "RepeatMasker File: Expected header - PASS\n";
   }
   
   my @row_split;
   my $rep_line_no = 4;
   while(<REPS>){
      chomp;
      my $rep_line = $_;
      $rep_line =~ s/^\s+//;
      @row_split = split " {1,}", $rep_line, -1;
      pop @row_split if ($row_split[-1] eq "*");
      #print scalar(@row_split)."\n";
      die "ERROR: $repFile line $rep_line_no - 6th field should be a coordinate (expect space delineation)\n" if (! $row_split[5] =~ /^\d+$/);
      die "ERROR: $repFile line $rep_line_no - 7th field should be a coordinate (expect space delineation)\n" if (! $row_split[6] =~ /^\d+$/);
     
      if(scalar(@row_split) != 15){
         push @rep_field_crunch, $rep_line_no;
      }else{
         if(! exists($chromosome_length{$row_split[4]})){
            $rep_chr_crunch{$row_split[4]}++;
         }else{
            if (($row_split[5] < 1) || ($row_split[6] > $chromosome_length{$row_split[4]})){
               push @rep_coord_crunch , $rep_line_no;
            }
         }
      }
      $rep_line_no++;
   }
   
   close (REPS) || die "Could not close $repFile\n";
   
   (scalar(keys(%rep_chr_crunch)) > 0) ? $failed{"rep_chr"} = "FAIL" : print STDERR "RepeatMasker File: Chromosome check - PASS\n";
   (scalar(@rep_field_crunch) > 0) ? $failed{"rep_fields"} = "FAIL" :  print STDERR "RepeatMasker File: Field number - PASS\n";
   (scalar(@rep_coord_crunch) > 0) ? $failed{"rep_coords"} = "FAIL" :  print STDERR "RepeatMasker File: Coordinate bounds - PASS\n";
}

# Check the tRNA annotation file - eg. data/drosophila/Drosophila_melanogaster.BDGP5.78.out 
   
my @tRNA_field_crunch;
my %tRNA_chr_crunch;
my @tRNA_coord_crunch;
my @tRNA_header_crunch;

if (length($tRNAFile) >0){
   print STDERR "\nCHECKING tRNAScan-SE File: $tRNAFile\n";
   my $tRNA_line_no = 0;
   
   open (TRNA,"gzip -cdf $tRNAFile |") || die "Could not open $tRNAFile\n";
   while(<TRNA>){
      my $fresh_line = $_;
      $tRNA_line_no++;
      
      if (($fresh_line =~ /^Seq/)||($fresh_line =~ /^Name/)||($fresh_line =~ /^---/)){
         push @tRNA_header_crunch, $tRNA_line_no;
         next;
      }

      $fresh_line=~s/\s+/\t/g;
      my @tRNA_line = split "\t",$fresh_line;
      #print scalar(@tRNA_line)."\n";
   
      if(scalar(@tRNA_line)!=9){
         push @tRNA_field_crunch, $tRNA_line_no;
      }else{
         die "ERROR: $tRNAFile line $tRNA_line_no - 3rd field should be a coordinate\n" if (! ($tRNA_line[2] =~ /^\d+$/));
         die "ERROR: $tRNAFile line $tRNA_line_no - 4th field should be a coordinate\n" if (! ($tRNA_line[3] =~ /^\d+$/));
         if (! exists($chromosome_length{$tRNA_line[0]})){
            $tRNA_chr_crunch{$tRNA_line[0]}++;
         }else{
            push @tRNA_coord_crunch, $tRNA_line_no if ($tRNA_line[2] < 1 || $tRNA_line[3] < 1 || $tRNA_line[2] > $chromosome_length{$tRNA_line[0]} || $tRNA_line[3] > $chromosome_length{$tRNA_line[0]});  
         }
      }
      
   }
   close (TRNA) || die "Could not close $tRNAFile\n";
   
   (scalar(@tRNA_header_crunch) > 0 ) ? $failed{"tRNA_header"} = "FAIL" : print STDERR "tRNAScan-SE File: Header absent - PASS\n";
   (scalar(@tRNA_field_crunch) > 0 ) ? $failed{"tRNA_fields"} = "FAIL" : print STDERR "tRNAScan-SE File: Field number - PASS\n";
   (scalar(keys(%tRNA_chr_crunch)) > 0 ) ? $failed{"tRNA_chr"} = "FAIL" : print STDERR "tRNAScan-SE File: Chromosome check - PASS\n";
   (scalar(@tRNA_coord_crunch) > 0) ? $failed{"tRNA_coords"} = "FAIL" : print STDERR "tRNAScan-SE File: Coordinate bounds - PASS\n";
}

# Check repeat element ncbi fasta files eg. data/drosophila_chr3R/copia_element.fa
my @ncbi_id_crunch;
my %ncbi_multi_crunch;

if (scalar(@all_ncbi)!=0){
   print STDERR "\nCHECKING NCBI FASTA files: $ncbi_seqs\n";
   for my $check_ncbi (@all_ncbi){
      open(NCBI, "< $check_ncbi") || die "Could not open $check_ncbi\n";
      my $ncbi_id = <NCBI>;
      push @ncbi_id_crunch , $check_ncbi if(! ($ncbi_id=~/^>/));
      while(<NCBI>){
         chomp;
         my $ncbi_line = $_;
         $ncbi_multi_crunch{$check_ncbi}++ if ($ncbi_line =~ /^>/);
      }
      close(NCBI) || die "Could not close $check_ncbi\n";
   }
   
   (scalar(@ncbi_id_crunch) > 0) ? $failed{"ncbi_id"} = "FAIL" : print "NCBI FASTA files: ID format - PASS\n"; 
   (scalar(keys(%ncbi_multi_crunch)) > 0) ? $failed{"ncbi_singleton"} = "FAIL" : print "NCBI FASTA files: Single sequence - PASS\n"; 

}



# Report the datails of all the failed tests
#$tRNAFile
#$miRFastaFile
#$miRGFFFile 
#$repFile 
#$genesFile 
#$genomeFasta 

if (scalar(keys(%failed))==0){
   print STDERR "\nFormat tests completed\n\n";
   exit (0);
}else{
   print STDERR "\n** FAILED TESTS **";
   local $, = " ";
   print STDERR "\nERROR $genomeFasta : Contains unexpected sequence classes -",@class_crunch if exists($failed{"genome_class"});
   print STDERR "\nERROR $genesFile : Genome version does not match Ensembl genome file supplied - $Ensembl_version vs. $Genes_GTF_version" if exists($failed{"gtf_genome"});
   print STDERR "\nERROR $genesFile : Chromosome sequences not found in genomic FASTA -",keys(%gene_chr_crunch) if exists($failed{"gtf_chr"});
   print STDERR "\nERROR $genesFile : Coordinates stray beyond chromosome boundaries - eg. line",$gene_coord_crunch[0] if exists($failed{"gtf_coords"});
   print STDERR "\nERROR $miRGFFFile : Genome version does not match Ensembl genome file supplied - $Ensembl_version vs. $miR_GFF_version" if exists($failed{"miR_genome"});
   print STDERR "\nERROR $miRGFFFile : Chromosome sequences not found in genomic FASTA -",keys(%miR_chr_crunch) if exists($failed{"miR_chr"});
   print STDERR "\nERROR $miRGFFFile : Coordinates stray beyond chromosome boundaries - eg. line",$miR_coord_crunch[0] if exists($failed{"miR_coords"});
   print STDERR "\nERROR $miRFastaFile : FASTA ID lines do not conform to FASTA format - eg. line",$miR_wrong_format[0] if exists($failed{"miR_fasta_format"});
   print STDERR "\nERROR $miRFastaFile : FASTA sequences contain unexpected nucleotides - eg. line",$miR_wrong_seq[0] if exists($failed{"miR_seqeunce"});
   print STDERR "\nERROR $miRFastaFile : Mature miRNAs are longer than expected - eg. line",$miR_wrong_length[0] if exists($failed{"miR_length"});
   print STDERR "\nERROR $repFile : RepeatMasker header not as expected" if exists($failed{"rep_header"});
   print STDERR "\nERROR $repFile : Chromosome sequences not found in genomic FASTA -",keys(%rep_chr_crunch) if exists($failed{"rep_chr"});
   print STDERR "\nERROR $repFile : Unexpected number of fields in a line - eg. line",$rep_field_crunch[0] if exists($failed{"rep_fields"});
   print STDERR "\nERROR $repFile : Coordinates stray beyond chromosome boundaries - eg. line",$rep_coord_crunch[0] if exists($failed{"rep_coords"});
   print STDERR "\nERROR $tRNAFile : Lines appear to resemble tRNAscan-SE header (run tRNAscan with -b specified) - eg. line",$tRNA_header_crunch[0] if exists($failed{"tRNA_header"});
   print STDERR "\nERROR $tRNAFile : Unexpected number of fields in a line (expect 9) - eg. line",$tRNA_field_crunch[0] if exists($failed{"tRNA_fields"});
   print STDERR "\nERROR $tRNAFile : Chromosome sequences not found in genomic FASTA -",keys(%tRNA_chr_crunch) if exists($failed{"tRNA_chr"});
   print STDERR "\nERROR $tRNAFile : Coordinates stray beyond chromosome boundaries - eg. line",@tRNA_coord_crunch if exists($failed{"tRNA_coords"});
   print STDERR "\nERROR NCBI FASTA files: FASTA ID expected in the first line of the file (check file is not compressed) -",@ncbi_id_crunch if exists($failed{"ncbi_id"});
   print STDERR "\nERROR NCBI FASTA files: Single sequence expected in each FASTA file (check for multiple IDs) -",keys(%ncbi_multi_crunch) if exists($failed{"ncbi_singleton"});
   print STDERR "\n\nPlease check file requirements and formats and resubmit\n\n";

   exit (1);
}

# my @ncbi_id_crunch;
# my %ncbi_multi_crunch;


#print Dumper \@class_crunch;
#print STDERR Dumper \@rep_field_crunch if $debug;
#print STDERR Dumper keys %rep_chr_crunch if $debug;


