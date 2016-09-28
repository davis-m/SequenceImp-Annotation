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

# This script will take a Ensembl GTF file and translate it into a set of
# introns, exons and ORFs for the annotation builder.

# TO DO:
# Add additional file format tests
# Check CDS and stop codon align (beware inter-exon boundaries inserting distance)
# Check sequential exon numbers

use strict;
use warnings;
use Getopt::Long;
#use Data::Dumper;

my $gtfFile = "";
my $outFile = "";
my $help = 0;
my %parsed_gtf;
my %class_info;
my %gene_info;
my %strand_info;
my %chr_info;
my $no_start_example = "";
my $no_stop_example = "";
my %feature_summary;
my $length_file = "";
my %chromosomes;

sub help {
   print <<EOH;
   --gtf=<FILE>   REQUIRED    The Ensembl GTF file to be parsed
   --len=<FILE>   REQUIRED    A chromosome length file to filter annotation to remove 
                              undesirable annotation.
   --help         OPTIONAL    Prints this page
EOH
}

if(!GetOptions(
      "gtf=s" => \$gtfFile,
      "len=s" => \$length_file,
      "help"  => \$help
   )){
   die "Option parsing has failed\n";
}

if($help){
   help();
   exit(0);
}

die "Need a GTF file (--gtf=<FILE>)\n" if(length($gtfFile)==0);
die "GTF file is emtpy or absent\n" if (! -s $gtfFile);
die "Need a chromosome length file (--len=<FILE>)\n" if(length($length_file)==0);
die "Length file is emtpy or absent\n" if (! -s $length_file);

open(LEN, "< $length_file") || die "Could not open the length file\n";
my $length_header = <LEN>;
while(<LEN>){
   chomp;
   my @len_info = split "\t", $_;
   die "Same chromosome twice in $length_file\n" if exists $chromosomes{$len_info[0]};
   $chromosomes{$len_info[0]}++;
}
close(LEN) || die "Could not close the length file\n";


open(GTF,"gzip -cd $gtfFile |") || die "Could not open the GTF file\n";

my $gtf_line;
my %extra_chr;
my $new_class_check = 0;
my $new_class_system = 0;

print STDERR "Parsing GTF file\n";
while(<GTF>){
   $gtf_line++;
   print STDERR "$gtf_line lines\n" if(! ($gtf_line%50000));
   chomp;

   if(/^#/){
      print STDERR "Skipping comment line: $_\n";
      next;
   }

   my @split_line = split "\t", $_;
   die "Expecting 9 fields per line in GTF: $gtf_line\n" if (scalar(@split_line) != 9);
   my($chr,$class,$feature,$start,$end,$score,$strand,$frame,$attribute) = @split_line;

   # Skip the gene lines in the GTF file. Focus on transcripts only.
   next if ($feature eq "gene");

   # Skip annotation for which chromosomal information is not available (no sequence in seq file)
   if (! exists($chromosomes{$chr})){
      $extra_chr{$chr}++;
      next;
   }

   my $gene_id = "";
   if ($attribute=~/.*gene_id\s+"([\w+\.-]+)";.*/){
      $gene_id = $1;
   }else{
      die "Could not find the gene ID: $gtf_line\n$attribute\n";
   }

   my $transcript_id = "";
   if($attribute=~/.*transcript_id\s+"([\w+\.-]+)";.*/){
      $transcript_id = $1; 
   }else{
      die "Could not find the transcript ID: $gtf_line\n$attribute\n";
   }
  
   # In later versions of Ensembl, position of transcript biotype info moved. Must work with both systems.
   if ($new_class_check == 0){
      $new_class_check = 1;
      if($attribute=~/.*transcript_biotype\s+"([\w+\.-]+)";.*/){
         $new_class_system = 1;
         $class=$1;
         print STDERR "Configuring pipeline to use 'transcript_biotype' attribute to locate transcript biotype\n";
      }else{
         print STDERR "Configuring pipeline to use 'source' field to locate transcript biotype\n";
      }
   }elsif($new_class_system > 0){
      if($attribute=~/.*transcript_biotype\s+"([\w+\.-]+)";.*/){
         $class=$1;   
      }else{
         die "Expecting all transcript/exon lines to have a transcript_biotype: $gtf_line\n$attribute\n";
      }
   }

   $chr_info{$chr}->{$transcript_id}++;
   die "Same feature seen twice: $transcript_id: $feature: $start\n" if exists($parsed_gtf{$transcript_id}{$feature}{$start});
   $parsed_gtf{$transcript_id}->{$feature}->{$start} = $end;
   $gene_info{$transcript_id}{$gene_id}++;
   $class_info{$transcript_id}{$class}++;
   $strand_info{$transcript_id}{$strand}++;
}
# Check strands
close (GTF) || die "Could not close the GTF file\n";

my @missing_chrs = map {(!exists($chr_info{$_})) ? $_ : ()} keys(%chromosomes);
{
   local $, = " ";
   print STDERR "Contigs missing from the annotation:",@missing_chrs,"\n";
}

if(scalar(keys %extra_chr)>0){
   local $, = " ";
   print STDERR "Extra chromosomes in annotation removed by filter:",keys %extra_chr,"\n";
}

my $CDSsum = 0;
my $no_start = 0;
my $no_stop = 0;


my %all_trans_per_chr;

print STDERR "\nArranging transcript annotation per chromosome\n";
for my $chromosome (keys %chr_info){
   print STDERR "Chromosome: $chromosome\n";
   for my $transcript (keys %{$chr_info{$chromosome}}){
     
      die "Same transcript on 2 chromosomes: $transcript\n" if exists($all_trans_per_chr{$transcript});
      $all_trans_per_chr{$transcript}++;

      # Check annotation consistency

      my @these_classes = keys(%{$class_info{$transcript}});
      die "Multiple classes for a single transcript\n" if (scalar(@these_classes)!=1);
      my $this_class = $these_classes[0];
      $feature_summary{$this_class}++;

      my @these_genes = keys(%{$gene_info{$transcript}});
      die "Multiple genes for a single transcript\n" if (scalar(@these_genes)!=1);
      my $this_gene = $these_genes[0];

      my @these_strands = keys(%{$strand_info{$transcript}});
      die "Multiple strands for a single transcript\n" if (scalar(@these_strands)!=1);
      my $this_strand = $these_strands[0];
      die "Unexpected strand definition: $transcript $this_strand\n" if ($this_strand ne "+" && $this_strand ne "-");


      # Strand specific UTR designation
      my $UTR_upstream = "5UTR";
      my $UTR_downstream = "3UTR";

      if($this_strand eq "-"){
         $UTR_upstream = "3UTR";
         $UTR_downstream = "5UTR";
      }

      # Sort out exons and introns
      my @ordered_exons = sort {$a <=> $b} keys %{$parsed_gtf{$transcript}{"exon"}};
      
      my $intron_start = -1;
      my $intron_end = -1;

      my $anyCDS = 0;
      my $CDSstart = -1;   # Assumes one CDS per transcript. Assume CDS annotation is a subset of exons.
      my $CDSend = -1;

      # Work out the start and stop of the ORF - Will depend on start and stop coordinates
      if(exists($parsed_gtf{$transcript}{"CDS"})){
         $anyCDS = 1;
         $CDSsum++;
         $no_start++ if (!exists($parsed_gtf{$transcript}{"start_codon"}));
         $no_start_example = $transcript if (!exists($parsed_gtf{$transcript}{"start_codon"}));
         $no_stop++ if (!exists($parsed_gtf{$transcript}{"stop_codon"}));
         $no_stop_example = $transcript if (!exists($parsed_gtf{$transcript}{"stop_codon"}));
         
#         print STDERR "$transcript:\n" if (!exists($parsed_gtf{$transcript}{"stop_codon"}));
#         print STDERR Dumper(\$parsed_gtf{$transcript}) if (!exists($parsed_gtf{$transcript}{"stop_codon"}));

         # Altered 010714 - Will use the start and stop codon coordinates for each ORF as stop codon is not included in the annotated CDS coordinates
         my @ordered_CDS_coords = sort {$a <=> $b} keys %{$parsed_gtf{$transcript}{"CDS"}};
#
#         $CDSstart = $ordered_CDS_coords[0];
#         $CDSend = $parsed_gtf{$transcript}{"CDS"}{$ordered_CDS_coords[-1]}[0];

#         die "CDS but no start codon: $transcript\n" if (! exists($parsed_gtf{$transcript}{"start_codon"}));  # Some ORFs are missing start codons
#         die "CDS but no stop codon: $transcript\n" if (! exists($parsed_gtf{$transcript}{"stop_codon"}));    # Some ORFs are missing stop codons
         
         #my @ordered_start_coords;
         #if (exists($parsed_gtf{$transcript}{"start_codon"})){
         #   @ordered_start_coords = sort {$a <=> $b} keys %{$parsed_gtf{$transcript}{"start_codon"}};
         #}
        

         # Calculate and report the number of ORFs with and without start and stop codons.
         
         die "CDS longer than the exonic annotation: $transcript\n" if (($ordered_CDS_coords[0] < $ordered_exons[0]) || ($parsed_gtf{$transcript}{"CDS"}{$ordered_CDS_coords[-1]} > $parsed_gtf{$transcript}{"exon"}{$ordered_exons[-1]}));

         # Does stop codon match ORF end?

         my @ordered_stop_coords;
         # Stop codon may span an exon-exon boundary

         if (exists($parsed_gtf{$transcript}{"stop_codon"})){
            @ordered_stop_coords = sort {$a <=> $b} keys %{$parsed_gtf{$transcript}{"stop_codon"}};
         }

         if ($this_strand eq "+"){
            $CDSstart = $ordered_CDS_coords[0];
            if (exists($parsed_gtf{$transcript}{"stop_codon"})){
               $CDSend = $parsed_gtf{$transcript}{"stop_codon"}{$ordered_stop_coords[-1]};
            }else{
               $CDSend = $parsed_gtf{$transcript}{"CDS"}{$ordered_CDS_coords[-1]};
            }
         }elsif($this_strand eq "-"){
            if (exists($parsed_gtf{$transcript}{"stop_codon"})){
               $CDSstart = $ordered_stop_coords[0];
            }else{
               $CDSstart = $ordered_CDS_coords[0];
            }
            $CDSend = $parsed_gtf{$transcript}{"CDS"}{$ordered_CDS_coords[-1]};
         }else{
            die "Strand unexpected: $transcript - $this_strand\n";
         }
      }

      for my $this_start (@ordered_exons){
        
         my $this_end = $parsed_gtf{$transcript}{"exon"}{$this_start};


         if($intron_start > 0){
            $intron_end = $this_start -1;
            print STDOUT "$chromosome\t$intron_start\t$intron_end\t",$this_strand,"\t$transcript\t",$this_class,"_intron\t",$this_gene,"\n";
         }
         $intron_start = $this_end +1; 
         if($anyCDS){
            if($this_start > $CDSend){ # 3'UTR
               print STDOUT "$chromosome\t",$this_start,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_$UTR_downstream\t",$this_gene,"\n"; 
            }elsif($this_end < $CDSstart){ # 5'UTR
               print STDOUT "$chromosome\t",$this_start,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_$UTR_upstream\t",$this_gene,"\n"; 
            }elsif($this_start >= $CDSstart && $this_end <= $CDSend){ # ORF                                    
               print STDOUT "$chromosome\t",$this_start,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_ORF\t",$this_gene,"\n"; 
            }elsif($this_start < $CDSstart && $this_end <= $CDSend && $this_end >= $CDSstart){ # 5'-ORF
               print STDOUT "$chromosome\t",$CDSstart,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_ORF\t",$this_gene,"\n"; 
               print STDOUT "$chromosome\t",$this_start,"\t",$CDSstart-1,"\t",$this_strand,"\t$transcript\t",$this_class,"_$UTR_upstream\t",$this_gene,"\n"; 
            }elsif($this_start >= $CDSstart && $this_start <= $CDSend && $this_end > $CDSend){ # 3'-ORF
               print STDOUT "$chromosome\t",$this_start,"\t",$CDSend,"\t",$this_strand,"\t$transcript\t",$this_class,"_ORF\t",$this_gene,"\n"; 
               print STDOUT "$chromosome\t",$CDSend+1,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_$UTR_downstream\t",$this_gene,"\n"; 
            }elsif($this_start < $CDSstart && $this_end > $CDSend){ # Single exon gene
               print STDOUT "$chromosome\t",$this_start,"\t",$CDSstart -1,"\t",$this_strand,"\t$transcript\t",$this_class,"_$UTR_upstream\t",$this_gene,"\n"; 
               print STDOUT "$chromosome\t",$CDSstart,"\t",$CDSend,"\t",$this_strand,"\t$transcript\t",$this_class,"_ORF\t",$this_gene,"\n"; 
               print STDOUT "$chromosome\t",$CDSend+1,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_$UTR_downstream\t",$this_gene,"\n"; 
            }else{
               die "Unexpected ORF-to-exon combination\n";
            }
         }else{
            print STDOUT "$chromosome\t",$this_start,"\t",$this_end,"\t",$this_strand,"\t$transcript\t",$this_class,"_exon\t",$this_gene,"\n"; 
         }
      }
   }
}

# NOTE: The GTF files from Ensembl do not necessarily seem to annotate the predicted intron sequences in tRNAs.

# Transcript summary
print STDERR "\nTranscripts classes identified\n";
print STDERR "$_\t$feature_summary{$_}\n" for (keys %feature_summary);
my $total_transcripts=0;
$total_transcripts+=$feature_summary{$_} for (keys %feature_summary);
print STDERR "Total\t$total_transcripts\n";

# CDS summary
print STDERR "\nCDSs:\t$CDSsum\nCDSs with no start codon:\t$no_start (eg. $no_start_example)\nCDSs with no stop codon:\t$no_stop (eg. $no_stop_example)\n";

print STDERR "GTF processing completed\n";



