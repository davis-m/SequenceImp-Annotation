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

use strict;
use warnings;

use Getopt::Long;

BEGIN {
   my $path = $0;
   $path =~ s|[^/]*$||;
   if ($path =~ /\/?bin\/$/){
      $path =~ s/\/?bin\/$//;
   }else{
      $path =~ s/\/$//;
   }
   #print "$path\n";
   $ENV{IMPANNOT_ROOT} = $path;
}

## How to check it is a fasta file?
# Use gff2 not gff3 -  Need to discriminate between the 2
# Add individual NCBI sequences to a previous annotation
# Check hierarchy. Duplicates?
# Report overlaps in miRNA annotation

print STDERR "\n### Welcome to the SequenceImp Annotation pipeline ###\n\n";

my $miRBase_gff_pl = "$ENV{IMPANNOT_ROOT}/bin/miRBase_gff3_parse.pl";
my $miRBase_GR_R = "$ENV{IMPANNOT_ROOT}/bin/miRBase_GR_convert.R";
my $repFormat_pl = "$ENV{IMPANNOT_ROOT}/bin/repMasker_format.pl";
my $gtfTranslate_pl = "$ENV{IMPANNOT_ROOT}/bin/gtf_translate.pl";
my $genomic_annot_conv_R = "$ENV{IMPANNOT_ROOT}/bin/genomic_annotation_GR_convert.R";
my $repeat_length_pl = "$ENV{IMPANNOT_ROOT}/bin/repeat_length.pl";
my $repeat_setup_sh = "$ENV{IMPANNOT_ROOT}/bin/repeat_setup.sh";
my $default_hierarchy_file = "$ENV{IMPANNOT_ROOT}/default_config/default_hierarchy.txt";
my $pipeline_version_file = "$ENV{IMPANNOT_ROOT}/VERSION";
my $annotation_check_pl = "$ENV{IMPANNOT_ROOT}/bin/annotation_check.pl";

my $version = 0;
my $species = "";
my $genome_file ="";
my $genes_file = "";
my $tRNA_file = "";
my $rep_file = "";
my $hier_file = "";
my $miR_fasta_file = "";
my $miR_gff_file = "";
my $seq_files = "";
my $outDir = "";
my $a_release = 0;
my $add_to = 0;
my $skip_gtf_chr_check = 0;

my $help = 0;
my $debug = 0;

#sub required {
#   print <<EOH
#
#Required options:
#--version
#--species
#--genome
#--genes
#--miRNA_dat
#--miRNA_gff
#--out
#
#EOH
#}

sub help {
   print <<EOH;
This pipeline is intended to allow SequenceImp users to build customised
annotation sets to suit their own goals. The annotation sets generated will
be organised in such a way as to conform to the pipeline's input criteria.

The user can provide a minimal set of annotation files (genome, ensembl genes
gtf and the appropriate miRBase annotation files) or can expand this selection
to include tRNAs, Repeats and ncbi derived fasta sequences for detailed
analysis. This allows the user to customise the pipeline to their own purposes,
only limited by the genome annotation available in Ensembl and miRBase.

Usage:
$0 [options]
Options:

   --annot-version=<INTEGER>              Required: A user defined version number that 
                                    will be appended to uniquely identify the 
                                    annotation set within the user's Annotation
                                    directories. Different species can possess the
                                    same version number in order to group them.
                                    
   --genome-name=<STRING>           Required: A user defined species name to label
                                    the annotation.

   --genome-fasta=<STRING>          Required: The fasta genome sequence for the 
                                    species, downloaded from Ensembl 
                                    (http://www.ensembl.org) or Ensembl Genomes 
                                    (http://www.ensemblgenomes.org). This can be
                                    gzipped. This should not be one of
                                    the masked genomes that are available.

   --ensembl-genes=<STRING>         Required: The gene annotation gtf file
                                    downloaded from Ensembl (http://www.ensembl.org) 
                                    or Ensembl Genomes (http://www.ensemblgenomes.org). 

   --trnascan-tRNAs=<STRING>        Optional: tRNAscan-SE output for the genomic fasta
                                    file in tabular format. tRNAscan-SE can be found 
                                    here: http://selab.janelia.org/tRNAscan-SE/.
   
   --repeatmasker-repeats=<STRING>  Optional: repeat sequences for the genomic fasta 
                                    file identified by RepeatMasker in tabular format.
                                    RepeatMasker can be found here: 
                                    http://www.repeatmasker.org.

   --annot-hierarchy=<STRING>       Optional: the hierarchy for annotation classes 
                                    identified in the genes, tRNAs and repeats files.
                                    A default hierarchy is supplied but users should
                                    note this will not contain all annnotations classes
                                    for all species and ensembl/database releases.

   --mirbase-fasta=<STRING>          Required: the miRBase .fasta file corresponding to
                                    the gff3 file supplied below. The miRBase
                                    annotation database can be found here:
                                    http://www.mirbase.org.         

   --mirbase-gff=<STRING>            Required: the miRBase .gff3 file relevant to the 
                                    species and genome version in question. The miRBase
                                    annotation database can be found here:
                                    http://www.mirbase.org.
 
   --ncbi-seqs=<STRING>             Optional: Where these are supplied they will be 
                                    organised and formatted for use in the features:
                                    repeat step of the pipeline. Multiple fasta files
                                    can be supplied as a comma seperated string.

   --out=<STRING>                   Required: The directory into which the annotation
                                    will be arranged. Multiple 'species' and 'versions'
                                    can be stored within the same output directory
                                    allowing a user to build up an annotation repository
                                    for multiple projects.

   --add <FLAG>                     If specified this will trigger the pipeline to assume
                                    the output directory specified (--out) is already 
                                    an annotation created by this pipeline. All new 
                                    annotation will be added to the pre-exiting ENSEMBL, 
                                    MIRBASE, REPEATS and LOGS directories if available or 
                                    these will be created, while an existing VERSION
                                    file will be replaced. New files will be indexed by the 
                                    --genome and --version specified. If --add is not 
                                    specified an annotation directory will be placed in the 
                                    --out directory, indexed with the process ID.
                                    CAUTION: If flag is used, care should be taken here, 
                                    to ensure the pipeline will not overwrite pre-existing 
                                    files in the directory, where this is not expected. 
                                    This will occur automatically, so a new unique
                                    --version number should be given.

   --skip-ens-gtf-chr-check <FLAG>  In some instances the Ensembl gene GTF will contain
                                    gene annotation that may originate from chromosomes
                                    not seen in the Ensembl genome sequence file provided.
                                    These genes may correspond to genome patches etc. In 
                                    these instances the pipeline will remove this annotation
                                    while parsing the Ensembl gene GTF file. If this is
                                    expected, this flag will trigger the pipeline to skip
                                    default checks that ensure the gene GTF and sequence FASTA 
                                    chromosomes correspond.

   --help
EOH
}

if
(!GetOptions(
   "annot-version=i"        => \$version,
   "genome-name=s"          => \$species,
   "genome-fasta=s"         => \$genome_file,
   "ensembl-genes=s"        => \$genes_file,
   "trnascan-tRNAs=s"       => \$tRNA_file,
   "repeatmasker-repeats=s" => \$rep_file,
   "annot-hierarchy=s"      => \$hier_file,
   "mirbase-fasta=s"        => \$miR_fasta_file,
   "mirbase-gff=s"          => \$miR_gff_file,
   "ncbi-seqs=s"            => \$seq_files,
   "out=s"                  => \$outDir,
   "add"                    => \$add_to,
   "release"                => \$a_release,
   "skip-ens-gtf-chr-check" => \$skip_gtf_chr_check,
   "help"                   => \$help,
   "debug"                  => \$debug
   )
){
   print "Failed to parse the commmand line options\n";
   exit(1);
}

if ($help){
   help();
   exit(0);
}

# Ressolve between pipes and die statement
sub mysystem {
   my $cline = shift;
   my $retval = system("bash", "-c", qq{set -o pipefail; $cline});
   if ($retval) {
      print STDERR "Error while executing [$cline]";
   }
   return $retval;
}

# Check the pipeline version information
die "Pipeline version file does not seem to exist or is empty: $pipeline_version_file\nPlease check the pipeline installation\n" if (! -s $pipeline_version_file);
open(VERSION,"< $pipeline_version_file") || die "Could not open the pipeline version file\n";
my @version_lines = <VERSION>;
close(VERSION) || die "Could not close the pipeline version file\n";
die "Version file should contain a single line\n" if (scalar(@version_lines)!=1);

# Check options
die "Require an annotation version to be specified (see --annot-version)\n" if ($version == 0);
die "Require a species to be specified (see --genome-name)\n" if (length($species)==0);
die "Require a genome file to be specified (see --genome-fasta)\n" if (length($genome_file)==0);
die "Require a gene annotation file to be specified (see --ensembl-genes)\n" if(length($genes_file)==0);
die "Require a miRBase .fasta file to be specified (see --mirbase-fasta)\n" if (length($miR_fasta_file)==0);
die "Require a miRBase .gff file to be specified (see --mirbase-gff)\n" if (length($miR_gff_file)==0);
die "Require an output directory in which annotation is being stored (see --out)\n" if (length($outDir)==0);

die "Genome file does not appear to exist (see --genome-fasta)\n" if (! -s $genome_file);
die "Annotation file does not seem to exist (see --ensembl-genes)\n" if (! -s $genes_file);
die "miRBase .dat file does not seem to exist (see --mirbase-fasta)\n" if (! -s $miR_fasta_file);
die "miRBase .gff file doese not seem to exits (see --mirbase-gff)\n" if (! -s $miR_gff_file);
die "Output directory does not appear to exist (see --out)\n" if (! -d $outDir);

die "Specified version (--annot-version) must not be negative\n" if ($version < 0);

if (length($tRNA_file) != 0){
   die "tRNAscan-SE output file does not seem to exist (see --trnascan-tRNAs)\n" if (! -s $tRNA_file);
}

if (length($rep_file) != 0){
   die "Repeat file does not seem to exist (see --repeatmasker-repeats)\n" if (! -s $rep_file);
}

if (length($hier_file) != 0){
   die "Hierarchy file does not appear to exist (see --annot-hierarchy)\n" if (! -s $hier_file);
}else{
   $hier_file = $default_hierarchy_file;
   print STDERR "No hierarchy specified. Using the default hierarchy file:\t$default_hierarchy_file\n";
   die "Could not locate the default hierarchy file\n" if (! -s $hier_file);
}

my @indiv_ncbis;
if (length($seq_files) != 0){
   @indiv_ncbis = split ",", $seq_files, -1;
   for my $an_ncbi (@indiv_ncbis){
      die "Empty entry in ncbi file list (see --ncbi_seqs)\n" if length($an_ncbi)==0;
      die "$an_ncbi does not appear to exist (see --ncbi_seqs)\n" if (! -s $an_ncbi);
   }
}

# Create a specific output directory for annotation if --add not provided
if(! $add_to){
   $outDir = "$outDir/Annotation_$$";
   print STDERR "\nCreating a directory for the annotation: $outDir\n";
   die "Intended Annotation directory $outDir exists. Try again?\n" if (-d $outDir);
   die "Could not make Annotation directory $outDir\n" if (mysystem("mkdir -p $outDir"));
}else{
   print STDERR "\nAdding the annotation to the existing directory\n";
}

### Create a log file to record important information concerning the assembly of the annotation
sub log_file_rec{
   my $rec_file = shift;
   my $split_type = shift;
   if($split_type eq "no_split"){
      if ($rec_file =~ /([^\/]+)$/){
         return($1);
      }else{
         die "Could not find the file name to be recorded in the log for: $rec_file\n";
      }
   }else{
      my @all_file_split = split $split_type, $rec_file;
      my @file_return;
      for (@all_file_split){
         if (/([^\/]+)$/){
            push @file_return, $1;
         }else{
            die "Could not find the file name to be recorded in the log for: $rec_file\n";
         }        
      }
      my $collapsed_file_return =join ", ", @file_return;
      return ($collapsed_file_return);
   }
}

die "Could not make LOG directory in $outDir\n" if (mysystem("mkdir -p $outDir/LOG"));

my $logfile = "$outDir/LOG/$species"."_$version.log";
print STDERR "\nCreating log file to record annotation assembly parameters: $logfile\n";

die "Log file exists: $logfile\nThere is a risk that annotation may be overwritten.\nTo continue regardless first remove this log file.\n" if -e $logfile;

open (LOG, ">> $logfile") || die "Could not open the log file\n";
### List software versions and dependencies used and files specified for run
print LOG "Pipeline parameters\n\nVersion: ".$version_lines[0]."\n";

print LOG "\nOther software versions:\n\n";
print LOG "R version used for data compilation:\n\n";

# Check R version is compatible

open (RVERS, "R --version |") || die "Could not locate R\n";
my @R_compatible = (2,15.1);

my $found_version = 0;
my $first_no = 0;
my $second_no = 0;
my $third_no = 0;

while(<RVERS>){
   chomp;
   my $vers_line = $_;
   if($vers_line=~/^R version (\d+)\.(\d+.\d+)\s+/){
      die "Unable to identify R version being used - can't parse R header - multiple versions found\n" if ($found_version);
      $found_version = 1;
      $first_no = $1;
      $second_no = $2;
      if (($first_no > $R_compatible[0])||(($first_no == $R_compatible[0])&&($second_no >= $R_compatible[1]))){
         next;
      }else{
         die "R version found: [$vers_line]\nR version ".$R_compatible[0].".".$R_compatible[1]." or later required\n";
      }
   }
}
#print "[$first_no] [$second_no]";
close (RVERS) || die "Could not close the connection to R\n";
die "Unable to identify R version being used\n" if (! $found_version);
print LOG "R version: ".$first_no.".".$second_no."\n\n";

#die "Could not locate R version\n" if mysystem("R --version >> $logfile");

print LOG "Bowtie version used for data compilation:\n\n";
die "Could not locate Bowtie\n" if mysystem("bowtie --version >> $logfile");

print LOG "\n\nCommand line options:\n\n";
print LOG "Species annotated: $species\n";
print LOG "Annotation version created: $version\n";
print LOG "Genome file: ",log_file_rec($genome_file,"no_split"),"\n";
print LOG "Gene file: ",log_file_rec($genes_file,"no_split"),"\n";
print LOG "tRNA file: ",log_file_rec($tRNA_file,"no_split"),"\n" if (length($tRNA_file) != 0);
print LOG "Repeat file: ",log_file_rec($rep_file,"no_split"),"\n" if (length($rep_file) != 0);
print LOG "Hierarchy used: ",log_file_rec($hier_file,"no_split"),"\n";
print LOG "miRNA fasta file: ",log_file_rec($miR_fasta_file,"no_split"),"\n";
print LOG "miRNA gff: ",log_file_rec($miR_gff_file,"no_split"),"\n";
print LOG "NCBI sequences: ",log_file_rec($seq_files,","),"\n\n" if (length($seq_files) != 0);
#print LOG "";

print LOG "\nAnnotation hierarchy used by annotation pipeline:\n\n";

close (LOG) || die "Could not close the log file\n";

### Write hierarchy to the log file.

die "Could not write the annotation hierarchy to the log file: $logfile\n" if mysystem("cat $hier_file >> $logfile"); 

### Create a version file with last edit.
my $annotation_version_file = "$outDir/VERSION";
open (VERS,"> $annotation_version_file") || die "Could not open the version file: $annotation_version_file\n";
my @alltime = localtime();
if($a_release){
   print VERS "SEQIMP\t",(sprintf("%02d", $alltime[5] % 100)),"-",$alltime[7]; # Adapted from perldoc
}else{
   print VERS "Last modified\t",(sprintf("%02d", $alltime[5] % 100)),"-",$alltime[7];
}
close(VERS) || die "Could not close the version file\n";

### Begin by creating directory structure into which new annotation will be deposited

# 1) Make base Annotation directories if required.
print STDERR "\nCreating annotation directory tree:\t$species\tVersion: $version\n";
die "Could not make ENSEMBL directory in $outDir\n" if (mysystem("mkdir -p $outDir/ENSEMBL"));
die "Could not make MIRBASE directory in $outDir\n" if (mysystem("mkdir -p $outDir/MIRBASE"));
die "Could not make REPEATS directory in $outDir\n" if (mysystem("mkdir -p $outDir/REPEATS"));

#my $ensVersDir = "$outDir/ENSEMBL/Version_$version"; # 150714: Edited to remove redundant directory level.
#my $ensSpecDir = "$ensVersDir/$species"."_$version"; 
my $ensSpecDir = "$outDir/ENSEMBL/$species"."_$version";

my $ensBowDir = "$ensSpecDir/BowIndex";
my $ensGRDir = "$ensSpecDir/GRangesObjects";
my $genomeDir = "$ensSpecDir/Genome";

#my $miRVersDir = "$outDir/MIRBASE/Version_$version";
#my $miRSpecDir = "$miRVersDir/$species"."_$version"; # 150714: Edited to remove redundant directory level.
my $miRSpecDir = "$outDir/MIRBASE/$species"."_$version";

my $repDir = "$outDir/REPEATS/$species"."_$version";
my $repBowDir = "$repDir/BowIndex";

#### Ensure that same species and version not run before to prevent overwriting results

#die "\nERROR: Species annotation already performed for this version number. Please delete: $ensSpecDir\n" if (-e $ensSpecDir);
#die "\nERROR: Species annotation already performed for this version number. Please delete: $miRSpecDir\n" if (-e $miRSpecDir);

# 2) Make Ensembl directories.
print STDERR "Making Ensembl directories\n";
#die "Could not make the ENSEMBL version directory: $ensVersDir\n" if (mysystem("mkdir -p $ensVersDir"));
die "Could not make the ENSEMBL species directory: $ensSpecDir\n" if (mysystem("mkdir -p $ensSpecDir"));
die "Could not make the ENSEMBL Bowtie directory: $ensBowDir\n" if (mysystem("mkdir -p $ensBowDir"));
die "Could not make the ENSEMBL GRanges directory: $ensGRDir\n" if (mysystem("mkdir -p $ensGRDir"));
die "Could not make the ENSEMBL Genome directory: $genomeDir\n" if (mysystem("mkdir -p $genomeDir"));

# 3) Make MIRBASE directories.
print STDERR "Making miRNA directories\n";
#die "Could not make the MIRBASE version directory: $miRVersDir\n" if (mysystem("mkdir -p $miRVersDir"));
die "Could not make the MIRBASE species directory: $miRSpecDir\n" if (mysystem("mkdir -p $miRSpecDir"));

# 4) Make REPEATS directory.
print STDERR "Making repeat directories\n";
die "Could not make the REPEATS version directory: $repDir\n" if (mysystem("mkdir -p $repDir"));
die "Could not make the REPEATS Bowtie directory: $repBowDir\n" if (mysystem("mkdir -p $repBowDir"));

##### ORGANISING THE GENOMIC SEQUENCE FILES
print STDERR "\nPreparing the genome\n";
print STDERR "Parsing genome file to calculate chromosome sizes and classes\n";

my $genome_parsed = "$genomeDir/$species"."_$version.Reference.fa";
my $length_file = "$genomeDir/$species"."_$version.chrom_len.tab";
my $class_file = "$genomeDir/$species"."_$version.chrom_class.tab";

open (SEQOUT, "> $genome_parsed") || die "Could not open the output file\n";
open (SEQ, "gzip -cd $genome_file |") || die "Could not open input file\n";
open (LENOUT, "> $length_file") || die "Could not open the chromosome length file\n";
open (CLASSOUT ,"> $class_file") || die "Could not open the chromosome class file\n";

my %id_collection;
my $seqlength = 0;         
my $id = "";

print LENOUT "Chr\tLength\n";
print CLASSOUT "Chr\tClass\n";

while(<SEQ>){
   
   chomp;
   my $current = $_;
   if ($current =~ /^>(\S+).*/) {
      
      # Record the length of the last chromosome
      print LENOUT "$id\t$seqlength\n" if (length($id)!=0);
      $seqlength = 0;

      $id = $1;

      die "Duplicate sequence name ($id) in FASTA file: $genome_file\n" if exists($id_collection{$id});
      $id_collection{$id}++;      

      # Record the class of the next chromosome

      if($current=~/^>\S+\s+dna:(\S+)\s+.*/){
         print CLASSOUT "$id\t$1\n";
      }else{
         die "Could not discern the class of chromosome $id.\nPlease check the FASTA file header format matches the Ensembl convention (eg. >2L dna:chromosome ... )\nand ensure that the Ensembl genome used is not masked\n";
      }

      print SEQOUT "$current\n";
   }else{
      $seqlength += length($current);
      print SEQOUT "$current\n";
   }
}

print LENOUT "$id\t$seqlength\n" if (length($id)!=0);

close (SEQOUT) || die "Can't close the parsed sequence file\n";
close (SEQ) || die "Can't close the genomic fasta file\n";
close (LENOUT) || die "Can't close the chromosome length file\n";
close (CLASSOUT) || die "Can't close the chromosome class file\n";

### Check genomic features files for compatability, corresponding chromosome names etc.
#$genome_file
#$genes_file
#$tRNA_file
#$rep_file
#$hier_file
#$miR_fasta_file
#$miR_gff_file
#$seq_files
#$length_file
#$class_file

my $tRNAOpt = length($tRNA_file) != 0 ? "--tRNAs=$tRNA_file" : "";
my $repOpt = length($rep_file) != 0 ? "--reps=$rep_file" : "";
my $ncbiOpt = length($seq_files) != 0 ? "--ncbi=$seq_files" : "";
my $checkdebugOpt =  ($debug != 0) ? "--debug" : "";
my $skipCheckOpt = ($skip_gtf_chr_check!=0) ? "--no_gtf_chr_check" : "";

my $annotCheckCall = "$annotation_check_pl --fasta=$genome_file --chrlen=$length_file --chrclass=$class_file --miR_fasta=$miR_fasta_file --miR_gff=$miR_gff_file --genes=$genes_file $tRNAOpt $repOpt $ncbiOpt $skipCheckOpt $checkdebugOpt";
print STDERR "Annotation file check:\n$annotCheckCall\n\n" if $debug;
die "Annotation check has failed\n" if system($annotCheckCall);

### Organising miRNA annotation data

print STDERR "\nOrganising the miRNA Annotation data\n";

my $parsed_miRNA_table = "$miRSpecDir/$species"."_$version.miRNAs.tab";

#./miRBase_gff3_parse.pl --gff3=data/mouse/mmu.gff3 --fasta=mature.fa.gz > ig.txt

print STDERR "Parsing the miRBase gff3 and fasta files\n";
my $gff_parse_call = "$miRBase_gff_pl --gff3=$miR_gff_file --fasta=$miR_fasta_file > $parsed_miRNA_table";
print STDERR "\nParse GFF/FASTA File:\n$gff_parse_call\n" if $debug;
die "Could not parse the miRBase gff3 and fasta files\n" if mysystem($gff_parse_call);

print STDERR "\nConverting miRNA coordinates into GenomicRanges objects\n";
my $miRBaseGRs = "$miRSpecDir/$species"."_$version".".miRNA_coordinates.RData";
my $miRGRCall = "R --vanilla --slave --args --outFile=$miRBaseGRs --miR-table=$parsed_miRNA_table --chrLen=$length_file < $miRBase_GR_R";
print STDERR "\nConvert MicroRNA Annotation Call:\n$miRGRCall\n" if $debug;
die "Could not convert miRBase derived coordinates into a GenomicRanges object\n" if mysystem($miRGRCall);

### Organising repeat and gene annotation (from Ensembl, tRNAScan and RepeatMasker)

print STDERR "\nPreparing the genomic annotation\n";
my $formattedReps = "";
if(length($rep_file) > 0){
   print STDERR "Formatting standard RepeatMasker output for pipeline\n";
   $formattedReps = "$ensGRDir/$species"."_$version.reps_formatted.tab.gz";
   my $formatCall = "gzip -cdf $rep_file | $repFormat_pl | gzip -1 > $formattedReps";
   print STDERR "\nReformat RepeatMasker Output:\n$formatCall\n" if $debug;
   die "Could not reformat the RepeatMasker standard output file\n" if mysystem($formatCall);
}

print STDERR "Formatting Ensembl GTF annotation for the pipeline\n\n";

my $converted_gene_gtf_file = "$ensGRDir/$species"."_$version.ensembl_gene.txt.gz";
my $gtf_trans_call = "$gtfTranslate_pl --gtf=$genes_file --len=$length_file | gzip -c1 > $converted_gene_gtf_file";
print STDERR "\nReformat Ensembl Gene GTF:\n$gtf_trans_call\n" if $debug;
die "Could not parse the Ensembl GTF file for the pipeline\n" if mysystem("$gtf_trans_call");


my $repeat_file_option = "";
my $tRNA_file_option = "";

if (length($tRNA_file) != 0){
   $tRNA_file_option = "--tRNAs=$tRNA_file";
}

if (length($rep_file) != 0){
   die "Could not find the formatted repeat file: $formattedReps\n" if (! -s $formattedReps);
   $repeat_file_option = "--repeats=$formattedReps";
}

print STDERR "\nConverting the genomic annotation to GRanges objects\n";
my $annotation_convert_call = "R --vanilla --slave --args --outDir=$ensGRDir/ --genes=$converted_gene_gtf_file $tRNA_file_option $repeat_file_option --hierarchy=$hier_file --lengths=$length_file --species=$species --version=$version< $genomic_annot_conv_R";
print STDERR "\nConvert Annotation to GRanges Object Call:\n$annotation_convert_call\n" if $debug;
die "Could not convert the annotation supplied into GenomicRanges format\n" if mysystem($annotation_convert_call);

### Building a Bowtie index

print STDERR "Building a Bowtie index\n";

my $buildBowCall = "bowtie-build $genome_parsed $ensBowDir/$species"."_$version";
print STDERR "\nBuild Bowtie Index Call:\n$buildBowCall\n" if $debug;
die "Could not build the bowtie index: $buildBowCall\n" if (mysystem($buildBowCall));

### Collapsing the genome file
print STDERR "\nCompressing the genome file for storage\n";
my $genCompCall = "gzip -f9 $genome_parsed";
print STDERR "\nCompressing Genome Call:\n$genCompCall\n" if $debug;
die "Could not compress the genome file: $genCompCall\n" if (mysystem($genCompCall));

### Convert repeat sequences to Bowtie indexes
if (length($seq_files) != 0){
   print STDERR "\nPreparing Bowtie indecies for NCBI sequences\n";
   
   my %rep_index;
   for my $indiv_rep  (@indiv_ncbis){
      my $rep_basename = "";
     
      print STDERR "\nCreating Bowtie index corresponding to $indiv_rep\n";

      die "Empty entry in ncbi file list (see --ncbi_seqs)\n" if length($indiv_rep)==0;
      die "$indiv_rep does not appear to exist (see --ncbi_seqs)\n" if (! -s $indiv_rep);
      
      if($indiv_rep =~ /.*\.fasta[\.gz]*$/ || $indiv_rep =~ /.*\.fa[\.gz]*$/){
         if($indiv_rep =~ /.*?([\w+\.-]+)\.[fasta]+[\.gz]*/){
            $rep_basename = $1;
            $rep_basename =~ s/[^\w_]/_/g;
            die "Seen the same repeat twice\n" if exists($rep_index{$rep_basename});
            $rep_index{$rep_basename}++;
         }else{
            die "Unable to capture the base name of repeat fasta file: $indiv_rep\n";
         }
      }else{
         die "Repeat files are expected to have '.fa', '.fa.gz', '.fasta' or '.fasta.gz' suffixes\n";
      }
   
      my $repeat_setup_call = "$repeat_setup_sh -i $indiv_rep -o $repBowDir -r $rep_basename -s $species -l $repeat_length_pl -v $version";
      print STDERR "\nCall to create a repeat Bowtie index:\n$repeat_setup_call\n" if $debug;
      
      print "$rep_basename\t$repBowDir\n";
      die "Could not create a repeat Bowtie index for: $rep_basename\n" if mysystem($repeat_setup_call);    
   }

}

print STDERR "\nCompleted the assembly of the annotation for $species version $version\n";

