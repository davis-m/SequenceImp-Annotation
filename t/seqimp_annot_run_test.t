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

###################################################################
# This script will run prerelease tests on the annotation pipeline
###################################################################

# TO DO:
# Annotation class RData file test

use strict;
use warnings;
use Test::More;
use Data::Dumper;
use Cwd 'abs_path';
use File::Spec;
use File::Basename;
use File::Path;
use IPC::Open3;
use Env;

my $SEQIMP_ANNOT = "imp_annot.pl";
my $SEQIMP_GTF_TRANSLATE = "gtf_translate.pl";
my $FIND_REF;

my ($test_script,$test_directory) = fileparse(abs_path($0));
my $SEQIMP_ANNOT_ROOT = abs_path(File::Spec->catfile($test_directory, ".." ));
my $RESOURCES_PATH = File::Spec->catfile($test_directory,"resources");
my $TEST_ANNOT_PATH = File::Spec->catfile($SEQIMP_ANNOT_ROOT,"../../data-for-tests");
my $TEST_DIR = File::Spec->catfile($SEQIMP_ANNOT_ROOT, "test-data");
my $TEST_RUN =  File::Spec->catfile($TEST_DIR,"run");
my $TEST_MISC =  File::Spec->catfile($TEST_DIR,"misc");
my $EXPECTED = File::Spec->catfile($RESOURCES_PATH,"expected-results.txt");
my $HIERARCHY = File::Spec->catfile($RESOURCES_PATH,"test-hierarchy.txt");
my $GRCHECK = abs_path(File::Spec->catfile($test_directory, "seqimp_annot_GR_test.R"));

# Add seqimp_annot bin to PATH
$ENV{PATH} .= ":".abs_path(File::Spec->catfile($SEQIMP_ANNOT_ROOT,"bin"));

# Cleanup test directories
sub cleanup {
   rmtree($TEST_DIR);
}

# Prepare test directories
sub prepare {
   mkpath($TEST_DIR);
   mkpath($TEST_RUN);
   mkpath($TEST_MISC);

   ok(-e $EXPECTED, "Test expected results");
   ok(-e $HIERARCHY, "Test hierarchy file");
}

# Load expected results file
sub load_expected{
   my %expected;
   open(EXP, "< $EXPECTED") || die "Could not open the file containing the expected test results\n";
   while(<EXP>){
      chomp;
      my ($key,$pattern) = split "\t";
      $expected{$key} = $pattern;
   }
   close(EXP) || die "Could not close the expected results file\n";
   return(\%expected);
}

# Run a test set
# Edited from Matloob's test scripts
sub runStep {
   my ($test_name, $command) = @_;

   my $pid = open3(\*IN, \*OUT, \*ERR, $command);
   close(IN) || die "STDIN feed would not close\n";
   my @outlines = <OUT>;
   my @errlines = <ERR>;
   waitpid($pid, 0);
   close(OUT) || die "STDOUT feed would not close\n";  
   close(ERR) || die "STDERR feed would not close\n";  
 
   is ($?, 0, "$test_name exit code");
   print "ERROR\n@errlines\n" if ($?);
}

# Test required programs are accessible
sub test_dependencies{
   runStep("R executable","R --version");
   runStep("perl executable","perl -v");
   runStep("bowtie executable","bowtie --version");
}

# File check - search for pattern
sub pattern_check {
   die "Need 3 arguments for pattern_check()\n" if  (scalar(@_)!=3);
   my $pattern = shift;
   my $file = shift;
   my $test_name = shift;

   open(TEST ,"gzip -cdf $file |") || die "Could not open $file\n";
   my @file = <TEST>;
   close(TEST) || die "Could not close $file\n";

   ok(grep(/$pattern/,@file), $test_name);
}

# File check - line number
sub line_check {
   die "Need 3 arguments for pattern_check()\n" if  (scalar(@_)!=3);
   my $line_number = shift;
   my $file = shift;
   my $test_name = shift;

   open(TEST ,"gzip -cdf $file |") || die "Could not open $file\n";
   my @file = <TEST>;
   close(TEST) || die "Could not close $file\n";

   ok(scalar(@file)==$line_number, $test_name);
}

# Test a miRBase GRanges object element
sub test_GR_object {
   my $R_call;
   
   die "Need 5 arguments for miR_GR_element()\n" if (scalar(@_)!=5);
   my $command_comment = $_[0];
   my @command_options = split ";", $_[1];
   my $annotFile = $_[2];
   my $type_of_test = $_[3];
   my $format = $_[4];
   if ($type_of_test eq "element"){
      if($format eq "miRBase"){
         die "miRBase GRanges element test requires 8 options\n" if (scalar (@command_options)!=8);
         $R_call = "R --vanilla --args --testtype=$type_of_test --RDataObject=$annotFile --GRtype=miRBase --element=$command_options[0] --chr=$command_options[7] --start=$command_options[5] --stop=$command_options[6] --strand=$command_options[4] --sequence=$command_options[3] --mature=$command_options[1] --precursor=$command_options[2] < $GRCHECK";
      }elsif($format eq "Ensembl"){
         die "Ensembl GRanges element test requires 6 options\n" if (scalar (@command_options)!=6);
         $R_call = "R --vanilla --args --testtype=$type_of_test --RDataObject=$annotFile --GRtype=Ensembl --class=$command_options[0] --orientation=$command_options[1] --chr=$command_options[2] --start=$command_options[3] --stop=$command_options[4] --strand=$command_options[5] < $GRCHECK";
      }else{
         die "Need a recognised GR test format\n";
      }
   }elsif($type_of_test eq "number"){
      if($format eq "miRBase"){
         die "miRBase GRanges length test requires 1 option\n" if (scalar (@command_options)!=1);
         $R_call = "R --vanilla --args --testtype=$type_of_test --RDataObject=$annotFile --GRtype=miRBase --number=$command_options[0] < $GRCHECK";
      }else{
         die "Need a recognised GR test format\n";
      }
   }elsif($type_of_test eq "chromosome"){
      if($format eq "miRBase"){
         die "miRBase GRanges chromosomal test requires 3 options\n" if (scalar (@command_options)!=3);
         $R_call = "R --vanilla --args --testtype=$type_of_test --RDataObject=$annotFile --GRtype=miRBase --chrcount=$command_options[0] --chosenchr=$command_options[1] --chrlength=$command_options[2]  < $GRCHECK";
      }elsif($format eq "Ensembl"){
         die "Ensembl GRanges chromosomal test requires 5 options\n" if (scalar (@command_options)!=5);
         $R_call = "R --vanilla --args --testtype=$type_of_test --RDataObject=$annotFile --GRtype=Ensembl --class=$command_options[0] --orientation=$command_options[1] --chrcount=$command_options[2] --chosenchr=$command_options[3] --chrlength=$command_options[4]  < $GRCHECK";
      }else{
         die "Need a recognised GR test format\n";
      }
   }else{
      die "Unrecognised GRanges test type\n";
   }
   runStep($command_comment,$R_call); 
}

# Test fragment dataset.
sub test_dm_frag{
   cleanup();
   prepare();
   
   my $annot_version = 1;
   my $genome_name = "fly";
   my $test_name = "Dm frags";

   my $resources_3R =  File::Spec->catfile($TEST_ANNOT_PATH,"drosophila_frags");
   my $genome_fasta =  File::Spec->catfile($resources_3R,"Drosophila_melanogaster.BDGP5.dna_frags_1millish.toplevel.fa.gz");
   my $ensembl_gtf =   File::Spec->catfile($resources_3R,"Drosophila_melanogaster.BDGP5.78_frags_millish.gtf.gz");
   my $mirbase_gff = File::Spec->catfile($resources_3R,"dme_frags_millish.gff3");
   my $mirbase_fasta =   File::Spec->catfile($resources_3R,"mature_frags_millish.fa.gz");
   my $repmasker_out = File::Spec->catfile($resources_3R,"dm3_frags_millish.fa.out.gz");
   my $trnascan_out = File::Spec->catfile($resources_3R,"Drosophila_melanogaster.BDGP5.78_frags_millish.out");
   my @ncbi_fastas = (File::Spec->catfile($resources_3R,"gypsy_element.fa"),File::Spec->catfile($resources_3R,"copia_element.fa"));
   my $old_format_gtf = File::Spec->catfile($resources_3R, "Drosophila_melanogaster_76_biotype_test.gtf.gz");

   ok(-e $genome_fasta, $test_name.": Test genomic FASTA");
   ok(-e $ensembl_gtf, $test_name.": Test Ensembl GTF");
   ok(-e $mirbase_fasta, $test_name.": Test miRBase FASTA");
   ok(-e $mirbase_gff, $test_name.": Test miRBase GFF");
   ok(-e $repmasker_out, $test_name.": Test RepeatMasker File");
   ok(-e $trnascan_out, $test_name.": Test tRNAScan-SE File");
   ok(-e $old_format_gtf, $test_name.": Test alternative Ensembl GTF");

   my $ncbicount = 0;
   foreach(@ncbi_fastas){
      $ncbicount++;
      ok(-e , $test_name.": Test NCBI FASTA ".$ncbicount);
   }

   my $ncbi_set = join "," , @ncbi_fastas;

   my $callFrags = "$SEQIMP_ANNOT --annot-version=$annot_version --genome-name=$genome_name --annot-hierarchy=$HIERARCHY --genome-fasta=$genome_fasta --ensembl-genes=$ensembl_gtf --mirbase-fasta=$mirbase_fasta --mirbase-gff=$mirbase_gff --out=$TEST_RUN --repeatmasker=$repmasker_out --trnascan-tRNAs=$trnascan_out --ncbi-seqs=$ncbi_set --add";
   #print("\n\n".$callFrags."\n\n");
   runStep("Drosophila fragment test", $callFrags);
   
   my $version_file =File::Spec->catfile($TEST_RUN,"VERSION");
   my $log_dir =  File::Spec->catfile($TEST_RUN,"LOG");
   my $log_file = File::Spec->catfile($log_dir,$genome_name."_".$annot_version.".log");

   my $mirbase_dir =File::Spec->catfile($TEST_RUN,"MIRBASE");
   my $mirbase_sample_dir =File::Spec->catfile($mirbase_dir,$genome_name."_".$annot_version);
   my $mirbase_coord_file =File::Spec->catfile($mirbase_sample_dir,$genome_name."_".$annot_version.".miRNA_coordinates.RData");
   my $mirbase_miRNA_file =File::Spec->catfile($mirbase_sample_dir,$genome_name."_".$annot_version.".miRNAs.tab");

   my $repeats_dir = File::Spec->catfile($TEST_RUN,"REPEATS");
   my $repeats_sample_dir =File::Spec->catfile($repeats_dir,$genome_name."_".$annot_version);
   my $repeats_bowindex_dir =File::Spec->catfile($repeats_sample_dir,"BowIndex");
   my $copia_index = File::Spec->catfile($repeats_bowindex_dir,$genome_name."_".$annot_version.".copia_element.rev.2.ebwt"); 
   my $gypsy_index = File::Spec->catfile($repeats_bowindex_dir,$genome_name."_".$annot_version.".gypsy_element.1.ebwt"); 
   my $gypsy_length = File::Spec->catfile($repeats_bowindex_dir,$genome_name."_".$annot_version.".gypsy_element.length.txt"); 
   
   my $ensembl_dir = File::Spec->catfile($TEST_RUN,"ENSEMBL"); 
   my $ensembl_sample_dir = File::Spec->catfile($ensembl_dir, $genome_name."_".$annot_version);
   my $bowtie_index_dir = File::Spec->catfile($ensembl_sample_dir,"BowIndex");
   my $bowtie_index = File::Spec->catfile($bowtie_index_dir, $genome_name."_".$annot_version.".4.ebwt");

   my $genome_dir = File::Spec->catfile($ensembl_sample_dir,"Genome");
   my $reference_file = File::Spec->catfile($genome_dir, $genome_name."_".$annot_version.".Reference.fa.gz");
   my $length_file = File::Spec->catfile($genome_dir, $genome_name."_".$annot_version.".chrom_len.tab");
   my $class_file = File::Spec->catfile($genome_dir, $genome_name."_".$annot_version.".chrom_class.tab");

   my $granges_dir = File::Spec->catfile($ensembl_sample_dir,"GRangesObjects");
   my $annotation_classes = File::Spec->catfile($granges_dir,$genome_name."_".$annot_version.".annotation.classes.RData");
   my $genes_formatted = File::Spec->catfile($granges_dir,$genome_name."_".$annot_version.".ensembl_gene.txt.gz");
   my $alternative_genes_formatted = File::Spec->catfile($TEST_MISC,"Pre_77_format.txt.gz");
   my $reduced_GRanges = File::Spec->catfile($granges_dir,$genome_name."_".$annot_version.".reduced.GRanges.RData");
   my $formatted_reps = File::Spec->catfile($granges_dir,$genome_name."_".$annot_version.".reps_formatted.tab.gz");
   my $unreduced_GRanges = File::Spec->catfile($granges_dir,$genome_name."_".$annot_version.".unreduced.GRanges.RData");


   ok(-s $version_file,$test_name.": Test version file");
   SKIP: {
      skip "Version file not found",1,if !(-s $version_file);
      pattern_check(${$FIND_REF}{Dm_frags_version_format},$version_file, $test_name.": Test version file format");
   }

   ok(-d $log_dir,$test_name.": Test LOG directory");
   ok(-s $log_file,$test_name.": Test LOG file");
   SKIP: {
      skip "LOG file not found",1,if !(-s $log_file);
      pattern_check(${$FIND_REF}{Dm_frags_log_hierarchy},$log_file, $test_name.": Test hierarchy logged");
   }

   ok(-d $mirbase_dir,$test_name.": Test MIRBASE directory");
   ok(-d $mirbase_sample_dir,$test_name.": Test miRBase sample directory");
   ok(-s $mirbase_miRNA_file,$test_name.": Test miRBase miRNA file");
   SKIP: {
      skip "miRBase table not found",2,if !(-s $mirbase_miRNA_file);
      pattern_check(${$FIND_REF}{Dm_frags_miRBase_table_header},$mirbase_miRNA_file, $test_name.": Test miRBase table header");
      pattern_check(${$FIND_REF}{Dm_frags_miRBase_table_MIMAT0000365},$mirbase_miRNA_file, $test_name.": Test miRBase table MIMAT0000365 check");
   } 
   ok(-s $mirbase_coord_file,$test_name.": Test miRBase coordinate GRanges file");
   SKIP:{
      skip "miRBase coordinate GRanges file not found",3,if !(-s $mirbase_coord_file);
      test_GR_object($test_name.": Test miRBase coordinate GRanges object for number of elements", ${$FIND_REF}{Dm_frags_miRBase_GRanges_number}, $mirbase_coord_file ,"number","miRBase");
      test_GR_object($test_name.": Test miRBase coordinate GRanges element - MIMAT0035232", ${$FIND_REF}{Dm_frags_miRBase_GRanges_MIMAT0035232}, $mirbase_coord_file ,"element","miRBase");
      test_GR_object($test_name.": Test miRBase coordinate GRanges chromosome - 3L", ${$FIND_REF}{Dm_frags_miRBase_GRanges_3L_length}, $mirbase_coord_file ,"chromosome","miRBase");
   }

   ok(-d $repeats_dir,$test_name.": Test REPEATS directory");
   ok(-d $repeats_sample_dir,$test_name.": Test repeat element sample directory");
   ok(-d $repeats_bowindex_dir,$test_name.": Test repeat element Bowtie index directory");
   ok(-s $copia_index, $test_name.": Test NCBI bowtie index 1");
   ok(-s $gypsy_index, $test_name.": Test NCBI bowtie index 2");
   ok(-s $gypsy_length, $test_name.": Test NCBI element length file");
   SKIP: {
      skip "NCBI element length file not found",1,if !(-s $gypsy_length);
      pattern_check(${$FIND_REF}{Dm_frags_gypsy_length},$gypsy_length, $test_name.": Test NCBI element length check");
   } 

   ok(-d $ensembl_dir,$test_name.": Test ENSEMBL directory");
   ok(-d $ensembl_sample_dir,$test_name.": Test Ensembl sample directory");
   ok(-d $bowtie_index_dir,$test_name.": Test Bowtie index directory");
   ok(-s $bowtie_index,$test_name.": Test genome Bowtie index");

   ok(-d $genome_dir,$test_name.": Test genome directory");
   ok(-s $reference_file,$test_name.": Test filtered reference file");
   SKIP: {
      skip "Reference sequence file not found", 1, if !(-s $reference_file);
      line_check(${$FIND_REF}{Dm_frags_genome_reference_lines},$reference_file, $test_name.": Test fragment reference file line number");
   }
   ok(-s $length_file,$test_name.": Test reference length file");
   SKIP: {
      skip "Reference length file not found",4,if !(-s $length_file);
      pattern_check(${$FIND_REF}{Dm_frags_genome_length_header},$length_file, $test_name.": Test fragment length header");
      pattern_check(${$FIND_REF}{Dm_frags_genome_length_3L},$length_file, $test_name.": Test 3L fragment length");
      pattern_check(${$FIND_REF}{Dm_frags_genome_length_3R},$length_file, $test_name.": Test 3R fragment length");
      line_check(${$FIND_REF}{Dm_frags_genome_length_lines},$length_file, $test_name.": Test fragment length file line number");
   }
   
   ok(-s $class_file,$test_name.": Test reference class file");
   SKIP: {
      skip "Reference class file not found",3,if !(-s $class_file);
      pattern_check(${$FIND_REF}{Dm_frags_genome_class_header},$class_file, $test_name.": Test fragment class header");
      pattern_check(${$FIND_REF}{Dm_frags_genome_class_3L},$class_file, $test_name.": Test 3L fragment class");
      line_check(${$FIND_REF}{Dm_frags_genome_class_lines},$class_file, $test_name.": Test fragment class file line number");
   }
   
   ok(-d $granges_dir,$test_name.": Test GRanges directory");
   ok(-s $annotation_classes,$test_name.": Test annotation class R object");
   ok(-s $genes_formatted,$test_name.": Test formatted gene file");
   SKIP: {
      skip "Reformatted gene file not found",9,if !(-s $genes_formatted);
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0037213_5UTR_I},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0037213 5'UTR I");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0037213_5UTR_II},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0037213 5'UTR II");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0037213_intron},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0037213 intron");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0037213_ORF},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0037213 ORF");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0037213_3UTR},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0037213 3'UTR");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0053294_exon},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0053294 exon");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0024277_3UTR},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0024277 3'UTR");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0024277_ORF_I},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0024277 ORF I");
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0024277_ORF_II},$genes_formatted, $test_name.": Test reformatted gene file - FBgn0024277 ORF II");
   }

   my $alternative = "$SEQIMP_GTF_TRANSLATE --gtf=$old_format_gtf --len=$length_file | gzip > $alternative_genes_formatted";
   #print("\n\n".$callFrags."\n\n");
   runStep("Reformatting alternative gene GTF file", $alternative);
   ok(-s $alternative_genes_formatted,$test_name.": Test alternative formatted Ensembl GTF file");
   SKIP: {
      skip "Reformatted alternative GTF not found",1,if !(-s $alternative_genes_formatted);
      pattern_check(${$FIND_REF}{Dm_frags_ensembl_FBgn0037213_ORF_Alternative},$alternative_genes_formatted,$test_name.": Test alternative reformatted gene file - FBgn0037213 ORF")
   }

   ok(-s $formatted_reps,$test_name.": Test formatted repeat file");
   SKIP: {
      skip "Reformatted repeat file not found",2,if !(-s $formatted_reps);
      pattern_check(${$FIND_REF}{Dm_frags_repeatmasker_formatted_element},$formatted_reps,$test_name.": Test reformatted repeat element");
      line_check(${$FIND_REF}{Dm_frags_repeatmasker_formatted_lines},$formatted_reps, $test_name.": Test reformatted repeat line number");
   }

   ok(-s $unreduced_GRanges,$test_name.": Test unreduced annotation GRanges file");
   ok(-s $reduced_GRanges,$test_name.": Test reduced annotation GRanges file");
   SKIP: {
      skip "Reduced annotation file not found",12,if !(-s $reduced_GRanges);
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - tRNA", ${$FIND_REF}{Dm_frags_reduced_GRanges_tRNA}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - antisense intron I", ${$FIND_REF}{Dm_frags_reduced_GRanges_intron_antisense_I}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - antisense intron II", ${$FIND_REF}{Dm_frags_reduced_GRanges_intron_antisense_II}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - sense intron I", ${$FIND_REF}{Dm_frags_reduced_GRanges_intron_sense_I}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - sense intron II", ${$FIND_REF}{Dm_frags_reduced_GRanges_intron_sense_II}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - sense ORF I", ${$FIND_REF}{Dm_frags_reduced_GRanges_ORF_sense_I}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - sense ORF II", ${$FIND_REF}{Dm_frags_reduced_GRanges_ORF_sense_II}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - antisense ORF", ${$FIND_REF}{Dm_frags_reduced_GRanges_ORF_antisense}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - antisense LTR", ${$FIND_REF}{Dm_frags_reduced_GRanges_LTR_antisense}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - sense 5'UTR", ${$FIND_REF}{Dm_frags_reduced_GRanges_5UTR_sense}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges element - antisense 5'UTR", ${$FIND_REF}{Dm_frags_reduced_GRanges_5UTR_antisense}, $reduced_GRanges ,"element", "Ensembl");
      test_GR_object($test_name.": Test Ensembl coordinate GRanges chromosome - 3R", ${$FIND_REF}{Dm_frags_reduced_GRanges_3R_length}, $reduced_GRanges , "chromosome", "Ensembl");
   } 
}

$FIND_REF = load_expected(); # Returns hash reference
test_dependencies();
test_dm_frag();
cleanup();
done_testing();

