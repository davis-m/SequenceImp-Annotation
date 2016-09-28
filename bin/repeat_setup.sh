#!/usr/bin/env bash

###########################################################################

#This annotation pipeline is part of the Kraken framework which aims to
#facilitate RNA sequence analysis in a streamlined and efficient manner.
#Copyright (C) 2011 2012 2013 2014 EMBL - European Bioinformatics Institute 

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program (Please see the COPYING file for details).
#If not, see <http://www.gnu.org/licenses/>.

#Please send bug reports to:
#kraken@ebi.ac.uk

###########################################################################

### repeat_setup.sh:
# Provided with selected repeat canonical sequences derived from NCBI this script will create the Bowtie indexes.
# It also requires the seq_imp compatible species name (see seqimp_config), the repeat type with which you wish
# the repeat files to be labeled and the directory into which the repeat annotation will be built.
# Consensus sequences should be selected to match those used in the literature.

set -e

outDir=
inFile=
repeatname=
species=
#repeatLen="/homes/matdavis/svn/rnagen/seqimp_annot/scripts/Repeats/repeat_length.pl"

function clean_up {
   if [[ $? == 0 ]]; then
      echo "[happy] Script succeeded"
   else
      echo "[grumpy] Script failed"
   fi
}

trap clean_up SIGTERM EXIT

while getopts :i:o:r:s:l:v:h opt
do
   case "$opt" in
      o)
         outDir=$OPTARG
         ;;
      i)
         inFile=$OPTARG
         ;;
      r)
         repeatname=$OPTARG
         ;;
      s) 
         species=$OPTARG
         ;;
      l)
         repeatLen=$OPTARG
         ;;
      v)
         version=$OPTARG
         ;;
      h)
         cat <<EOH

### WELCOME TO THE REPEAT ANNOTATION PIPELINE ###

This pipeline will take the specified repeat FASTA file and convert it to a
Bowtie index.  In addition it will organise an accompanying repeat length file
required for analysis. The pipeline should be passed the seqimp annotation base
directory as it will build a repeat annotation directory in here that will be
compatible with the seqimp pipeline. The repeat pipeline will convert a single
repeat FASTA at a time and store the results in the relevant species
subdirectory. Note that the species name specified for each repeat should be
the version of the name that matches the pipeline nomenclature. For each
species any number of uniquely named repeat sequences can be added.  Samples
submitted to the analysis pipeline will be compared to each in turn.

-i The reapeat FASTA file from which the Bowtie index should be built
-o The annotation directory into which the repeat annotation directories should 
   be built
-s The name of the species to which the repeat belongs. The species name should 
   match seqimp nomenclature as it will be used to identify the repeat files
-r The repeat type for which the bowtie index is being built. This should be 
   informative and unique as it will be used by the seqimp to name result files.
   This name must contain only alphanumeric characters and underscores.
   The sequence should be traceable in the NCBI repository based upon this name.
-l The path to the repeat length script.
-v The version of the annotation being created
-h Displays these options

EOH
         exit
         ;;
      
      :)
         echo "Flag $OPTARG needs an argument"
         exit 1
         ;;
      ?)
         echo "Flag $OPTARG unknown"
         exit 1
         ;;
      esac
done


if [[ -z $outDir ]]; then
   echo "Require a -o output directory specification"
   exit 1
fi
if [[ -z $inFile ]]; then
   echo "Require a -i input file specification"
   exit 1
fi
if [[ -z $repeatname ]]; then
   echo "Require a -r repeat type specification"
   exit 1
fi
if [[ -z $species ]]; then
   echo "Require a -s species specification"
   exit 1
fi
if [[ -z $repeatLen ]]; then
   echo "Require a -l repeat length script path"
   exit 1
fi
if [[ -z $version ]]; then
   echo "Require a -v annotation version"
   exit 1
fi

bowtie-build $inFile $outDir/$species'_'$version.$repeatname
$repeatLen $inFile $outDir/$species'_'$version.$repeatname.length.txt



