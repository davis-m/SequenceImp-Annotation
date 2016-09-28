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

# This script will convert the miRNA coordinates derived from the miRBase .gff3
# file along with the fasta sequences derived from the miRBase .fasta file into
# a genomicRanges object for calculating miRNA depth.

# Report the number and ID of miRNAs that overlap

library(GenomicRanges)

miRBaseGRs <- list()

##############################
### Command line arguments ###
##############################

args <- R.utils::commandArgs(asValues=TRUE)

if (!is.null(args$out)){
   rObjectsFile <- args$out
}else{
   stop("Require an output file: --out=<FILE>")
}

if (!is.null(args$"miR-table")){
   miRFile <- args$"miR-table"
}else{
   stop("Require a path to the miRBase miRNA information table: --miR-table=<FILE>")
}

if (!is.null(args$chrLen)){
   chrLenFile <- args$chrLen
}else{
   stop("Require a chromosome length file: --chrLen=<FILE>")
}

###################################
### Upload the miRNA info table ###
###################################
write("Reading in the miRNA information file",stderr())

miRInfo <- read.delim(miRFile, header = TRUE, sep="\t", as.is = TRUE, quote = "" )

###################################
### Organise Chromosome Lengths ###
###################################
write("Finding chromosome length file",stderr())

chrLenData <- read.table(chrLenFile, header = TRUE, sep = "\t", check.names = FALSE, as.is = TRUE )
chrLens <- as.numeric(chrLenData[,"Length"])
names(chrLens) <- as.character(chrLenData[,"Chr"])

###############################
### Function for making GRs ###
###############################

### mature miRNA GR
### check all chromosome names match the chromosome length file

makeMatureGR <- function(matCds, chrStuff) { 
   if(! all(matCds$Chr %in% names(chrStuff))){stop("miRBase and Ensembl chromosomes do not correspond")}
   matGR <- GRanges(
      seqnames = Rle(factor(matCds[,"Chr"], levels=names(chrStuff))),
      ranges   = IRanges(start=matCds[,"Start"], end=matCds[,"End"]),
      strand   = Rle(factor(matCds[,"Strand"])),
      seqlengths = chrStuff,
      mature = matCds[,"Name"],
      id = matCds[,"ID"],
      precursor = matCds[,"Precursor"],
      sequence = matCds[,"Sequence"]
   )
   return(matGR)
}

############################
### Making precursor GRs ###
############################
write("Converting mature coordinates to GRanges format",stderr())

miRBaseGRs[["mature"]] <- makeMatureGR(miRInfo,chrLens)

#######################################
### Writing GRanges Objects to file ###
#######################################

save(miRBaseGRs, file = rObjectsFile)








