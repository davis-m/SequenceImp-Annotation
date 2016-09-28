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

#### Converts API derived annotation into Granges class object.

library(GenomicRanges)

outFile <- ""
tRNAFile <- ""
geneFile <- ""
repeatFile <- ""
hierarchy <- ""
chrLenFile <- ""
species <- ""
tannot <- FALSE
rannot <- FALSE
version <- -1;

EnsemblAnnot <- list()

##############################
### Command line arguments ###
##############################

args <- R.utils::commandArgs(asValues=TRUE)


if(!is.null(args$out)){
   outDir <- args$out
}

if(!is.null(args$tRNAs)){
   tRNAFile <- args$tRNAs
}

if(!is.null(args$genes)){
   geneFile <- args$genes
}

if(!is.null(args$repeats)){
   repeatFile <- args$repeats
}

if(!is.null(args$hierarchy)){
   hierarchy <- args$hierarchy
}

if(!is.null(args$lengths)){
   chrLenFile <- args$lengths
}

if(!is.null(args$species)){
   species <- args$species
}

if(!is.null(args$version)){
   version <- as.integer(args$version)
}

if(version < 0){
   stop("Require a non-negative version number to be specified to define file names")
}

if(nchar(outDir)==0){
   stop("Require a file for the annotation GRanges output")
}

if(nchar(geneFile)==0){
   stop("Require a ensembl gene annotation file")
}

if(nchar(hierarchy)==0){
   stop("Require an annotation hierarchy file")
}

if(nchar(chrLenFile)==0){
   stop("Require an chromosome length file")
}

if(nchar(species)==0){
   stop("Require a specified species")
}

if(nchar(tRNAFile)!=0){
   tannot <- TRUE
}

if(nchar(repeatFile)!=0){
   rannot <- TRUE
}

###########################################
### Loading the Ensembl annotation data ###
###########################################
#stopifnot(FALSE)

write(paste("Reading gene file: ",geneFile),stderr())
geneData <- read.table(gzfile(geneFile), header = FALSE, sep="\t",quote="", check.names = FALSE, as.is = TRUE ,comment.char = "#")

write(paste("Reading chromosome length file: ",chrLenFile),stderr())
chrLenData <- read.table(gzfile(chrLenFile), header = TRUE, sep = "\t",quote="", check.names = FALSE, as.is = TRUE , comment.char = "")

# Only read the repeat and tRNA annotation if they have been specified and are present

if(rannot){
   write(paste("Reading repeat file: ",repeatFile),stderr())
#   repeatData <- read.table(gzfile(repeatFile), header = FALSE, sep = "\t",quote="", check.names = TRUE, as.is = TRUE , comment.char = "")
   repeatData <- read.table(gzfile(repeatFile), header = FALSE,skip=3, sep = "\t", fill=FALSE ,quote="", check.names = TRUE, as.is = TRUE , comment.char = "")
}else{
   write("Not handling additional RepeatMasker annotation",stderr())
   repeatData <- NULL
}

if (tannot){
   write(paste("Reading tRNA file: ",tRNAFile),stderr())
   tRNAData <- read.table(gzfile(tRNAFile), header = FALSE, sep = "",quote="", check.names = FALSE, as.is = TRUE , comment.char = "") # Note: seperate based on white space - tRNAscan-SE output delineated by tabs and spaces
}else{
   write("Not handling additional tRNAscan-SE annotation",stderr())
   tRNAData <- NULL
}


######################
### Load hierarchy ###
######################

write(paste("Reading in species hierarchy:",hierarchy),stderr())
hierarchyList <- scan(file=hierarchy,sep="\n", what="character")
if((! "unexpected_repeats" %in% hierarchyList) | (! "unexpected_genes" %in% hierarchyList)){
   stop("hierarchy must contain both 'unexpected_repeats' and 'unexpected_genes' classes.")
}

######################################
### Chromosome length reformatting ###
######################################

write("Reformatting chromosome lengths",stderr())
chrLens <- as.numeric(chrLenData[,"Length"])
names(chrLens) <- as.character(chrLenData[,"Chr"])

# For objects downstream remove any features that don't match the sequencess in the length file required for the GRanges objects.

chromCheck <- function(chrLengths, AnnotObject){
   if (any(! levels(factor(AnnotObject$Seqname)) %in% names(chrLengths))){
      excessChrs <- levels(factor(AnnotObject$Seqname))[which(! levels(factor(AnnotObject$Seqname)) %in% names(chrLengths))]
      write(paste("There are undefined chromosomes within the gene file:",paste(excessChrs, collapse=", "), sep=" "),stderr())
      write("WARNING: Features associated with these sequences are being removed from the annotation.",stderr())
      AnnotObject <- AnnotObject[! AnnotObject$Seqname %in% excessChrs,]
   }
   return(AnnotObject)
}
################################
## tRNA reformatting function ##
################################

### IS TRNASCAN-SE 1 BASED?
### tRNA introns?


tRNAFormat <- function(tRNADF){
   colnames(tRNADF) <- c("Seqname","Number","Start","End", "Class", "Anticodon", "Intron_start", "Intron_end", "Cove_score")
   tRNADF$Name <- paste("tRNAScan_",tRNADF$Seqname,"_",tRNADF$Number,sep="")
   tRNADFNeg <- tRNADF[tRNADF$Start > tRNADF$End,]
   tRNADFPos <- tRNADF[tRNADF$Start <= tRNADF$End,]
   
   tRNADFPos <- tRNAintronANDstrand(tRNADFPos,"+")
   tRNADFNeg <- tRNAintronANDstrand(tRNADFNeg,"-")
   return(rbind(tRNADFPos,tRNADFNeg))

}

tRNAintronANDstrand <- function(tRNAFrag,newStrand){
   tRNAFrag$Strand <- newStrand
   
   # Create classes for tRNA annotation - must differentiate from tRNA annotation already in Ensembl
   tRNAFrag[tRNAFrag$Class!="Pseudo","Class"] <- "tRNA"
   tRNAFrag$Class <- paste("tRNAScan_",tRNAFrag$Class,sep="")
   
   if (newStrand == "-"){
      temp <- tRNAFrag$Start
      tRNAFrag$Start <- tRNAFrag$End
      tRNAFrag$End <- temp
   }
   
   tRNAsAllInfo <- tRNAFrag[(tRNAFrag$Intron_start == 0 & tRNAFrag$Intron_end == 0) ,]   
   
   if(any(tRNAFrag$Intron_start != 0 | tRNAFrag$Intron_end != 0)){
      intronsPre <- tRNAFrag[(tRNAFrag$Intron_start != 0 | tRNAFrag$Intron_end != 0),]

      tRNAintronFirst <- intronsPre
      tRNAintronSecond <- intronsPre
      intronInfo <- intronsPre

      if (newStrand == "+"){
         tRNAintronFirst$End <- tRNAintronFirst$Intron_start -1
         tRNAintronSecond$Start <- tRNAintronFirst$Intron_end +1
         intronInfo$Start <- intronInfo$Intron_start
         intronInfo$End <- intronInfo$Intron_end
         intronInfo$Class <- paste(intronInfo$Class,"_Intron",sep="")
      }else if(newStrand == "-"){
         tRNAintronFirst$End <- tRNAintronFirst$Intron_end -1
         tRNAintronSecond$Start <- tRNAintronFirst$Intron_start +1
         intronInfo$Start <- intronInfo$Intron_end
         intronInfo$End <- intronInfo$Intron_start
         intronInfo$Class <- paste(intronInfo$Class,"_Intron",sep="")
      }else{
         stop("Unrecognised strand")
      }
      tRNAsAllInfo <- rbind(tRNAsAllInfo,tRNAintronFirst)
      tRNAsAllInfo <- rbind(tRNAsAllInfo,tRNAintronSecond)
      tRNAsAllInfo <- rbind(tRNAsAllInfo,intronInfo)
   }

   return(tRNAsAllInfo[,c("Class","Seqname","Start","End","Strand","Name","Anticodon")])
}

if(tannot){ 
   write("Formatting tRNA annotation",stderr())
   tRNAFormatted <- tRNAFormat(tRNAData)
   tRNAFormatted <- chromCheck(chrLens,tRNAFormatted)
}else{
   tRNAFormatted <- NULL
}
####### GOT TO HERE!!!

##################################
## Repeat reformatting function ##
##################################

repFormat <- function(repObject){
   repSlim <- repObject[,c("V11","V5","V6","V7","V9","V10")]   
   colnames(repSlim) <- c("Class","Seqname","Start","End","Strand","Name")
   if(any(repSlim$Strand != "+" & repSlim$Strand != "C" )){
      stop("Unrecognised strand in RepeatMasker output")
   }
   repSlim[(repSlim$Strand == "C"),"Strand"] <- "-"
   return(repSlim)
}


if (rannot){
   write("Formatting repeat annotation",stderr())
   repFormatted <- repFormat(repeatData)
   repFormatted <- chromCheck(chrLens,repFormatted)   
}else{
   repFormatted <- NULL
}


################################
## Gene reformatting function ##
################################

# 020714 - Old system depricated. Organising Ensembl GTF replaced with perl script for speed.
#genFormat <- function(geneTable){
#   colnames(geneTable) <- c("Seqname", "Class", "Fragment", "Start", "End", "Score", "Strand", "Frame", "Attribute")
#   geneTable$GeneID <- sub(".*gene_id\\s+\\\"(\\S+)\\\".*","\\1",geneTable$Attribute ,perl=TRUE) 
#   geneTable$TranscriptID <- sub(".*transcript_id\\s+\\\"(\\S+)\\\".*","\\1",geneTable$Attribute ,perl=TRUE)
#   geneTable$Exon_Number <- sub(".*exon_number\\s+\\\"(\\d+)\\\".*","\\1",geneTable$Attribute ,perl=TRUE)
#   
#   geneTable <- geneTable[geneTable$Fragment == "exon" | geneTable$Fragment == "CDS",c("Seqname", "Class", "Fragment", "Start", "End","Strand","GeneID","TranscriptID","Exon_Number")]   
#   return(geneTable)
#}

genFormat <- function(geneTable){
   colnames(geneTable) <- c("Seqname", "Start", "End", "Strand", "TranscriptID", "Class", "GeneID")
   return(geneTable)
}

#stopifnot(FALSE)

write("Formatting the gene annotation",stderr())
geneFormatted <- genFormat(geneData)
geneFormatted <- chromCheck(chrLens,geneFormatted)

## Clean up

rm(list=c("geneData","tRNAData","repeatData"))
gc()

###################################
## Conversion to GRanges Objects ##
###################################


###############################################
### Subroutines for forming Granges objects ###
###############################################

write("Beginning the converion of Ensembl features into GRanges objects",stderr())

makeGeneGR <- function(FormattedGs,chrLengths) {

   write("Creating GRanges objects from gene annotation",stderr())

   geneGR <-GRanges(
      seqnames = Rle(factor(FormattedGs[,"Seqname"], levels=names(chrLengths))),
      ranges   = IRanges(start=FormattedGs[,"Start"], end=FormattedGs[,"End"]),
      strand   = Rle(factor(FormattedGs[,"Strand"])),
      seqlengths = chrLengths,
      transcript = FormattedGs[,"TranscriptID"],
      classification = FormattedGs[,"Class"],
      gene = FormattedGs[,"GeneID"]
   ) 
  
   return(geneGR)
}


# 020714 - Altered - Formatting and break down of Ensembl GTF now handled by an additional Perl script as this is more efficient.
#makeGeneGR <- function(FormattedGs,chrLengths) {
#
#   transcriptList <- split(FormattedGs, FormattedGs$TranscriptID)
#
#   print("Creating GRanges objects from gene annotation, including break down into UTRs and CDSs")
#   print(paste("Transcript number:",length(transcriptList)))
#
#   geneGR <- lapply(transcriptList, function(x){
#   #lapply(ug, function(x){   
#      if (length(unique(x$Seqname))!=1){stop(paste("Transcript has multiple Seqnames:",x[1,"TranscriptID"]))}
#      if (length(unique(x$Strand))!=1){stop(paste("Transcript has multiple Strands:",x[1,"TranscriptID"]))}
#      if (length(unique(x$GeneID))!=1){stop(paste("Transcript has multiple Gene IDs:",x[1,"TranscriptID"]))}
#      if (length(unique(x$Class))!=1){stop(paste("Transcript has multiple Classes:",x[1,"TranscriptID"]))}
#      
#      GR <- GRanges(
#         seqnames = Rle(factor(x[,"Seqname"], levels=names(chrLengths))),
#         ranges   = IRanges(start=x[,"Start"], end=x[,"End"]),
#         strand   = Rle(factor(x[,"Strand"])),
#         seqlengths = chrLengths,
#         gene = x[,"GeneID"],      # Include gene name in gene EnsemblAnnot
#         transcript = x[,"TranscriptID"],
#         classification = x[,"Class"],
#         fragment = x[,"Fragment"]
#      ) 
#      
#      transcriptBreak <- split(GR,GR$fragment)
#      transcriptRange <- range(transcriptBreak$exon)
#      
#      if(!is.null(transcriptBreak$CDS)){
#         if(length(setdiff(transcriptBreak$CDS, transcriptBreak$exon))!=0){stop(paste("CDS extends beyond exons for:", x[1,"TranscriptID"]))} 
#
#         transcriptUTRs <- setdiff(transcriptBreak$exon, transcriptBreak$CDS) 
#
#         if(length(transcriptUTRs)!=0){
##            if(length(transcriptUTRs)>2){ ### Not necessarily at the end of transcript and can span multiple exons - Will now be distinguished by misc_UTR if mid CDS
##               stop(paste("More than 2 UTR sequences for a single transcript:",x[1,"TranscriptID"]))
##            }
#            
#            transcriptUTRs$gene <-  x[1,"GeneID"]
#            transcriptUTRs$transcript <- x[1,"TranscriptID"]
#            transcriptUTRs$classification <- x[1,"Class"]
#            transcriptUTRs$fragment <- "misc_UTR"
#            
#            if(x[1,"Strand"] == "+"){
#               if(any(start(transcriptUTRs) < min(start(transcriptBreak$CDS)))){   
#                  transcriptUTRs[start(transcriptUTRs) < min(start(transcriptBreak$CDS))]$fragment <- "5_UTR"
#               }
#               if(any(start(transcriptUTRs) > max(end(transcriptBreak$CDS)))){   
#                  transcriptUTRs[start(transcriptUTRs) > max(end(transcriptBreak$CDS))]$fragment <- "3_UTR"
#               }
#            }else if(x[1,"Strand"] == "-"){
#               if(any(start(transcriptUTRs) < min(start(transcriptBreak$CDS)))){   
#                  transcriptUTRs[start(transcriptUTRs) < min(start(transcriptBreak$CDS))]$fragment <- "3_UTR"
#               }
#               if(any(start(transcriptUTRs) > max(end(transcriptBreak$CDS)))){   
#                  transcriptUTRs[start(transcriptUTRs) > max(end(transcriptBreak$CDS))]$fragment <- "5_UTR"
#               }
#            }
#
#            transcriptInfo <- c(transcriptBreak$CDS,transcriptUTRs)
#         }else{
#            transcriptInfo <- transcriptBreak$CDS
#         }
#      }else{
#         transcriptInfo <- transcriptBreak$exon
#      }
#
#      transcriptIntrons <- setdiff(transcriptRange,transcriptBreak$exon)
#      
#      if(length(transcriptIntrons)>0){
#         transcriptIntrons$gene <-  x[1,"GeneID"]
#         transcriptIntrons$transcript <- x[1,"TranscriptID"]
#         transcriptIntrons$classification <- x[1,"Class"]
#         transcriptIntrons$fragment <- "Intron"
#
#         transcriptInfo <- c(transcriptInfo, transcriptIntrons)
#      }
#      return(transcriptInfo)
#   }) 
#   allGeneGR <- unlist(GRangesList(geneGR))
#   return(allGeneGR)
#}


makeRepeatGR <- function(FormattedRs, chrLengths) {
   repeatGR <- GRanges(
         seqnames = Rle(factor(FormattedRs[,"Seqname"], levels=names(chrLengths))),
         ranges   = IRanges(start=FormattedRs[,"Start"], end=FormattedRs[,"End"]),
         strand   = Rle(factor(FormattedRs[,"Strand"])),
         seqlengths = chrLengths,
         classification = FormattedRs[,"Class"],
         type = FormattedRs$Name
         )
   return(repeatGR)
}

maketRNAGR <- function(FormattedTs, chrLengths) {
   tRNAGR <- GRanges(
         seqnames = Rle(factor(FormattedTs[,"Seqname"], levels=names(chrLengths))),
         ranges = IRanges(start=FormattedTs[,"Start"], end=FormattedTs[,"End"]),
         strand = Rle(factor(FormattedTs[,"Strand"])),
         seqlengths = chrLengths,
         classification = FormattedTs[,"Class"],
         ID = FormattedTs[,"Name"],
         anticodon = FormattedTs[,"Anticodon"] 
        )
   return(tRNAGR)
}

################################################
### Creating GRanges objects for all classes ###
################################################

completeGeneGR <- makeGeneGR(geneFormatted,chrLens)

if(rannot){
   completeRepeatGR <- makeRepeatGR(repFormatted,chrLens)
}else{
   completeRepeatGR <- NULL
}

if(tannot){
   completetRNAGR <- maketRNAGR(tRNAFormatted,chrLens)
}else{
   completetRNAGR <- NULL
}

# Save as unreduced GR and then alter Gene classification to combine Classification and Fragment
unreducedGRList <- list("genes"=completeGeneGR,"repeats"=completeRepeatGR,"tRNAs"=completetRNAGR)
unreduced_outFile <- paste(outDir,"/",species,"_",version,".unreduced.GRanges.RData",sep="")
save(unreducedGRList, file = unreduced_outFile)

## Clean up

rm(list=c("unreducedGRList","tRNAFormatted","repFormatted","geneFormatted"))
gc()

##########################################################
### Streamline classification categories for reduction ###
##########################################################

#completeGeneGR$classification <- paste(completeGeneGR$classification,"_",completeGeneGR$fragment,sep="")
classes <- mcols(completeGeneGR)$classification
classes <- paste("Ensembl_",classes,sep="")
mcols(completeGeneGR) <- NULL
mcols(completeGeneGR)$classification <- classes

if(tannot){
   classes <- mcols(completetRNAGR)$classification
   mcols(completetRNAGR) <- NULL
   mcols(completetRNAGR)$classification <- classes
   completeGeneGR <- c(completeGeneGR,completetRNAGR)
}
if(rannot){
   classes <- mcols(completeRepeatGR)$classification
   classes <- paste("RepeatMasker_",classes,sep="")
   mcols(completeRepeatGR) <- NULL
   mcols(completeRepeatGR)$classification <- classes
}

##################################################
### Sorting biotype and repeat type categories ###
##################################################

# What if all gene types are expected?

geneTypes <- levels(factor(mcols(completeGeneGR)$classification))
repeatTypes <- NULL

if(rannot){
   repeatTypes <- levels(factor(mcols(completeRepeatGR)$classification))
   if (any(geneTypes %in% repeatTypes)){stop("Conflict between Gene and Repeat classifications")}

   repAnnotList <- split(completeRepeatGR,mcols(completeRepeatGR)$classification)
   
   if (any(!names(repAnnotList) %in% hierarchyList)){
      unexpRepNames <- which(!names(repAnnotList) %in% hierarchyList) 
     
      write(paste("RepeatMasker repeats not in hierarchy provided, so merged into 'unexpected_repeats':",paste(names(repAnnotList[unexpRepNames]), collapse=", ")),stderr())

      unexpRep <- unlist(repAnnotList[unexpRepNames])
      repAnnotList <- repAnnotList[-unexpRepNames]
      repAnnotList$unexpected_repeats <- unexpRep
      repeatTypes <- names(repAnnotList)
   }
}else{
   repeatTypes <- NULL
}

genAnnotList <- split(completeGeneGR,mcols(completeGeneGR)$classification)

if (any(!names(genAnnotList) %in% hierarchyList)){
   unexpGenNames <- which(!names(genAnnotList) %in% hierarchyList) 
      
   write(paste("Ensembl/tRNAScan annotations not in hierarchy provided, so merged into 'unexpected_genes':",paste(names(genAnnotList[unexpGenNames]), collapse=", ")),stderr())
   
   unexpGen <- unlist(genAnnotList[unexpGenNames])
   genAnnotList <- genAnnotList[-unexpGenNames]
   genAnnotList$unexpected_genes <- unexpGen
   geneTypes <- names(genAnnotList)
}

if(rannot){
   allAnnotList <- c(repAnnotList,genAnnotList)
}else{
   allAnnotList <- genAnnotList
}


rm(list=c("completeGeneGR","completeRepeatGR","completetRNAGR"))
gc()

#stopifnot(FALSE)

#################################################################################
### Save vector of RNA types and broad class to split Bowtie annotation plots ###
#################################################################################

# HERE

classRecords <- list()
classRecords[["GeneClasses"]] <- geneTypes
classRecords[["Repeats"]] <- repeatTypes

class_outFile <- paste(outDir, species,"_",version,".annotation.classes.RData", sep="")
save(classRecords, file = class_outFile)

#stopifnot(FALSE)
########################################
### Working with annotation overlaps ### #### Heirarchy should include all possible RNA types from all species.
########################################

write("Reducing GRanges objects to non-redundant annotation objects:",stderr())

allAnnotListSense <- GRangesList()
allAnnotListAntiSense <- GRangesList()
keptGRs <- GRanges()

# Remove annotation classes from heirarchy that are not required
hierarchyList <- hierarchyList[hierarchyList %in% names(allAnnotList)]

#########################################
### Functions for reducing annotation ###
#########################################

strandsANDreduce <- function(type,annotStuff){
   sense <- annotStuff[[type]]
   bothStrands <- ((!strand(sense) == '+')&(!strand(sense) == '-'))
   if(any(bothStrands)){strand(sense)[bothStrands] <- '+'}   # Assign double stranded features to '+' strand for sense strand
   senseStrandPos <- strand(sense)=='+'
   senseStrandNeg <- strand(sense)=='-'
   antisense <- sense
   if(any(senseStrandPos)){strand(antisense)[senseStrandPos] <- '-'}   # Swap strands for antisense
   if(any(senseStrandNeg)){strand(antisense)[senseStrandNeg] <- '+'}
   
   write(paste("Removing redundant sense annotation for:",type),stderr())

   tempsense <- removeRedundant(sense)
   if(length(tempsense)>0){
      allAnnotListSense[[type]] <<- tempsense
   }
   
   write(paste("Removing redundant antisense annotation for:",type),stderr())
   
   tempantisense <<- removeRedundant(antisense)  
   # Empty GRanges : invalid class “GRanges” object: 'seqlevels(seqinfo(x))' and 'levels(seqnames(x))' are not identical 
   if(length(tempantisense)>0){
      allAnnotListAntiSense[[type]] <<- tempantisense  
      # Empty GRanges : invalid class “GRanges” object: 'seqlevels(seqinfo(x))' and 'levels(seqnames(x))' are not identical
   }
#   write(paste("We have",nrow(elementMetadata(EnsemblAntiSense[[kind]])),"and",nrow(elementMetadata(EnsemblSense[[kind]])) ))
}

#annotationOverlaps <- GRangesList(lapply(antisenseGRs, intersect, reduce(unlist(GRs))))
removeRedundant <- function(gObject){
   if(!exists("keptGRs")){
      gObject <- reduce(gObject)      
      keptGRs <<- gObject      
   }else{
      gObject <- setdiff(gObject,keptGRs)   # Test case - returns what is in the first set but not the second
      keptGRs <<- union(keptGRs,gObject)    # Assigning to a global variable from a function - record current intervals accounted for
   }
   return(gObject)
} 

######################################
### Initiate the reduction process ###
######################################
#stopifnot(FALSE)
for(kind in hierarchyList){
   write(paste("Removing redundancy in",kind,"set"),stderr())
   strandsANDreduce(kind,allAnnotList)
}

EnsemblAnnotReduced <- list("sense" = allAnnotListSense,"antisense" = allAnnotListAntiSense)

write("Saving non-redundant annotation sets",stderr())

reduced_GRs_outFile <- paste(outDir,species,"_",version,".reduced.GRanges.RData",sep="")
save(EnsemblAnnotReduced, file = reduced_GRs_outFile)

write("Completed conversion of Ensembl annotation to GRanges objects",stderr())


#stopifnot(FALSE)

