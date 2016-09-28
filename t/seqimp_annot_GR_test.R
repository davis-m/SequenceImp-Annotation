library("GenomicRanges")

# Provided by OSPDavis - many thanks

magicLoad<-function(file,pattern=""){
   env<-new.env()
   load(file,envir=env)
   length(ls(env,pattern=pattern)) > 0 || stop("cannot find an object matching ", pattern)
   return(eval(as.name(ls(env,pattern=pattern)[1]),envir=env))
}

args <- R.utils::commandArgs(asValues=TRUE)

if(!is.null(args$"testtype")){
  which_test <- as.character(args$"testtype")
}else{
   stop("Require a GenomicRanges object test type")
}

if(!is.null(args$RDataObject)){
   DataFile <- as.character(args$RDataObject)
}else{
   stop("Require an RData object")
}

if(!is.null(args$GRtype)){
   GRtype <- as.character(args$GRtype)
}else{
   stop("Need an specified annotation origin")
}

GR_annotation <- magicLoad(DataFile)

if(GRtype == "miRBase"){
   
   if(is.null(GR_annotation[["mature"]])){stop("Missing the 'mature' annotation")} 
   GR_selected <- GR_annotation[["mature"]] 

}else if(GRtype == "Ensembl"){
  
   if(!is.null(args$orientation)){
      orientation <- as.character(args$orientation)
   }else{
      stop("Need a specified orientation (sense/antisense)")
   }
   
   if(!is.null(args$class)){
      annotclass <- as.character(args$class)
   }else{
      stop("Need an specified annotation class")
   }
   
   if(is.null(GR_annotation[[orientation]])){stop("Missing the strand specified")} 
   if(is.null(GR_annotation[[orientation]][[annotclass]])){stop("Missing Ensembl annotation class")} 
   
   GR_selected <- GR_annotation[[orientation]][[annotclass]]
   
}else{
   stop("Unrecognised annotation type specified")
}

if (which_test == "chromosome"){
   if(!is.null(args$chrcount)){
      chrcount <- as.integer(args$chrcount)
   }else{
      stop("Require an expected number of chromosomes")
   }
   
   if(!is.null(args$chosenchr)){
      chosenchr <- as.character(args$chosenchr)
   }else{
      stop("Require a selected chromosome to test")
   }
   
   if(!is.null(args$chrlength)){
      chrlen <- as.integer(args$chrlength)
   }else{
      stop("Require an expected chromosome length")
   }

   chrStuff <- seqlengths(GR_selected)
   if(length(chrStuff) != chrcount){stop("Unexpected number of chromosomes")}
   if(is.na(chrStuff[chosenchr])){stop("Chosen chromosome information is not found")}
   if(chrStuff[chosenchr] != chrlen){stop("Chosen chromosome length is not as expected")}

}else if (which_test == "number"){
   
   if(!is.null(args$number)){
      how_many <- as.integer(args$number)
   }else{
      stop("Require an expected number of elements")
   }
   
   if (length(GR_selected) != how_many){stop("Unexpected number of elements in GRanges object")}

} else if (which_test == "element"){
   
   if(!is.null(args$chr)){
      chr <- as.character(args$chr)
   }else{
      stop("Require a miRBase chr")
   }
   
   if(!is.null(args$start)){
      start <- as.integer(args$start)
   }else{
      stop("Require a miRBase start")
   }
   
   if(!is.null(args$stop)){
      stop <- as.integer(args$stop)
   }else{
      stop("Require a miRBase stop")
   }
   
   if(!is.null(args$strand)){
      strand <- as.character(args$strand)
   }else{
      stop("Require a miRBase strand")
   }

   if(GRtype == "miRBase"){   
   
      if(!is.null(args$element)){
         element <- as.character(args$element)
      }else{
         stop("Require a miRBase ID")
      }
      
      if(!is.null(args$sequence)){
         sequence <- as.character(args$sequence)
      }else{
         stop("Require a miRBase sequence")
      }
      
      if(!is.null(args$mature)){
         mature <- as.character(args$mature)
      }else{
         stop("Require a miRBase mature")
      }
      
      if(!is.null(args$precursor)){
         precursor <- as.character(args$precursor)
      }else{
         stop("Require a miRBase precursor")
      }

      if(length(GR_selected[mcols(GR_selected)$id == element]) != 1){stop("Expecting one element with miRBase ID")}
      element_of_interest <- GR_selected[mcols(GR_selected)$id == element]
      if(as.character(strand(element_of_interest)) != strand){stop("strand does not match")}
      if(as.character(seqnames(element_of_interest)) != chr){stop("chr does not match")}
      if(as.integer(start(element_of_interest)) != start){stop("start does not match")}
      if(as.integer(end(element_of_interest)) != stop){stop("stop does not match")}
      if(as.character(element_of_interest$sequence) != sequence){stop("sequence does not match")}
      if(as.character(element_of_interest$mature) != mature){stop("mature does not match")}
      if(as.character(element_of_interest$precursor) != precursor){stop("precursor does not match")}
   }else if(GRtype == "Ensembl"){
      if(sum(paste(as.character(seqnames(GR_selected)),as.integer(start(GR_selected)),as.integer(end(GR_selected)),as.character(strand(GR_selected)),sep=";")==paste(chr,start,stop,strand,sep=";"))!=1){stop("Expected coodinate info does not match expected element")}
   }else{
      stop("Confused about annotation type")
   }

}else{
   stop("Unrecognised element to test")
}



