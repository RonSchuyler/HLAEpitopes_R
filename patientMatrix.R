# patientMatrix.R
# The patientMatrix class holds data read directly from a file, listing patients, their alleles, and diagnosis.
# 

source("hla3.R");

# Class patientMatrix
setClass("patientMatrix", 
   representation(
      dataMat="data.frame", 
      dataFile="character", 
      control="logical", # control or affected
      allOrNone="logical", # If true, patient is ignored if any null alleles.
      patientCount="numeric"), # How many patients are counted.
   prototype(dataFile="../data/mhc.csv", control=FALSE, dataMat=data.frame())
);
# patientMatrix constructor
patientMatrix <- function(dataFile="../data/mhc.csv", control=FALSE){
   #print("patientMatrix initialize");
   this <- new("patientMatrix", dataFile=dataFile, control=control);
   this <- loadData(this);
   this@allOrNone <- FALSE;
   this@patientCount <- 0;
   this;
}

# patientMatrix loadData
if (!isGeneric("loadData")) {
     if (is.function("loadData")){ fun <- loadData; }
     else{ fun <- function(object) standardGeneric("loadData"); }
     setGeneric("loadData", fun);
}
setMethod("loadData", "patientMatrix", function(object){
   #print(sprintf("loadData: my dataFile:%s", object@dataFile));
   # Set up variables for this file format.
   if(any(grep("COPD", object@dataFile, fixed=TRUE))){
      # COPD datafile
      object <- loadDataCOPD(object);
   }
   else{
      # Assume datafile is in mhc.csv format.
      object <- loadDataMhc.csv(object);
   }
   object;
})

# patientMatrix loadDataMhc.csv
if (!isGeneric("loadDataMhc.csv")) {
     if (is.function("loadDataMhc.csv")){ fun <- loadDataMhc.csv; }
     else{ fun <- function(object) standardGeneric("loadDataMhc.csv"); }
     setGeneric("loadDataMhc.csv", fun);
}
setMethod("loadDataMhc.csv", "patientMatrix", function(object){
   #print(sprintf("loadDataMhc.csv w dataFile %s", object@dataFile));
   # dx = diagnosis: "Affected" | "Control"
   # "Control" = "Matched Control" + "Random Control"
   mhc_table <- read.table(object@dataFile, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);
   #dx_col_name <- "Dx";
   affectedDx <- "Affected";
   controlDx <- "Control"; 
   #controlDx <- "Matched Control"; 
   # Affected or control subjects?
   if(object@control){
      dx <- controlDx;
   }
   else{
      dx <- affectedDx;
   }
   # Indices of rows with this Dx.
   #rows <- grep(dx, mhc_table$Dx, fixed=TRUE);
   rows <- grep(dx, mhc_table$Dx, ignore.case=TRUE);
   object@dataMat <- mhc_table[rows,];
   object;
})

# patientMatrix loadDataCOPD
if (!isGeneric("loadDataCOPD")) {
     if (is.function("loadDataCOPD")){ fun <- loadDataCOPD; }
     else{ fun <- function(object) standardGeneric("loadDataCOPD"); }
     setGeneric("loadDataCOPD", fun);
}
setMethod("loadDataCOPD", "patientMatrix", function(object){
   #print(sprintf("loadDataCOPD w dataFile %s", object@dataFile));
   mhc_table <- read.table(object@dataFile, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);
   # Replace "DNA" with "HLA" in column names to be consistent with mhc.csv.
   colnames(mhc_table) <- gsub("DNA", "HLA", colnames(mhc_table));
   #dx_col_name <- "GOLD";
   affectedDx <- "4";
   controlDx <- "0"; 
   if(object@control){
      dx <- controlDx;
   }
   else{
      dx <- affectedDx;
   }
   rows <- grep(dx, mhc_table$GOLD, fixed=TRUE);
   object@dataMat <- mhc_table[rows,];
   # Check for '-' in second allele, which means same as first allele.
   for(row_index in 1:nrow(object@dataMat)){
      for(loci_index in c(6,8)){
         allele <- as.character(object@dataMat[row_index,loci_index]);
         if(allele == '-'){
            # Replace with first allele.
            object@dataMat[row_index, loci_index] <- 
                               object@dataMat[row_index, (loci_index-1)];
         }
      }
   }
   object;
})


# patientMatrix alleleCountHash
# Fill a list with allele counts at locus
# Its main use is to generate a list of all the allele names in a dataset.
# That list can then be used to get all sequence data necessary for an analysis.
if (!isGeneric("alleleCountHash")) {
     if (is.function("alleleCountHash")){ fun <- alleleCountHash; }
     else{ fun <- function(this, locus) standardGeneric("alleleCountHash"); }
     setGeneric("alleleCountHash", fun);
}
# Should be ok to use cross-loci:
# locus=
#     "DRB1" OK
#     c("HLA.DQB1.1", "HLA.DQB1.2") OK
#     c("HLA.DQB1.1", "HLA.DQB1.2", "HLA.DRB1.1", "HLA.DRB1.2") OK 
#     c("DRB1", "DQB1") NOT ok
setMethod("alleleCountHash", "patientMatrix", function(this, locus){
   #print(sprintf("patientMatrix alleleCountHash for locus:%s", paste(locus, collapse=' ')));
   # locus may be length 1 
   #  ex: "DRB1" will be converted to c("HLA.DRB1.1", "HLA.DRB1.2")
   # OR length > 1 ex: c("HLA.DRB1.1", "HLA.DRB1.2", ...)
   if(is.null(locus) || length(locus) == 0){
      stop("patientMatrix alleleCountHash no locus");
   }
   else if(length(locus) == 1){
      # expand it
      locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
      #print(sprintf("expanded locus to:%s",  paste(locus, collapse=' ')));
   }
   # locus should now contain the names of the rows that we want

   # Find which columns have these loci names.
   columns <- which(locus == colnames(this@dataMat));

   # Count the times each allele name appears in this@dataMat.
   countHash <- list();
   for(row_index in 1:nrow(this@dataMat)){
      for(loci_index in columns){
         allele <- as.character(this@dataMat[row_index,loci_index]); # DRB1*0301
         ## 1/13/2009
         ## Checking length is not necessary.  The names are checked against 
         ## the entries in AlleleImport.txt.
         # Check the number of characters after the '*'.
         # Fewer than 4 is ambiguous. 
         #splitAllele <- unlist(strsplit(allele, '*', fixed=TRUE));
         #if(length(splitAllele) < 2 || nchar(splitAllele[2]) < 4){
         #if(nchar(allele, type="chars") < 9){
            # not a complete allele name, maybe DRB1*02 or A*01
            # short names are not specific enough to get sequence data
            #next;
         #}
         ## Just check that this is not null.
         if(is.null(allele) || is.na(allele) || (length(allele) < 1) ||
             nchar(allele) == 0){
            next;
         }
         if(allele %in% names(countHash)){
            countHash[[allele]] <- countHash[[allele]] + 1;
         }
         else{
            countHash[[allele]] <- 1;
         }
         # A better way may be to read contents of AlleleImport file and
         # check directly if each allele in the data file is represented there.
         # if any for this patient (row) are too short skip all of them?
      }
   }
   return(countHash);
})

# patientMatrix patientsWithAllele
# Return two hashes of counts keyed on allele names.
# If oneORtwo is TRUE, first hash is count of patients with 1 or 2 copies.
# Else, first is a count of patients with one copy of the named allele.
# Second is a count of patients homozygous for the named allele.
# Note: this includes null and ambiguous alleles.  
# Not really what is intended; think hard before using this method.
if (!isGeneric("patientsWithAllele")) {
   if (is.function("patientsWithAllele")){ fun <- patientsWithAllele; }
   else{ fun <- function(this, locus, oneORtwo=TRUE) standardGeneric("patientsWithAllele"); }
   setGeneric("patientsWithAllele", fun);
}
# Don't use this cross-loci:
# locus=
#     "DRB1" OK
#     c("HLA.DQB1.1", "HLA.DQB1.2") OK
#     c("DRB1", "DQB1") NOT ok
#     c("HLA.DQB1.1", "HLA.DQB1.2", "HLA.DRB1.1", "HLA.DRB1.2") NOT ok
setMethod("patientsWithAllele", "patientMatrix", function(this, locus, oneORtwo=TRUE){
   if(is.null(locus) || length(locus) == 0){
      stop("patientMatrix patientsWithAllele no locus");
   }
   else if(length(locus) == 1){
      # expand it
      locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
   }
   # locus should now contain the names of the rows that we want
   # Find which columns have these loci names. Should be length=2.
   columns <- which(locus == colnames(this@dataMat));

   # Count the times each allele name appears in this@dataMat.
   heterozygousHash <- list();
   homozygousHash <- list();
   for(row_index in 1:nrow(this@dataMat)){
      # Assume allele1 is columns[1], allele2 is columns[2]
      allele1 <- as.character(this@dataMat[row_index, columns[1]]);
      allele2 <- as.character(this@dataMat[row_index, columns[2]]);
      if(allele1 == allele2){
         # Homozygous.
         if(allele1 %in% names(homozygousHash)){
            homozygousHash[[allele1]] <- homozygousHash[[allele1]] + 1;
         }
         else{
            homozygousHash[[allele1]] <- 1;
         }
         if(oneORtwo){
            if(allele1 %in% names(heterozygousHash)){
               heterozygousHash[[allele1]] <- heterozygousHash[[allele1]] + 1;
            }
            else{
               heterozygousHash[[allele1]] <- 1;
            }
         }
      }
      else{
         # Heterozygous, count both.
         for(allele in c(allele1, allele2)){
            if(allele %in% names(heterozygousHash)){
               heterozygousHash[[allele]] <- heterozygousHash[[allele]] + 1;
            }
            else{
               heterozygousHash[[allele]] <- 1;
            }
         }
      }
   }
   return(list(heterozygousHash, homozygousHash));
})

# Return a 2-column matrix from two lists.
#two_lists_as_mat <- function(list1, list2)

# Similar to what as.matrix should do.
#get_hash_values_as_matrix <- function(hash)


# patientMatrix moduleCountHash
# Fill a matrix with counts of patients that have each module at the 
# given position vector.
# If all
if (!isGeneric("moduleCountHash")) {
     if (is.function("moduleCountHash")){ fun <- moduleCountHash; }
     else{ fun <- function(this, locus, posV, seqMat) standardGeneric("moduleCountHash"); }
     setGeneric("moduleCountHash", fun);
}
# Don't use cross-loci:
# locus=
#     "DRB1" OK
#     c("HLA.DQB1.1", "HLA.DQB1.2") OK
#     c("HLA.DQB1.1", "HLA.DQB1.2", "HLA.DRB1.1", "HLA.DRB1.2") NOT OK 
#     c("DRB1", "DQB1") NOT ok
setMethod("moduleCountHash", "patientMatrix", 
   function(this, locus, posV, seqMat){
   #print(sprintf("patientMatrix moduleCountHash for locus:%s", paste(locus, collapse=' ')));
   # locus may be length 1 
   #  ex: "DRB1" will be converted to c("HLA.DRB1.1", "HLA.DRB1.2")
   # OR length > 1 ex: c("HLA.DRB1.1", "HLA.DRB1.2", ...)
   if(is.null(locus) || length(locus) == 0){
      stop("patientMatrix moduleCountHash no locus");
   }
   else if(length(locus) == 1){
      # expand it
      locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
      #print(sprintf("expanded locus to:%s",  paste(locus, collapse=' ')));
   }
   # locus should now contain the names of the rows that we want

   # Find which columns have these loci names.
   columns <- which(locus == colnames(this@dataMat));

   # Count the times each allele name appears in this@dataMat.
   countHash <- list();
   # Examine each row.  Get the allele names, then the modules, then count them.
   for(row_index in 1:nrow(this@dataMat)){
      allele1 <- as.character(this@dataMat[row_index,columns[1]]); # DRB1*0301
      module1 <- c();
      # Check if this allele is high resolution typed.
      # Short allele names are not specific enough to get sequence data.
      #if(!is.null(allele1) && !is.na(allele1) && nchar(allele1) >= 9 && 
      # Being in the names of seqMat is a good enough check.
      if(!is.null(allele1) && !is.na(allele1) && 
                                                allele1 %in% rownames(seqMat)){
         # Allele is ok, now check that there are no padding characters (zeros)
         # in the module positions.
         if(!any(seqMat[allele1, posV] == 0)){
            module1 <- paste(seqMat[allele1, posV], collapse='');
         }
         else{
            #module1 <- paste(seqMat[allele1, posV], collapse='');
            #print(sprintf("skipping %s", module1));
            module1<-c();
         }
      }

      allele2 <- as.character(this@dataMat[row_index,columns[2]]); # DRB1*0301
      module2 <- c();
      if(!is.null(allele2) && !is.na(allele2) && 
                                                allele2 %in% rownames(seqMat)){
         if(!any(seqMat[allele2, posV] == 0)){
            module2 <- paste(seqMat[allele2, posV], collapse='');
         }
         else{
            #module2 <- paste(seqMat[allele2, posV], collapse='');
            #print(sprintf("skipping %s", module2));
            module2<-c();
         }
            
      }

      # If allOrNone, any null module means skip both.
      if(this@allOrNone){
         if(is.null(module1) || is.null(module2)){
            module1 <- c();
            module2 <- c();
         }
      }

      # Count them.
      if(!is.null(module1)){
         if(module1 %in% names(countHash)){
            countHash[[module1]] <- countHash[[module1]] + 1;
         }
         else{
            countHash[[module1]] <- 1;
         }
      }
      if(!is.null(module2)){
         # One allele counts as half of a patient.  
         #this@patientCount <- this@patientCount + 0.5;
         if(is.null(module1) || module1 != module2){
            # Because we're counting patients, 
            # only add module2 if it is not the same as module1.
            if(module2 %in% names(countHash)){
               countHash[[module2]] <- countHash[[module2]] + 1;
            }
            else{
               countHash[[module2]] <- 1;
            }
         }
      }
   }

   return(countHash);

})


# patientMatrix countPatients
# Count patients and fill this@patientCount.  
# Verify each allele exists in the sequence matrix (which is now built from Alignments.zip) 
# this@allOrNone controls how patients are counted.
# If allOrNone=TRUE, if a patient has any untyped, missing or ambiguous 
# alleles, all alleles for that patient are ignored.  
# If allOrNone=FALSE, each allele counts as half of a patient if
# countHalves=TRUE.  If countHalves=FALSE the patient is counted if either
# allele is not null.
# Note: this is not a replacement method.  It is up to the caller to insert
# the returned value into the calling object. 
if (!isGeneric("countPatients")) {
   if (is.function("countPatients")){ fun <- countPatients; }
   else{fun <- function(this, locus, seqMat, logfile, countHalves=FALSE) standardGeneric("countPatients");}
   setGeneric("countPatients", fun);
}
# Don't use this cross-loci:
# locus=
#     "DRB1" OK
#     c("HLA.DQB1.1", "HLA.DQB1.2") OK
#     c("DRB1", "DQB1") NOT ok
#     c("HLA.DQB1.1", "HLA.DQB1.2", "HLA.DRB1.1", "HLA.DRB1.2") NOT ok
setMethod("countPatients", "patientMatrix", function(this, locus, seqMat, logfile, countHalves=FALSE){
   if(is.null(locus) || length(locus) == 0){
      stop("patientMatrix countPatients no locus");
   }
   else if(length(locus) == 1){
      # expand it
      locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
   }
   # locus should now contain the names of the rows that we want
   # Find which columns have these loci names. Should have length=2.
   columns <- which(locus == colnames(this@dataMat));

   # Examine each row.  Get the allele names, then verify in file or hash.
   for(row_index in 1:nrow(this@dataMat)){
      messageToLog <- c();
      allele1 <- as.character(this@dataMat[row_index,columns[1]]); # DRB1*0301
      # Being in the names of seqMat is a good enough check.
      if(is.null(allele1) || is.na(allele1) || nchar(allele1)==0){
         allele1 <- c(); # not interesting, don't bother logging
      }
      else{
         # This allele is not null, now check if we have sequence data for it.
         if(!allele1 %in% rownames(seqMat)){
            # It's not empty, but it's not in the sequence matrix.  
            # Report the fact that we have to skip this allele.
            messageToLog <- sprintf("allele1 %s not found", allele1);
            allele1 <- c(); 
         }
      }

      allele2 <- as.character(this@dataMat[row_index,columns[2]]); # DRB1*0301
      if(is.null(allele2) || is.na(allele2) || nchar(allele2)==0){
         allele2 <- c();
      }
      else{
         if(!allele2 %in% rownames(seqMat)){
            # This incident will be reported.
            if(is.null(messageToLog)){
               messageToLog <- sprintf("allele2 %s not found", allele2);
            }
            else{
               messageToLog <- sprintf("%s; allele2 %s not found", 
                                          messageToLog, allele2);
            }
            allele2 <- c();
         }
      }

      # If allOrNone, any null module means skip both.
      if(this@allOrNone){
         if(!is.null(allele1) && is.null(allele2)){
            if(is.null(messageToLog)){
               messageToLog <- sprintf("allele1 %s skipped because allele2 is null and allOrNone=TRUE", allele1, messageToLog);
            }
            else{
               messageToLog <- sprintf("allele1 %s skipped because %s and allOrNone=TRUE", allele1, messageToLog);
            }
            allele1 <- c();
         }
         else if(!is.null(allele2) && is.null(allele1)){
            if(is.null(messageToLog)){
               messageToLog <- sprintf("allele2 %s skipped because allele1 is null and allOrNone=TRUE", allele1, messageToLog);
            }
            else{
               messageToLog <- sprintf("allele2 %s skipped because %s and allOrNone=TRUE", allele1, messageToLog);
            }
            allele2 <- c();
         }
      }

      # Check for a message to log.
      if(!is.null(messageToLog)){
         logToFile(logfile, logmessage=messageToLog, timestamp=FALSE);
      }

      # Count them.
      if(countHalves == TRUE){
         if(!is.null(allele1)){
            # One allele counts as half of a patient.  
            this@patientCount <- this@patientCount + 0.5;
         }
         if(!is.null(allele2)){
            # One allele counts as half of a patient.  
            this@patientCount <- this@patientCount + 0.5;
         }
      }
      else{
         if(!is.null(allele1) || !is.null(allele2)){
            this@patientCount <- this@patientCount + 1;
         }
      }

   }
   return(this@patientCount);
})

# trimToLocus
# Drop all unnecessary columns from dataMat,
# leaving only columns for the given locus.
if (!isGeneric("trimToLocus")) {
     if (is.function("trimToLocus")){ fun <- trimToLocus; }
     else{ fun <- function(object, locus) standardGeneric("trimToLocus"); }
     setGeneric("trimToLocus", fun);
}
setMethod("trimToLocus", "patientMatrix", function(object, locus){
   #print(sprintf("trimToLocus: my dataFile:%s", object@dataFile));
   # locus may be length 1 
   #  ex: "DRB1" will be converted to c("HLA.DRB1.1", "HLA.DRB1.2")
   # OR length > 1 ex: c("HLA.DRB1.1", "HLA.DRB1.2", ...)
   if(is.null(locus) || length(locus) == 0 || length(locus) > 2){
      stop("patientMatrix trimToLocus bad locus");
   }
   else if(length(locus) == 1){
      # expand it
      locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
      #print(sprintf("expanded locus to:%s",  paste(locus, collapse=' ')));
   }
   # locus should now contain the names of the rows that we want

   # Which columns contain this locus?
   columns <- which(locus[1] == colnames(object@dataMat) | 
                     locus[2] == colnames(object@dataMat));
   
   # Take the slice of just this locus.
   object@dataMat <- object@dataMat[,columns];

   object;
})

# patternCount
# Count the frequency of individuals with and without the pattern.
# See makePattern.
if (!isGeneric("patternCount")) {
     if (is.function("patternCount")){ fun <- patternCount; }
     else{ fun <- function(object, patternList, locus, seqMat) standardGeneric("patternCount"); }
     setGeneric("patternCount", fun);
}
setMethod("patternCount", "patientMatrix", function(object, patternList, locus, seqMat){
   if(is.null(locus) || length(locus) == 0){
      stop("patientMatrix patternCount no locus");
   }
   else if(length(locus) == 1){
      # expand it
      locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
   }
   # locus should now contain the names of the rows that we want

   # Find which columns have these loci names.
   columns <- which(locus == colnames(object@dataMat));
   if(length(columns) != 2){
      stop("patientMatrix patternCount not properly formatted");
   }

   pattern <- patternList[[1]];
   positions <- patternList[[2]];

   countWith <- 0;
   countWithout <- 0;
   for(row_index in 1:nrow(object@dataMat)){
      allele1 <- as.character(object@dataMat[row_index,columns[1]]); # DRB1*0301
      module1 <- c();
      if(!is.null(allele1) && !is.na(allele1) && nchar(allele1) > 0 &&
                                                allele1 %in% rownames(seqMat)){
         # Allele is ok, now check that there are no padding characters (zeros)
         # in the module positions.
         if(!any(seqMat[allele1, positions] == 0)){
            module1 <- paste(seqMat[allele1, positions], collapse='');
         }
         else{
            #module1 <- paste(seqMat[allele1, positions], collapse='');
            #print(sprintf("skipping %s", module1));
            module1<-c();
         }
      }

      allele2 <- as.character(object@dataMat[row_index,columns[2]]); # DRB1*0301
      module2 <- c();
      if(!is.null(allele2) && !is.na(allele2) && nchar(allele2) > 0 &&
                                                allele2 %in% rownames(seqMat)){
         # Allele is ok, now check that there are no padding characters (zeros)
         # in the module positions.
         if(!any(seqMat[allele2, positions] == 0)){
            module2 <- paste(seqMat[allele2, positions], collapse='');
         }
         else{
            #module2 <- paste(seqMat[allele2, positions], collapse='');
            #print(sprintf("skipping %s", module2));
            module2<-c();
         }
      }
      if(is.null(module1) && is.null(module2)){
         next;
      }
      if( length(grep(pattern=pattern, module1)) > 0 ||
            length(grep(pattern=pattern, module2)) > 0 ){
         countWith <- countWith + 1;
      }
      else{
         countWithout <- countWithout + 1;
      }

   }
   return(c(countWith, countWithout));
})


# Take a string of the form: "70:Q,71:!R,73:A"
# and convert to a pattern like this: "Q[^R]A"
# Return a list, pattern first, positions vector second.
makePattern <- function(in_str){
  positions <- c();
  pattern <- c();
  # remove spaces
  in_str <- gsub("[[:space:]]", "", in_str);
  # split on commas
  pairs <- unlist(strsplit(in_str, split=",", fixed=TRUE));
  for(apair in pairs){
     # split on colons
     parts <- unlist(strsplit(apair, split=":", fixed=TRUE));
     # Check that the first part contains only digits,
     # and the second part a single character, optionally preceded with '!'. 
     if( length(grep("^[[:digit:]]+$", parts[1])) == 0 ||
         length(grep("^!?[[:alpha:]]$", parts[2])) == 0 ){
        return(NULL);
     }
     positions <- c(positions, as.numeric(parts[1]));
     # Check for negation:
     if( length(grep("^![[:alpha:]]$", parts[2])) == 1 ){
        parts[2] <- paste(sub(pattern="^!", replacement="[^", parts[2]), "]",
            sep='', collapse='');
     }
     pattern <- c(pattern, parts[2]);
  }
  return(list(paste(pattern, sep='', collapse=''), positions));
}



