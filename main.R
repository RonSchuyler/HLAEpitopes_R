# main.R
# Main Analysis class.
#

# Make dummy classes for default arguments.
setClassUnion("NumOrNULL", c("numeric", "NULL"));
setClassUnion("ListOrNULL", c("list", "NULL"));
setClassUnion("NumOrNULLOrList", c("numeric", "NULL", "list"));

source("patientMatrix.R"); # patient data frame, representation of file contents
source("modules.R"); # module names, positions, counts, p-values...
source("readfiles.R"); # allelesInFile

library("tools"); # for file_path_sans_ext()


# Class Analysis
setClass("Analysis",
   representation(
      locus="character", # ex: "DRB1" or c("HLA.DRB1.1", "HLA.DRB1.2") 
      dataFile="character", 
      # If positions is not NULL, only look at this position subset.
      # If it is a 2 element list it is a pattern (see makePattern).
      positions="NumOrNULLOrList", 
      polyPos="NumOrNULL", 
      seqMat="matrix", # Same as seqHash, but as a named matrix.
      seqHash="list", # Hash of 0-padded sequences keyed on allele name.
      distMat="matrix", 
      controlPatients="patientMatrix", 
      affectedPatients="patientMatrix", 
      moduleSet="modules", # not used
      cluster="data.frame", # not used
      controlAlleles="list",
      affectedAlleles="list",
      uiStatusString="character", # put output filenames here
      FDR="numeric"), # False discovery rate.
   prototype(locus="DRB1", dataFile="../data/mhc.csv", positions=c())

);

# Analysis constructor
Analysis <- function(
      locus="DRB1", 
      dataFile="../data/mhc.csv", 
      positions=c(),
      FDR=.05,
      groupsOfN=c(),
      #doPatientCounts=TRUE,
      doCluster=TRUE,
      #doAlleleFreqs=FALSE,
      outputToScreen=TRUE,
      outputToFile=FALSE,
      printAllelesWithModules=FALSE,
      awmToFile=FALSE,
      uiStatusString="[error...arrrgh...]",
      a.version="0"
      ){
   #print("Analysis contstuctor");

   # Create the new object.
   this <- new("Analysis", locus=locus, dataFile=dataFile, positions=positions,
                                                                     FDR=FDR);
   this@uiStatusString <- c("Running 1");

   # Read the data file into affected and control patient matrices.
   this@affectedPatients <- patientMatrix(dataFile=dataFile, control=FALSE);
   this@controlPatients <- patientMatrix(dataFile=dataFile, control=TRUE);
   this@affectedPatients <- trimToLocus(this@affectedPatients, locus);
   this@controlPatients <- trimToLocus(this@controlPatients, locus);

   #if(doPatientCounts){
      # Get affected and control allele counts, counting patients.
      #this@affectedAlleles <- alleleCountHashPatients(this@affectedPatients, 
                                                   #locus=this@locus);
      #this@controlAlleles <- alleleCountHashPatients(this@controlPatients, 
                                                   #locus=this@locus);
   #}
   #else{
      # Get affected and control allele frequencies as hash lists.
      # Used below to get the padded sequences list.
      this@affectedAlleles <- alleleCountHash(this@affectedPatients, 
                                                   locus=this@locus);
      this@controlAlleles <- alleleCountHash(this@controlPatients, 
                                                   locus=this@locus);
   #}

   # This is the alleles file with sequence data for each allele.
   alleleFile <- "../AlleleImport2.txt";
   # Get hash and matrix of aligned padded sequences.
   this@seqHash <- get_padded_seqs(this@affectedAlleles, this@controlAlleles,
                                                         file_name=alleleFile);
   this@seqMat <- get_hash_values_as_matrix(this@seqHash);

   if(is.list(positions) && length(positions) == 2){
      #inputFileName <- substr(dataFile,9,nchar(dataFile));
      #fileBaseName <- substr(inputFileName, 1, (nchar(inputFileName)-4));
      inputFileName <- basename(dataFile);
      fileBaseName <- file_path_sans_ext(inputFileName);
      #toFile <- sprintf("../output/%s_%s_pattern_%s.txt", fileBaseName, this@locus, paste(unlist(positions), collapse="."));
      timestamp <- format(Sys.time(), format="%Y.%m.%d.%H%M%S");
      toFile <- sprintf("../output/%s_%s_pattern_%s.txt", fileBaseName, this@locus, timestamp);
      toFile <- gsub("^","not",toFile,fixed=TRUE);
      toFile <- gsub("[","",toFile,fixed=TRUE);
      toFile <- gsub("]","",toFile,fixed=TRUE);
      #print(sprintf("patternfile: %s", toFile));
      # This is a pattern.  Don't do the normal full analysis.
      doPatternAnalysis(this@affectedPatients, this@controlPatients, 
            pattern=positions, locus=this@locus, seqMat=this@seqMat, toFile=toFile);
      this@uiStatusString <- c("Finished.","Output file:",normalizePath(toFile));
      return(this);
   }

   # Get a list of polymorphic positions in this dataset.
   this@polyPos <- get_polys(this@seqMat);
   #print("poly:");
   #print(this@polyPos);

   # If positions is provided, the actual set of positions to check is
   # the intersection of positions and the set of polymorphic positions.
   # Otherwise, check all polymorphic positions.
   if(is.null(this@positions)){
      this@positions <- this@polyPos;
   }
   else{
      this@positions <- intersect(this@polyPos, this@positions);
   }

   # Get the count of modules.
   #this@moduleSet <- modules(this@affectedAlleles, this@controlAlleles, 
                           #this@seqHash, this@positions, FDR=this@FDR, 
                              #groupsOfN=groupsOfN);

   # Parse the input file name to build output file names, if necessary.
   # Previous: Assume input filename is "../data/*.csv"; 6/25/19 Don't assume.
   #inputFileName <- substr(dataFile,9,nchar(dataFile));
   #fileBaseName <- substr(inputFileName, 1, (nchar(inputFileName)-4));
   inputFileName <- basename(dataFile);
   fileBaseName <- file_path_sans_ext(inputFileName);
   toClusterFile <- c();
   logFile <- c(); 
   if(is.null(groupsOfN)){
      gon <- "All";
   }
   else{
      gon <- as.character(groupsOfN);
      groupsOfN <- eval(parse(text=groupsOfN));
   }

   if(identical(this@positions, this@polyPos)){
      pos <- "All";
   }
   else{
      pos <-  sprintf(paste(this@positions,collapse='.'));
   }
print("running");
   this@uiStatusString <- c("Running 2");
   if(outputToFile){ # Build filenames if necessary.
      # This should always be true, now that there is no screen output.
      # All modules output will go to this file.
      #toFile <- sprintf("../output/%s_%s_output%sof%s.txt", fileBaseName, this@locus, gon, pos);
      timestamp <- format(Sys.time(), format="%Y.%m.%d.%H%M%S");
      toFile <- sprintf("../output/%s_%s_output%s_%s.txt", fileBaseName, this@locus, gon, timestamp);
      if(awmToFile){
         # This includes all variable positions in file name.
         #awmFileName <- sprintf("../output/AWM_%s_%s_output%sof%s.txt", fileBaseName, this@locus, gon, pos);
         # Name w/o positions if too long:
         #awmFileName <- sprintf("../output/AWM_%s_%s_outputGroupsOf%s.txt", fileBaseName, this@locus, gon);
         # add timestamp
         timestamp <- format(Sys.time(), format="%Y.%m.%d.%H%M%S");
         # 10/22/2018: adding alleles to results file, so will ignore this anyway...
         awmFileName <- sprintf("../output/AWM_%s_%s_output%s_%s.txt", fileBaseName, this@locus, gon, timestamp);
      }
      # Counts of every module at each position will be output to this file.
      # This information is redundant with that in toFile, 
      # but in a different format.
      # Could be added to printModules call. Not used yet.
      #countOutFile <- sprintf("../output/%s_%s_Counts%sof%s.txt", fileBaseName,
         #this@locus, gon, pos);

      #toFile <- sprintf("../output/%s_output.csv", fileBaseName);
      if(doCluster){
         toClusterFile <- sprintf("../output/%s_%s_clusters%sof%s.txt", 
            fileBaseName, this@locus, gon, pos);
         #toClusterFile <- sprintf("../output/%s_clusters.csv", fileBaseName);
      }
   }
   else{
      toFile <- c();
      # This shouldn't happen now that there is no screen output.
      this@uiStatusString <- ""; 
   }

   this@uiStatusString <- c("Running 3");
   # Alleles that are skipped will be logged to this file.
   #skippedAllelesFile <- sprintf("../output/%s_SkippedAlleles_%s.txt", 
   #   fileBaseName, this@locus);
   timestamp <- format(Sys.time(), format="%Y.%m.%d.%H%M%S");
   logFile <- sprintf("../output/%s_%s_output%s_%s.log", fileBaseName,
         this@locus, gon, timestamp);

   logToFile(logFile, logmessage=sprintf("HLA Epitopes version %s", a.version), firsttime=TRUE);
   # Log the name of the input data file, the locus, and the allele file used.
   logToFile(logFile, logmessage=sprintf("Data file: %s", inputFileName), echo=TRUE);
   logToFile(logFile, logmessage=sprintf("Locus: %s", this@locus), echo=TRUE);
   logToFile(logFile, logmessage=sprintf("Positions:%s", paste(this@positions,collapse=' ')), echo=TRUE);
   logToFile(logFile, logmessage=sprintf("Groups of %s",groupsOfN), echo=TRUE);

   # Count the number of counted individuals.
   this@affectedPatients@patientCount <- countPatients(this@affectedPatients, 
                          this@locus, this@seqMat, logfile=logFile);
   this@controlPatients@patientCount <- countPatients(this@controlPatients, 
                          this@locus, this@seqMat, logfile=logFile);
   logToFile(logFile, logmessage="Finished checking alleles.", echo=TRUE);

   #print('', quote=FALSE);
   logToFile(logFile, sprintf("Counted %.1f affected and %.1f control subjects.", 
      this@affectedPatients@patientCount, this@controlPatients@patientCount), echo=TRUE);

   logToFile(logFile, logmessage="Testing epitopes...", echo=TRUE);

   this@uiStatusString <- c("Running 4");
   timediff <- printModules(toScreen=outputToScreen, toFile=toFile, doCluster=doCluster, clusterFile=toClusterFile, seqMat=this@seqMat, seqHash=this@seqHash, locus=this@locus, affectedPatients=this@affectedPatients, controlPatients=this@controlPatients, FDR=this@FDR, allPositions=this@positions, groupsOfN=groupsOfN, printAllelesWithModules=printAllelesWithModules, awmToFile=awmToFile);
   logToFile(logFile, logmessage=sprintf("Finished. Duration: %s", format(timediff,digits=3)), echo=TRUE);

   # Update the status string to be displayed in UI.
   this@uiStatusString <- c("Output files:",normalizePath(toFile));
   this@uiStatusString <- c(this@uiStatusString,normalizePath(logFile));
   if(doCluster){
      this@uiStatusString <- c(this@uiStatusString,normalizePath(toClusterFile));
   }
   if(awmToFile){
      #this@uiStatusString <- c(this@uiStatusString,normalizePath(awmFileName));
   }
   this@uiStatusString <- c(this@uiStatusString,sprintf("Finished. Running time: %s",format(timediff,digits=3)));

   #this@moduleSet <- modules2(locus=this@locus, affectedPatients=this@affectedPatients, controlPatients=this@controlPatients, FDR=this@FDR, allPositions=this@positions, groupsOfN=groupsOfN);

   this;
}

# deprecated. This is not used. 
getSeqMat <- function(locus="DRB1", dataFile="../data/mhc.csv", alleleFile="../AlleleImport.txt"){
   print("getSeqMat start");
   affectedPatients <- patientMatrix(dataFile=dataFile, control=FALSE);
   controlPatients <- patientMatrix(dataFile=dataFile, control=TRUE);
   affectedAlleles <- alleleCountHash(affectedPatients, 
                                                   locus=locus);
      controlAlleles <- alleleCountHash(controlPatients, 
                                                   locus=locus);
   seqHash <- get_padded_seqs(affectedAlleles, controlAlleles, file_name=alleleFile);
   seqMat <- get_hash_values_as_matrix(seqHash);
   print("getSeqMat done");
   return(seqMat);
}


# How many patients have each allele?
# h=1 : heterozygous
# h=2 : homozygous
# h=3 : either
# Return can be passed to fdra (in modules.R) for FDR adjustment.
countPatientsWithAlleles <- function(locus="DRB1", dataFile="../data/mhc.csv", h=1, alpha=0.05){
   if(h != 1 && h != 2 && h!= 3){
      stop("h must be 1 (for heterozygous) or 2 (for homozygous) or 3 (for either)");
   }
   locus <- c(sprintf("HLA.%s.1", locus), sprintf("HLA.%s.2", locus));
   affected <- patientMatrix(dataFile=dataFile, control=FALSE);
   control <- patientMatrix(dataFile=dataFile, control=TRUE);
   affectedList <- patientsWithAllele(affected, locus=locus, oneORtwo=FALSE);
   controlList <- patientsWithAllele(control, locus=locus, oneORtwo=FALSE);
   if(h < 3){
      acMat <- two_lists_as_mat(affectedList[[h]], controlList[[h]]);
   }
   else{
      bothAff <- add2Lists(affectedList[[1]], affectedList[[2]]);
      bothCon <- add2Lists(controlList[[1]], controlList[[2]]);
      acMat <- two_lists_as_mat(bothAff, bothCon);
   }
   acMat <- cbind(acMat, 0);
   acMat <- cbind(acMat, 0);
   for(rowi in 1:nrow(acMat)){
      acMat[rowi, 3] <- signif(pvalue(acMat[rowi, 1], acMat[rowi, 2]), 3);
      if(acMat[rowi, 3] <= alpha){
         acMat[rowi, 4] <- 1;
      }
   }
   colnames(acMat) <- c("Affected", "Control", "        p     ", "accepted");
   #colnames(acMat) <- c("Affected", "Control");
   return(acMat);
}

# Only run when a pattern is supplied in the Positions to Include textbox.
doPatternAnalysis <- function(affectedPatients, controlPatients, patternList,
                              locus, seqMat, toFile){
   print("doPatternAnalysis start");
   affectedCounts <- patternCount(affectedPatients, patternList, locus, seqMat);
   controlCounts <- patternCount(controlPatients, patternList, locus, seqMat);
   mat <- matrix(c(affectedCounts, controlCounts), nrow=2, byrow=TRUE);
   num <- sum(mat);
   sr <- rowSums(mat);
   sc <- colSums(mat);
   E <- outer(sr, sc, "*")/num;
   pvalue <- c();
   method <- c();
   if(any(E < 10)){
      # Use Fisher's exact test.
      method <- "Fisher's exact test";
      fe <- fisher.test(mat, conf.int=FALSE);
      pvalue <- fe$p.value;
   }
   else{
      # Use chi square test.
      method <- "Pearson's Chi-squared test";
      cs <- chisq.test(mat, correct=FALSE);
      pvalue <- cs$p.value;
   }
   print("", quote=FALSE);
   print("Pattern:", quote=FALSE);
   print(patternList[[1]], quote=FALSE);
   print("", quote=FALSE);
   print("Positions:", quote=FALSE);
   print(patternList[[2]], quote=FALSE);
   print("", quote=FALSE);
   print("", quote=FALSE);
   print("         Affected  Control", quote=FALSE);
   print(sprintf("With     %6i %6i", affectedCounts[1], controlCounts[1]), quote=FALSE);
   print(sprintf("Without  %6i %6i", affectedCounts[2], controlCounts[2]), quote=FALSE);
   print("", quote=FALSE);
   print(sprintf("%s p-value:", method), quote=FALSE);
   print(pvalue, quote=FALSE);
   print("", quote=FALSE);

   patfile <- file(toFile, open="wt");
   writeLines("Pattern:", con=patfile);
   writeLines(gsub("^","not",patternList[[1]],fixed=TRUE), con=patfile);
   writeLines("", con=patfile);
   writeLines("Positions:", con=patfile);
   writeLines(paste(patternList[[2]],collapse=" "), con=patfile);
   writeLines("", con=patfile);
   writeLines("", con=patfile);
   writeLines("         Affected  Control", con=patfile);
   writeLines(sprintf("With     %6i %6i", affectedCounts[1], controlCounts[1]), con=patfile);
   writeLines(sprintf("Without  %6i %6i", affectedCounts[2], controlCounts[2]), con=patfile);
   writeLines("", con=patfile);
   writeLines(sprintf("%s p-value:", method), con=patfile);
   writeLines(format(pvalue,digits=3), con=patfile);
   close(patfile);
   print("doPatternAnalysis done");
}

   



