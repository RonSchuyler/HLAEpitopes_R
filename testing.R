#
# testing.R
# 
# Scratch file for manual testing and development. 
# Not production code.
#


 locus="DQA1"
 dataSet <-"../Data/aaT1DGCscrubbedhighresdata.csv"
 dataFile=dataSet;
  inputFileName <- substr(dataFile,9,nchar(dataFile));
  fileBaseName <- substr(inputFileName, 1, (nchar(inputFileName)-4));


 positions <- "8:85"
 positionsTextbox <- eval(parse(text=positions))
   a <- Analysis(locus=locus, dataFile=dataSet, positions=positionsTextbox, 
      FDR=0.05, groupsOfN=2, doCluster=FALSE,  outputToScreen=FALSE,
      #outputToFile=FALSE,
      outputToFile=TRUE,
      printAllelesWithModules=FALSE,
      awmToFile=FALSE);


affectedPatients <-    patientMatrix(dataFile=a@dataFile, control=FALSE);
controlPatients    <- patientMatrix(dataFile=a@dataFile, control=TRUE);

affectedPatients <- trimToLocus(affectedPatients, locus);
controlPatients <- trimToLocus(controlPatients, locus);

affectedAlleles <- alleleCountHash(affectedPatients,
                                                   locus=locus);
controlAlleles <- alleleCountHash(controlPatients,
                                                   locus=locus);



   alleleFile <- "../AlleleImport2.txt";

   seqHash <- get_padded_seqs(affectedAlleles, controlAlleles,
                                                         file_name=alleleFile);
   seqMat <- get_hash_values_as_matrix(seqHash);


slotNames(affectedPatients)
#[1] "dataMat"      "dataFile"     "control"      "allOrNone"    "patientCount"


# TODO:
# need to check doPatternAnalysis()
# add all loci in alignments?
# ParseAlignments : build seqHash, replace *->0, .->0, - -> reference

