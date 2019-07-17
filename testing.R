#
# testing.R
# 
# Scratch file for manual testing and development. 
# Not production code.
#


 locus="DQA1"
 locus="DRB1"
 dataSet <-"../Data/aaT1DGCscrubbedhighresdata.csv";
 dataFile=dataSet;
  inputFileName <- substr(dataFile,9,nchar(dataFile));
  fileBaseName <- substr(inputFileName, 1, (nchar(inputFileName)-4));


 positions <- "45:53"
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

affectedAlleles <- alleleCountHash(affectedPatients, locus=locus);
controlAlleles <- alleleCountHash(controlPatients, locus=locus);



   alleleFile <- "../AlleleImport2.txt";

   seqHash <- get_padded_seqs(affectedAlleles, controlAlleles, file_name=alleleFile);
   seqMat <- get_hash_values_as_matrix(seqHash);


locusToGet <- unique(sub("\\*.*$", replacement="", x=names(seqHash))); # should only be one, but allow for multiple?


slotNames(affectedPatients)
#[1] "dataMat"      "dataFile"     "control"      "allOrNone"    "patientCount"


# TODO:
# need to check doPatternAnalysis()
# add all loci in alignments?
# ParseAlignments : build seqHash, replace *->0, .->0, - -> reference   DONE
# see get_padded_seqs()

seqHashNEW <- get_padded_seqs_from_alignments(affectedAlleles, controlAlleles, zipfile="../Alignments.zip");
allAllelesList <- get_padded_seqs_from_alignments(affectedAlleles, controlAlleles, zipfile="../Alignments.zip");
seqMatNEW <- get_hash_values_as_matrix(seqHashNEW);



#modules:
#moduleCountHash(controlPatients, locus, posV=posSet, seqMat); # this will skip a module containing 0's; OK
#fixClusterCounts(clusterSet, apm=affectedPatients@dataMat, cpm=controlPatients@dataMat, posV=posSet, seqHash=seqHash);
# fixClusterCounts() is ok because no modules sent to it from moduleCountHash() will include 0
#which_has_module(module=allModules[rowi,1], posV=allModules[rowi,2], seq_mat=seqMat)
#which_has_module() is ok because no modules sent to it from moduleCountHash() will include 0

aline=" DQA1*01:01:01:01         MIL NKALLLGALA ";

for(ln in 1:11){
   aline <- alignment[ln];
    spl <- unlist(strsplit(aline, split=" "));
   alleleName <- spl[2];
   sequence <- paste(spl[3:length(spl)], collapse="");
}



# finding closest matching allele: eg DRB1*07:01 -> use DRB1*07:01:01:01  ?

#"DRB1\*07:01" DRB_prot.txt | awk -F" " '{print $1}'|sort -u| awk '{printf "\"%s\", ", $1}'
d07 <- c("DRB1*07:01:01:01", "DRB1*07:01:01:02", "DRB1*07:01:01:03", "DRB1*07:01:01:04", "DRB1*07:01:02", "DRB1*07:01:03", "DRB1*07:01:04", "DRB1*07:01:05", "DRB1*07:01:06", "DRB1*07:01:07", "DRB1*07:01:08", "DRB1*07:01:09", "DRB1*07:01:10", "DRB1*07:01:11", "DRB1*07:01:12", "DRB1*07:01:13", "DRB1*07:01:14", "DRB1*07:01:15", "DRB1*07:01:16", "DRB1*07:01:17", "DRB1*07:01:18", "DRB1*07:01:19", "DRB1*07:01:20", "DRB1*07:01:21", "DRB1*07:01:22")
d07 <- c("DRB1*07:01:01:01", "DRB1*07:01:01:02", "DRB1*07:01:01:03", "DRB1*07:01:01:04", "DRB1*07:01:02", "DRB1*07:01:03", "DRB1*07:01:04", "DRB1*07:01:05", "DRB1*07:01:06", "DRB1*07:01:07", "DRB1*07:01:08", "DRB1*07:01:09", "DRB1*07:01:10", "DRB1*07:01:11", "DRB1*07:01:12", "DRB1*07:01:13", "DRB1*07:01:14", "DRB1*07:01:15", "DRB1*07:01:16", "DRB1*07:01:17", "DRB1*07:01:18", "DRB1*07:01:19", "DRB1*07:01:20", "DRB1*07:01:21", "DRB1*07:01:22", "DRB1*07:03", "DRB1*07:04", "DRB1*07:05", "DRB1*07:06", "DRB1*07:07", "DRB1*07:08", "DRB1*07:09");


ds <- "DRB1*07:01"

ds %in% d07

grep(ds, d07, fixed=TRUE)

# see "not found" in patientMatrix.R

