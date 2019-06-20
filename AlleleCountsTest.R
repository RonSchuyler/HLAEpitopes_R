# AlleleCountsTest.R
# Functions to compute allele based tests.
# This is the "traditional" analysis, for comparison with the epitope 
# modules approach.
# See Calculating Gene Frequencies.doc
# Ron Schuyler 12/12/2008
#

# Main method is alleleTest at bottom.




# Use methods for reading data files and counting alleles from patientMatrix.R
source("patientMatrix.R");
# Also depends on hla3.R and hla2.R and must be in same directory.


# Merge matrices by row names.
# If a rowname does not exist in one of the matricies, it gets c(0, sum(col1)),
# or c(0, total_ind) if supplied.
# total_ind are the total number of individuals in that catagory.
cbindByName <- function(mat1, mat2, total_ind1=c(), total_ind2=c()){
   #allNames <- unique(c(rownames(mat1), rownames(mat2)));
   allNames <- sort(unique(c(rownames(mat1), rownames(mat2))));
   result <- matrix(0, nrow=length(allNames), ncol=4);
   rownames(result) <- allNames;
   if(is.null(total_ind1)){
      mat1total <- sum(mat1[,1]);
   }
   else{
      mat1total <- total_ind1;
   }
   if(is.null(total_ind2)){
      mat2total <- sum(mat2[,1]);
   }
   else{
      mat2total <- total_ind2;
   }
   for(aname in allNames){
      if(aname %in% rownames(mat1)){
         result[aname, 1:2] <- mat1[aname,];
      }
      else{
         result[aname, 1:2] <- c(0, mat1total);
      }
      if(aname %in% rownames(mat2)){
         result[aname, 3:4] <- mat2[aname,];
      }
      else{
         result[aname, 3:4] <- c(0, mat2total);
      }
   }
   return(result);
}


# Merge the values from two hashes into a 2 column matrix.
# Use hash keys for row names.
# Keys that don't exist in one matrix are given a 0 in that column.
# Result rows are ordered by name.
# Values for the first hash will be in column 1.
twoHashesToMatrix <- function(h1, h2){
   allKeys <- sort(unique(c( names(h1), names(h2))));
   result <- matrix(0, nrow=length(allKeys), ncol=2);
   rownames(result) <- allKeys;
   for(row_i in 1:nrow(result)){
      thisKey <- allKeys[row_i];
      # If the key exists in this hash, put the value into the result matrix.
      # Otherwise leave it as 0.
      if(thisKey %in% names(h1)){
         result[row_i, 1] <- h1[[thisKey]];
      }
      if(thisKey %in% names(h2)){
         result[row_i, 2] <- h2[[thisKey]];
      }
   }
   return(result);
}


# Compute probability of seeing this result by chance for each allele.
# Test used is fisher exact or chi square, as appropriate.
# Results are adjusted using FDR correction for multiple comparisons.
# Expected input is a matrix where row names are alleles:
#     affected gene frequency, not agf, control gf, not cgf.
testGFCounts <- function(countMatrix){
   result <- cbind(rownames(countMatrix), countMatrix, 
                        matrix(0,nrow=nrow(countMatrix),ncol=3));
   ci <- ncol(countMatrix) + 1;
   # What is in each column?
   methodCol <- ci + 1; # method: fisher exact or chi square
   testResultCol <- ci + 2; # pvalue from the test
   adjPcol <- ci + 3; # adjusted pvalue
   for(row_i in 1:nrow(countMatrix)){
      m <- matrix(as.numeric(countMatrix[row_i, ]), nrow=2, byrow=TRUE);
      n <- sum(m);
      sr <- rowSums(m);
      sc <- colSums(m);
      E <- outer(sr, sc, "*")/n; # Expected values.
      if(any(E < 10)){
         # Use Fisher's exact test.
         fe <- fisher.test(round(m), conf.int=FALSE);
         result[row_i, methodCol] <- "FisherExact";
         result[row_i, testResultCol] <- fe$p.value;
      }
      else{
         # Use chi square test.
         cs <- chisq.test(m, correct=FALSE);
         result[row_i, methodCol] <- "ChiSquare";
         result[row_i, testResultCol] <- cs$p.value;
      }
   }
   # Adjust for multiple comparisons using Benjamini and Yekutieli 
   # false discovery adjustment.
   result[, adjPcol] <- p.adjust(as.numeric(result[, testResultCol]), 
                                                               method="BY");
   colnames(result) <- c("Allele", "Affected GF Count", "Affected GF NOT Count", "Control GF Count", "Control GF NOT Count", "method", "p-value", "adjusted p-value");
   return(result);
}


# Calculate gene frequencies from phenotype frequencies.
# Input is a 2 column matrix where each row is a pair of alleles (one person).
# Output is a 2 column matrix where rows are alleles:
# column 1 = gene frequency count
# column 2 = not gene frequency count (col1 - sum(col1))
# row names are alleles.
calcGFCount <- function(dataMat){
   individual_count <- 0; # count the number of individuals
   allele_count <- 0; # count the number of alleles
   allele_hash <- list(); # hash of allele counts
   for(row_i in 1:nrow(dataMat)){
      a1 <- dataMat[row_i,1]; # allele 1
      a2 <- dataMat[row_i,2]; # allele 2
      if((a1 == '' || a1 == '-') && (a2 == '' || a2 == '-')){
         # If both alleles are blank or '-', skip this individual.
         next;
      }
      individual_count <- individual_count + 1;
      # Assume homozygosity.  If one allele is blank or '-', use the other.
      if(a1 == '' || a1 == '-'){
         a1 <- a2;
      }
      if(a2 == '' || a1 == '-'){
         a2 <- a1;
      }
      # Count the alleles.
      if(a1 != ''){
         allele_count <- allele_count + 1;
         if(a1 %in% names(allele_hash)){
            allele_hash[[a1]] <- allele_hash[[a1]] + 1;
         }
         else{
            allele_hash[[a1]] <- 1;
         }
      }
      if(a2 != ''){
         allele_count <- allele_count + 1;
         if(a1 != a2){
            if(a2 %in% names(allele_hash)){
               allele_hash[[a2]] <- allele_hash[[a2]] + 1;
            }
            else{
               allele_hash[[a2]] <- 1;
            }
         }
      }
   } # for each row
   oi <- order(names(allele_hash)); # order indicies
   j <- unlist(allele_hash)[oi]; # phenotype frequencies
   #k <-(1-sqrt(1-j/individual_count)) * allele_count;
   GFCount <-(1-sqrt(1-j/individual_count)) * individual_count * 2;
   # Rounding off to 1 decimal place makes the results match the sample data.
   # Don't round yet.
   #GFCount <- round(GFCount, 1); 
   NOTGFCount <- sum(GFCount) - GFCount;
   cbind(GFCount, NOTGFCount);
}
    

# Count individuals with each allele. 
# Argument is the dataMat (which is a data.frame) from a patientMatrix.
# Return a list: the first element is the hash: 
#     keys are alleles, 
#     values are counts of the number of individuals who carry that allele.
# The second element is the total number of counted individuals.
# Counts 1 or 2 copies, so sum may not match total individuals.
countIndwAllele <- function(patMat){
   individual_count <- 0; # count the number of individuals
   allele_count <- 0; # count the number of alleles
   allele_count_hash <- list(); # hash of allele counts
   for(row_i in 1:nrow(patMat)){
      a1 <- patMat[row_i,1]; # allele 1
      a2 <- patMat[row_i,2]; # allele 2
      if((a1 == '' || a1 == '-') && (a2 == '' || a2 == '-')){
         # If both alleles are blank or '-', skip this individual.
         next;
      }
      individual_count <- individual_count + 1;
      # Count the alleles.
      # Don't count blanks.
      if(a1 != '' && a1 != '-'){
         allele_count <- allele_count + 1;
         if(a1 %in% names(allele_count_hash)){
            allele_count_hash[[a1]] <- allele_count_hash[[a1]] + 1;
         }
         else{
            allele_count_hash[[a1]] <- 1;
         }
      }
      # Don't count blanks.
      if(a2 != a1 && a2 != '' && a2 != '-'){
         allele_count <- allele_count + 1;
         # Don't count homozygotes twice
         if(a1 != a2){
            if(a2 %in% names(allele_count_hash)){
               allele_count_hash[[a2]] <- allele_count_hash[[a2]] + 1;
            }
            else{
               allele_count_hash[[a2]] <- 1;
            }
         }
      }
   } # for each row

   # Make the hash a matrix, with counts in col1 and col2=individual_count-col1
   col1 <- unlist(allele_count_hash);
   col2 <- individual_count - col1;
   return_matrix <- cbind(col1, col2);

   return(list(return_matrix, individual_count));
}


# Main method for allele-based test.
# Input file should be in ../data directory.
# Result is alleleTestResult_inputFileName in ../data directory.
alleleTest <- function(infile="mhc.csv", locus="DRB1"){
   dataFile <- sprintf("../data/%s", infile);
   print(sprintf("alleleTest using data file %s", dataFile), quote=FALSE);
   print(sprintf("  using locus %s", locus), quote=FALSE);
   # Load data from data file at the given locus.
   affectedPatients <- patientMatrix(dataFile=dataFile, control=FALSE);
   controlPatients <- patientMatrix(dataFile=dataFile, control=TRUE);
   # Drop unneeded columns/rows.
   affectedPatients <- trimToLocus(affectedPatients, locus);
   controlPatients <- trimToLocus(controlPatients, locus);
   # Calculate gene frequencies.
   affectedGFcounts <- calcGFCount(affectedPatients@dataMat);
   controlGFcounts <- calcGFCount(controlPatients@dataMat);
   # Merge affected and control gene frequency counts.
   GFcountMatrix <- cbindByName(affectedGFcounts, controlGFcounts);
   # Compute probabilities.
   result <- testGFCounts(GFcountMatrix);
   #colnames(result) <- c("Allele", "Affected GF Count", "Affected GF NOT Count", "Control GF Count", "Control GF NOT Count", "method", "p-value", "adjusted p-value");
   result[,2:5] <- round(as.numeric(result[,2:5]),1); # round gene frequencies
   outfile <- sprintf("../data/alleleTestResult_%s_%s", locus, infile);
   write.table(result, file=outfile, quote=FALSE,
      row.names=FALSE, eol="\r\n", sep=",");

   # Now do the phenotype analysis.
   # Count the affected and control individuals with each allele.
   # Each _iwa is a 2 column matrix, counts w and counts wo allele.
   affected_iwa_and_total <- countIndwAllele(affectedPatients@dataMat);
   control_iwa_and_total <- countIndwAllele(controlPatients@dataMat);
   affected_iwa <- affected_iwa_and_total[[1]];
   affected_total_individuals <- affected_iwa_and_total[[2]];
   control_iwa <- control_iwa_and_total[[1]];
   control_total_individuals <- control_iwa_and_total[[2]];
   # Merge the two allele count hashes into a matrix.
   iwa_mat <- cbindByName(affected_iwa, control_iwa, 
                     affected_total_individuals, control_total_individuals);
   # Compute probabilities.
   pheno_result <- testGFCounts(iwa_mat);
   colnames(pheno_result) <- c("Allele", "Affected Individuals With", "Affected Individuals Without", "Control Individuals With", "Control Individuals Without", "method", "p-value", "adjusted p-value");
   pheno_outfile <- sprintf("../data/pheno_alleleTestResult_%s_%s", 
                                                            locus, infile);
   write.table(pheno_result, file=pheno_outfile, quote=FALSE,
      row.names=FALSE, eol="\r\n", sep=",");

}


