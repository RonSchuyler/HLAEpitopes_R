
#getSigPos <- function(fileName="../output/short.tab"){
getSigPos <- function(fileName="../output/mhc_DRB1_output4ofAll.csv"){
   con <- file(fileName, open="r");
   line <- readLines(con, 1);
   count <- 0;
   indexCount <- rep(0, 500);
   while(1){
      line <- readLines(con, 1);
      if(length(line) == 0){
         break;
      }
      count <- count + 1;
      elements <- unlist(strsplit(line, "\t", fixed=TRUE));
      ##print(elements[2]);
      positions <- as.numeric(unlist(strsplit(elements[2], ",", fixed=TRUE)));
      #print(indexCount[positions]);
      if(elements[7] == 1 && elements[3] > elements[4]){ 
         #print(positions);
         indexCount[positions] <- indexCount[positions] + 1;
      }
      #print(indexCount);
   }
   close(con);
   i <- which(indexCount > 0);
   plot(indexCount[i], main=sprintf("%i nonzero", length(i)));

   return(indexCount);
}


# This is looking for any alleles that share contain a susceptible module.
# This is including shared modules that occur by chance and is not what
# is intended.  It should not be used as is.
allelesInFile <- function(fileName="../output/short.tab", seqMat, alpha=c(),
   susceptible=TRUE){
   con <- file(fileName, open="r");
   line <- readLines(con, 1);
   count <- 0;
   nrow_seqMat <- nrow(seqMat);
   alleleNames <- rownames(seqMat);
   allelesHash <- list();
   firstHash <- list();
   bedone=FALSE;
   while(1){
      line <- readLines(con, 1);
      if(length(line) == 0){
         break;
      }
      count <- count + 1;
      elements <- unlist(strsplit(line, "\t", fixed=TRUE));
      #print(elements[2]);
      modules <- unlist(strsplit(elements[1], ",", fixed=TRUE));
      positions <- as.numeric(unlist(strsplit(elements[2], ",", fixed=TRUE)));
      affectedCount <- as.numeric(unlist(strsplit(elements[3],",",fixed=TRUE)));
      controlCount <- as.numeric(unlist(strsplit(elements[4],",",fixed=TRUE)));
      accepted <- as.numeric(unlist(strsplit(elements[7],",",fixed=TRUE)));
      pAdj <- as.numeric(unlist(strsplit(elements[8],",",fixed=TRUE)));

      if( (is.null(alpha) && accepted==1) || (!is.null(alpha) && pAdj<alpha) ){
         if( (susceptible && (affectedCount > controlCount)) ||
            (!susceptible && (controlCount > affectedCount)) ){ 
            # This module is accepted, get alleles containing it.
            for(modStr in unlist(strsplit(modules, split=" ", fixed=TRUE))){
               mod <- unlist(strsplit(modStr, split=c(), fixed=TRUE));
               for(rowi in 1:nrow_seqMat){
                  if(identical(seqMat[rowi, positions], mod)){
                     #allelesHash[[alleleNames[rowi]]] <- 1;
                     allelesHash[[alleleNames[rowi]]] <- sprintf("%s %s",
                        modStr, paste(positions, collapse=','));
                     if(!(alleleNames[rowi] %in% names(firstHash))){
                        firstHash[[alleleNames[rowi]]] <- sprintf("%s %s",
                           modStr, paste(positions, collapse=','));
                     }
                     #if(length(allelesHash) ==57){
                     if(FALSE){
                        print(sprintf("reached all at count:%i", count));
                        close(con);
                        return();
                     }
                     #print(sprintf("%s in %s", modStr, alleleNames[rowi]), 
                        #quote=FALSE);
                  }
               }

            }
         }
      }
   }
   close(con);
   filePathParts <- unlist(strsplit(fileName, split='/', fixed=TRUE));
   if(length(filePathParts) == 1){
      outputFile <- sprintf("alleles_%s", filePathParts);
   }
   else{
      outputFile <- sprintf("%s/alleles_%s", 
         paste(filePathParts[1:(length(filePathParts)-1)], collapse='/'),
            filePathParts[length(filePathParts)]);
   }
   con <- file(outputFile, open="w");
   for(allele in names(allelesHash)){
      #print(allele);
      line <- sprintf("%s %s %s", allele, firstHash[[allele]], allelesHash[[allele]]);
      writeLines(line, con);
      #writeLines(allele, con);
   }
   close(con);
}

# Open a connection to write alleles with module to a file.
# Output file name is constructed from fileName.
openAWMfile <- function(fileName){
   filePathParts <- unlist(strsplit(fileName, split='/', fixed=TRUE));
   if(length(filePathParts) == 1){
      outputFile <- sprintf("AWM_%s", filePathParts);
   }
   else{
      outputFile <- sprintf("%s/AWM_%s", 
         paste(filePathParts[1:(length(filePathParts)-1)], collapse='/'),
            filePathParts[length(filePathParts)]);
   }
   con <- file(outputFile, open="w");

}

# convertNom2to3
# Convert the data file using nomenclature v2 to v3:
#   v3                v2
# DRB4*01:03:01:02N DRB4*01030102N
# DRB4*02:01N DRB4*0201N
# DRB5*01:01  DRB5*0101
# DRB5*01:02  DRB5*0102
# DRB5*01:08N DRB5*0108N
# DRB5*02:02  DRB5*0202
# If the allele name ends in a letter then DO NOT truncate.
# If the allele has more than 1 colon then truncate everything 
# including and to the right of the second colon.
# Only works for mhc.csv format files, not COPD format.
convertNom2to3 <- function(filename="../data/mhc.csv", col=c(31,32)){
   mhc_table <- read.table(filename, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);
   for(col_i in col){
      for(row_i in 1:nrow(mhc_table)){
         allele <- mhc_table[row_i, col_i];
         if(nchar(allele) <= 7){
            next;
         }
         ss <- seq(8,nchar(allele),by=2);
         mhc_table[row_i, col_i] <- paste(substring(allele,1,7), 
                                  paste(substring(allele,ss,(ss+1)), collapse=":"), sep=":");
      }
   }
   outfile <- sprintf("%s%s%s", 
         substring(filename,1,(nchar(filename)-4)),
         "_nom3",
         substring(filename,(nchar(filename)-3)));

   write.table(mhc_table, file=outfile, quote=FALSE, sep=",", row.names=FALSE);

}

convertAlleleImport <- function(newAlleleFile="../AlleleImport2.txt"){
   # Match convfile col2 to AlleleImport col2, 
   # replace AlleleImport col2 with col1 from alleleConversion
   if(newAlleleFile == "../AlleleImport.txt"){
      stop("Don't overwrite original AlleleImport file! Pick a different output file name.");
   }

   notConvertedCount <- 0; # count skipped alleles not in conversion file
   convertedCount <- 0;
   usedFromConversionFile <- c();

   filename <- "../AlleleImport.txt";
   header_line <- readLines(filename,1);

   outfile <- file(newAlleleFile, open="wt"); # open and truncate
   writeLines(header_line, outfile);

   ccol <- 2; # column to convert
   allele_table <- read.table(filename, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);

   convfile <- "../alleleConversion.txt"; # Col1: new name; col2: old name
   conv_table <- read.table(convfile, header=TRUE, as.is=TRUE, strip.white=TRUE);
   for(row_i in 1:nrow(allele_table)){
      w_row <- which(allele_table[row_i,ccol] == conv_table[,2]);
      if(length(w_row)==0){
         print(sprintf("%s not in conversion table, skipping",allele_table[row_i,ccol]),quote=FALSE);
         notConvertedCount <- notConvertedCount + 1;
         next;
      }
      if(length(w_row) > 1){
         stop("%s has more than one entry in conversion file", allele_table[row_i,ccol]);
      }
      allele_table[row_i,ccol] <- conv_table[w_row,1];
      writeLines(sprintf("%i,\"%s\",\"%s\",%i,\"%s\",%i",allele_table[row_i,1],allele_table[row_i,2],allele_table[row_i,3],allele_table[row_i,4],allele_table[row_i,5],allele_table[row_i,6]),outfile);
      convertedCount <- convertedCount + 1;
      usedFromConversionFile <- c(usedFromConversionFile, w_row);
   }
   close(outfile);
   print(sprintf("converted %i allele names, skipped %i", convertedCount, notConvertedCount),
      quote=FALSE);
   return(usedFromConversionFile);
}

         



