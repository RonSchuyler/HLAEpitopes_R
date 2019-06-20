#
# 10/1/07
# Functions to read data files and create hashes.
# 

# Read file mhc.csv
# Return a hash of alleleName:count
get_allele_counts_slow <- function(loci="HLA.DRB1", dx="Control", file_name="../data/mhc.csv"){
   # dx = diagnosis: "Affected" | "Control"
   # "Control" = "Matched Control" + "Random Control"
   mhc_table <- read.table(file_name, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);
   # Indices of rows with this Dx.
   rows <- grep(dx, mhc_table$Dx);
   countHash <- list();
   for(row_index in rows){
      for(loci_index in 1:2){
         locus_name <- paste(loci, ".", loci_index, sep=''); # HLA.DRB1.1
         loci_col <- which(names(mhc_table) == locus_name);
         allele <- as.character(mhc_table[row_index,loci_col]); # DRB1*0301
         #print(sprintf("allele:%s: row:%i,col:%i",
          #  allele, row_index, loci_index));
         if(allele %in% names(countHash)){
            countHash[[allele]] <- countHash[[allele]] + 1;
         }
         else{
            countHash[[allele]] <- 1;
         }
      }
   }
   return(countHash);
}

# Get allele counts for RA dataset.
# See print_allele_counts_RA
get_allele_counts <- function(loci="HLA.DRB1", dx="Control", file_name="../data/mhc.csv"){
   mhc_table <- read.table(file_name, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);
   rows <- grep(dx, mhc_table$Dx);
   countHash <- list();
   locus1 <- paste(loci, ".", 1, sep='');
   locus1_col <- which(names(mhc_table) == locus1);
   for(col in (locus1_col:(locus1_col+1))){
      #print(sprintf("column %i",col));
      for(row_index in rows){
         allele <- as.character(mhc_table[row_index, col]);
         if(nchar(allele, type="chars") < 9){
            # not a complete allele name, maybe DRB1*02
            # short names are not specific enough to get sequence data
            next;
         }
         #print(sprintf("%i,%i:  %s",row_index, col, allele));
         if(allele %in% names(countHash)){
            countHash[[allele]] <- countHash[[allele]] + 1;
         }
         else{
            countHash[[allele]] <- 1;
         }
      }
   }
   return(countHash);
}

# See print_allele_counts_copd
# control dx=0, affected dx=4
get_allele_counts_copd <- function(loci="DR", dx=0, file_name="../data/COPD_Cleaned_Up.csv"){
   #print(sprintf("get_a_c_copd: loci=%s, dx:%i, fn:%s",loci,dx,file_name));
   mhc_table <- read.table(file_name, header=TRUE, sep=',', as.is=TRUE, strip.white=TRUE);
   rows <- grep(dx, mhc_table$GOLD);
   countHash <- list();
   if(loci=="DR"){
      col1 <- 5;
   }
   else{
      col1 <- 7; # DQ
   }
   for(col in col1:(col1+1)){
      for(row_index in rows){
         allele <- as.character(mhc_table[row_index, col]);
         if(allele == '-'){
            # dash means same as previous column
            allele <- as.character(mhc_table[row_index, (col-1)]);
         }
         if(nchar(allele, type="chars") < 9){
            # not a complete allele name, maybe DRB1*02
            # short names are not specific enough to get sequence data
            next;
         }
         if( allele %in% names(countHash)){
            countHash[[allele]] <- countHash[[allele]] + 1;
         }
         else{
            countHash[[allele]] <- 1;
         }
      }
   }
   return(countHash);
}

# Print contents of affected and control count lists.
print_allele_counts <- function(aff, con){
   alleles <- sort(union(names(aff), names(con)));
   for(allele in alleles){
      if(allele %in% names(aff)){
         ac <- aff[[allele]];
      }
      else{
         ac <- 0;
      }
      if(allele %in% names(con)){
         cc <- con[[allele]];
      }
      else{
         cc <- 0;
      }
      pv <- pvalue(ac,cc);
      print(sprintf("%s %i %i   p=%g",allele,ac,cc,pv));
      #print(sprintf("%s %i %i",allele,ac,cc));
   }
}

# Return a 2-column matrix from two lists.
two_lists_as_mat <- function(list1, list2){
   allNames <- sort(union(names(list1), names(list2)));
   #mat <- matrix(0, nrow=length(allNames), ncol=2);
   mat <- c();
   currentRow <- 0;
   for(nm in allNames){
      currentRow <- currentRow+1;
      c1 <- 0;
      c2 <- 0;
      if(nm %in% names(list1)){
         #mat[currentRow,1] <- list1[[nm]];
         c1 <- list1[[nm]];
      }
      if(nm %in% names(list2)){
         c2 <- list2[[nm]];
         #mat[currentRow,2] <- list2[[nm]];
      }
      vec <- matrix(c(c1, c2), nrow=1);
      rownames(vec) <- nm;
      #names(vec) <- nm;
      #names(mat[currentRow,]) <- nm;
      mat <- rbind(mat, vec);
   }
   return(mat);
}

# Same results as:
#tmp<-p_correct(get_module_counts(get_polys(get_seq_mat(dataset="RA")),dataset="RA"))
print_allele_counts_RA <- function(loci="DR"){
   if(loci == "DR"){
      print("Allele counts for RA DR:");
      control <- get_allele_counts(loci="HLA.DRB1", dx="Control", file_name="../data/mhc.csv");
      affected <- get_allele_counts(loci="HLA.DRB1", dx="Affected", file_name="../data/mhc.csv");
   }
   else{
      print("Allele counts for RA DQ:");
      control <- get_allele_counts(loci="HLA.DQB1", dx="Control", file_name="../data/mhc.csv");
      affected <- get_allele_counts(loci="HLA.DQB1", dx="Affected", file_name="../data/mhc.csv");
   }
   #print_allele_counts(affected, control);
   allele_counts <- two_lists_as_mat(affected, control);
   p_correct(allele_counts);
   return(allele_counts);
}


print_allele_counts_copd <- function(loci="DR"){
   if(loci == "DR"){
      print("Allele counts for copd DR:");
      control <- get_allele_counts_copd(loci="DR", dx=0, file_name="../data/COPD_Cleaned_Up.csv");
      affected <- get_allele_counts_copd(loci="DR", dx=4, file_name="../data/COPD_Cleaned_Up.csv");
   }
   else{
      print("Allele counts for copd DQ:");
      control <- get_allele_counts_copd(loci="DQ", dx=0, file_name="../data/COPD_Cleaned_Up.csv");
      affected <- get_allele_counts_copd(loci="DQ", dx=4, file_name="../data/COPD_Cleaned_Up.csv");
   }
   #print_allele_counts(affected, control);
   allele_counts <- two_lists_as_mat(affected, control);
   p_correct(allele_counts);
   return(allele_counts);
}



# Get a list (hash) of aligned 0-padded sequences keyed by alleleName,
# given countHashes from get_allele_counts.  Actual counts are not used,
# the hash names are just used as a list of allele names to get.
# see get_hash_values_as_matrix()
get_padded_seqs <- function(affectedCounts, controlCounts, file_name="../AlleleImport.txt"){
   #print(sprintf("get_padded_seqs: %s", file_name));
   alleleNames <- union(names(affectedCounts), names(controlCounts));
   fileHandle <- file(file_name, open="r");
   header_line <- readLines(fileHandle,1);
   unpadded <- list();
   maxLen <- 0;
   while(1){
      line <- readLines(fileHandle,1);
      if(length(line) == 0){
         break;
      }
      element <- unlist(strsplit(line,",",fixed=TRUE));
      element <- gsub("\"","",element,fixed=TRUE);
      if(element[2] %in% alleleNames){
         # if the allele name is one we are looking for, add it to the list.
         offset <- as.numeric(element[4]);
         seq <- unlist(strsplit(element[5],split=""));
         if(offset >= 0){
            seq <- c(rep(0,offset),seq);
         }
         else{
            seq <- seq[(1 + offset*(-1)):length(seq)];
         }
         # list (hash) indexed by alleleName, seq is aligned at start,
         # possibly with different lengths
         unpadded[[element[2]]] <- seq; 
         if(length(seq) > maxLen){
            maxLen <- length(seq);
         }
      }
   }
   padded <- list();
   for(alleleName in names(unpadded)){
      # pad it and add it
      u_seq <- unpadded[[alleleName]];
      padded[[alleleName]] <- c(u_seq, rep(0,(maxLen-length(u_seq))));
   }
   close(fileHandle);
   return(padded);
}

# Similar to what as.matrix should do.
get_hash_values_as_matrix <- function(hash){
   mat <- c();
   for(key in sort(names(hash))){
      #value <- matrix(unlist(strsplit(hash[[key]], split="")),nrow=1);
      value <- hash[[key]];
      value <- matrix(value,nrow=1);
      rownames(value) <- key;
      if(is.null(mat)){
         mat <- matrix(value, nrow=1, ncol=length(value));
         rownames(mat) <- key;
      }
      else{
         mat <- rbind(mat, value);
      }
   }
   return(mat);
}
 

# posV (position vector) set of positions to check
# controlCounts and affectedCounts are hashes of counts, keyed on allele name
# paddedSeqs is a hash of 0-padded sequences, keyed on allele name
# Example setup:
# posV <- pocket1B;
# controlCounts <- get_allele_counts_copd(loci="DR", dx=0, file_name="../data/COPD_Cleaned_Up.csv");
# affectedCounts <- get_allele_counts_copd(loci="DR", dx=4, file_name="../data/COPD_Cleaned_Up.csv");
# paddedSeqs <- get_padded_seqs(affectedCounts, controlCounts);
#   OR
# controlCounts <- get_allele_counts(loci="HLA.DQB1", dx="Control", file_name="../data/mhc.csv");
# affectedCounts <- get_allele_counts(loci="HLA.DQB1", dx="Affected", file_name="../data/mhc.csv");
# paddedSeqs <- get_padded_seqs(affectedCounts, controlCounts);
get_position_counts <- function(posV, controlCounts, affectedCounts, paddedSeqs, doPrint=FALSE){
   allNames <- union(names(controlCounts), names(affectedCounts));
   affected <- list();
   control <- list();
   posV <- sort(posV);
   for(allele in allNames){
      if(allele %in% names(paddedSeqs)){
         if(any(paddedSeqs[[allele]][posV] == 0)){
            # 0 is the padding character, meaning this allele doesn't 
            # have complete sequence data, so we don't count it.
            next;
         }
         # module is the string of amino acids at the posV set of positions
         # for this allele.  Multiple alleles may have the same module.
         module <- paste(paddedSeqs[[allele]][posV], collapse='');
         if(module %in% names(affected)){
            if(allele %in% names(affectedCounts)){
               affected[[module]] <- affected[[module]] + affectedCounts[[allele]];
            }
         }
         else{
            if(allele %in% names(affectedCounts)){
               affected[[module]] <- affectedCounts[[allele]];
            }
         }
   
         if(module %in% names(control)){
            if(allele %in% names(controlCounts)){
               control[[module]] <- control[[module]] + controlCounts[[allele]];
            }
         }
         else{
            if(allele %in% names(controlCounts)){
               control[[module]] <- controlCounts[[allele]];
            }
         }
      }
   }

   positionsString <- sprintf(paste(posV,collapse=','));
   loop_count <- 0;
   pvaluesVec <- c();
   #allInfoList <- c(); # list of all info, same order as pvaluesVec
   # m=module key for counts
   for(m in sort(union(names(control), names(affected)))){
      if(m %in% names(control)){
         count_c <- control[[m]];
      }
      else{
         count_c <- 0;
      }
      if(m %in% names(affected)){
         count_a <- affected[[m]];
      }
      else{
         count_a <- 0;
      }
      pv <- pvalue(count_a, count_c);
      loop_count <- loop_count + 1;
      #allInfoList[[loop_count]] <- sprintf("(%s): %s  %i:%i   p=%g", 
            #positionsString, m, count_a, count_c, pv);
      if(doPrint){
         print(sprintf("%s:%i:%i   p=%g", m, count_a, count_c, pv));
      }
      pvaluesVec <- c(pvaluesVec, pv);
   }
   #return(list(affected, control, pvaluesVec, allInfoList));
   return(list(affected, control, pvaluesVec));

}

      
# Pocket residue positions for MHC I.
# Note that some positions may belong to more than one pocket.
# (Up to 84 is domain 1, 97 and up is domain 2.)
# (From HistoCheck website.)
pocketA <- c(5, 7, 59, 63, 66, 70, 99, 159, 163, 167, 171);
pocketB <- c(7, 9, 24, 25, 34, 45, 63, 66, 67, 70, 99);
pocketC <- c(9, 70, 73, 74, 97);
pocketD <- c(99, 114, 156, 159, 160);
pocketE <- c(97, 114, 133, 147, 152, 156);
pocketF <- c(77, 80, 81, 84, 116, 123, 143, 146, 147);
pocketE <- c(97, 114, 133, 147, 152, 156);

# Pocket residue positions for MHC II beta.      
# From Fu: "Pocket 4 of the HLA-DR(a,B1*0401) Molecule Is A Major Determinant 
# of T Cell Recognition of Peptide"
pocket1B <- c(85,86,89,90);
pocket4B <- c(13,70,71,74,78);
pocket6B <- c(11,13);
pocket7B <- c(28,47,61,67,71);
pocket9B <- c(9,57);

# Pocket residue positions for MHC II beta (DRB1) w/i 5A of peptide.
# From Floudas: "A Predictive Method for the Evaluation of Peptide Binding 
# in Pocket 1 of HLA-DRB1 via Global Minimization of Energy Interactions"
# More positions that contact peptide or TCR are listed on HistoCheck,
# but not grouped into pockets: 30, 32, 37, 38, 56, 64, 65, 66, 68, 69, 77, 81
pocket1B <- c(82, 85, 86, 89, 90); # +82
pocket4B <- c(13, 26, 70, 71, 74, 78); # +26
pocket6B <- c(11, 13, 71); # +71
pocket7B <- c(28, 47, 61, 67, 71);
pocket9B <- c(9, 57, 60, 61); # +60, 61   

# Additional positions for DQ?
# DQ 9:37
# DQ 4:28
# From Baas:"Peptide binding motifs and specificities for HLA-DQ molecules"
# DQ 1:86,87
# DQ 2:77
# DQ 3:74
# DQ 4:13,26,28,71,74
# DQ 5:70,71,74
# DQ 6:9,30,70
# DQ 7:30,47,67,70,71
# DQ 9:9,30,37,38,57,59?

pocketResiduesClassI <- c(pocketA,pocketB,pocketC,pocketD,pocketE,pocketF);
pocketResiduesClassII <- c(pocket1B,pocket4B,pocket6B,pocket7B,pocket9B);


source("hla2.R");
test_subsets <- function(posV=pocketResiduesClassII,dataset="copd", loci="DR"){
   # loci="HLA.DRB1"
   # loci="HLA.DQB1"
   posV <- sort(posV);
   print("From positions (posV):");
   print(posV);
   if(dataset == "copd"){
      cc <- get_allele_counts_copd(dx=0, loci=loci); # control
      ac <- get_allele_counts_copd(dx=4, loci=loci); # affected
   }
   else{
      if(loci == "DR"){
         loci <- "HLA.DRB1";
      }
      else if(loci == "DQ"){
         loci <- "HLA.DQB1";
      }

      cc <- get_allele_counts(loci=loci, dx="Control");
      ac <- get_allele_counts(loci=loci, dx="Affected");
   }

   padded_seq_hash <- get_padded_seqs(ac, cc);
   padded_seq_mat <- get_hash_values_as_matrix(padded_seq_hash);
   # Find the polymorphic positions in this dataset.
   my_polys <- get_polys(padded_seq_mat); 
   test_pos <- intersect(my_polys, posV);
   print("intersection of posV and polymorphic positions in this dataset:");
   print(test_pos);
   return(test_pos);
}


# N of X
# Return all possible combinations of n elements of vector x.
# The returned list will have choose(length(x),n) vectors.
nofx <- function(n, x){
   if(is.null(n) || is.na(n) || length(n)==0 || n==0){
      n <- length(x);
   }
   if(n > length(x)){
      warning(sprintf("nofx: n=%i > length(x)=%i", n, length(x)), 
         immediate.=TRUE);
      n <- length(x);
   }

   x <- sort(x);
   Qlist <- list();
   len <- length(x);
   posMat <- matrix(0,nrow=n,ncol=(len-n+1));
   jj <- c();
   for(i in 1:n){
      jj[i] <- 1;
      for(j in i:(len-n+i)){
         #print(sprintf("i:%i, j:%i",i,j));
         posMat[i,(j-i+1)] <- x[j];
         #group[i] <- x[j];
         #print(sprintf("%s",group));
      }
   }
   #print("posMat:");
   #print(posMat);


   Q <- c(); # results
   done <- FALSE;
   while(!done){
      save <- TRUE;
      for(pos in 1:n){
         if(pos > 1 && jj[pos] < jj[pos-1]){
            save <- FALSE;
            break;
         }
         Q[pos] <- posMat[pos, jj[pos]];
      }
      if(save){
         #print(Q);
         Qlist <- c(Qlist, list(Q));
      }
      for(ni in n:1){
         jj[ni] <- jj[ni] + 1;
         if(jj[ni] > (len-n+1)){
            if(ni == 1){
               done <- TRUE;
               break;
            }
            jj[ni] <- jj[ni-1];
            #jj[ni] <- jj[ni-1] + 1;
            #print(sprintf("set jjni: ni:%i, jjni:%i",ni,jj[ni]));
         }
         else{
            break;
         }
      }
   }
   return (Qlist);
}

   
# Test all combinations of n elements of the position vector posSet.
# Only polymorhic positions are used.
# posSet is optional.  If it is supplied, the values actually used are the
# intersection of posSet and all polymorphic positions in this dataset.
# If posSet is not supplied, all polymorphic positions are used.
testSetsofN <- function(n, posSet=c(), dataset="copd", loci="DR", ac=c(), cc=c(), FDR=.1){
   if(is.null(ac) || is.null(cc)){
      if(dataset=="copd"){
         cc <- get_allele_counts_copd(dx=0, loci=loci);
         ac <- get_allele_counts_copd(dx=4, loci=loci);
      }
      else{
         if(loci=="DR"){
            loci <- "HLA.DRB1";
         }
         else if(loci=="DQ"){
            loci <- "HLA.DQB1";
         }
         cc <- get_allele_counts(loci=loci, dx="Control");
         ac <- get_allele_counts(loci=loci, dx="Affected");
      }
   }

   padded_seq_hash <- get_padded_seqs(ac, cc);
   padded_seq_mat <- get_hash_values_as_matrix(padded_seq_hash);
   # Find the polymorphic positions in this dataset.
   my_polys <- get_polys(padded_seq_mat); 

   if(!is.null(posSet)){
      posSet <- intersect(my_polys, posSet);
   }
   else{
      posSet <- my_polys;
   }

   #print("print allele counts:");
   #print_allele_counts(ac,cc);
   setList <- nofx(n, posSet);
   #print("");
   allPvalues <- c();
   for(listI in 1:length(setList)){
      #print(sprintf("combination %i, %s", listI, paste(setList[[listI]], collapse=' ')));
      acl <- get_position_counts(posV=setList[[listI]], controlCounts=cc, 
               affectedCounts=ac, paddedSeqs=padded_seq_hash);
      allPvalues <- c(allPvalues, acl[[3]]);
   }
   isSig <- fdr_dep(allPvalues, FDR); # returns 3 column matrix
   sumSig <- sum(isSig[,3]==1);
   if(sumSig > 0){
      print(sprintf("n=%i:  %i of %i significant p-values at fdr=%.3f", n, sumSig, nrow(isSig), FDR));
   }
   return(allPvalues);
}


test_pocket <- function(posV, loci="DR", dataset="copd", affectedCounts=c(),
   controlCounts=c(), silent=TRUE){
   if(is.null(controlCounts) || is.null(affectedCounts)){
      if(dataset=="copd"){
         cc <- get_allele_counts_copd(dx=0, loci=loci);
         ac <- get_allele_counts_copd(dx=4, loci=loci);
      }
      else{
         if(loci=="DR"){
            loci <- "HLA.DRB1";
         }
         else if(loci == "DQ"){
            loci <- "HLA.DQB1";
         }
         cc <- get_allele_counts(loci=loci, dx="Control");
         ac <- get_allele_counts(loci=loci, dx="Affected");
      }
   }
   padded_seq_hash <- get_padded_seqs(ac, cc);
   padded_seq_mat <- get_hash_values_as_matrix(padded_seq_hash);
   # Find the polymorphic positions in this dataset.
   my_polys <- get_polys(padded_seq_mat); 

   posV <- intersect(posV, my_polys);
   acl <- get_position_counts(posV, controlCounts=cc, affectedCounts=ac,
                        paddedSeqs=padded_seq_hash,doPrint=!silent);
   pvalues <- acl[[3]];
   fdr_i <- fdr_index(pvalues, .1);
   if(!is.null(fdr_i)){
      print("Significant");
      for(i in 1:length(fdr_i)){
         print(acl[[4]][fdr_i[i]]);
      }
   }
   else{
      print("not significant");
   }
   return(pvalues);
}


wrapper <- function(loci="DR", dataset="copd", posV=c(), FDR=.1, maxN=5){
   allPV <- c();
   if(dataset=="copd"){
      cc <- get_allele_counts_copd(dx=0, loci=loci);
      ac <- get_allele_counts_copd(dx=4, loci=loci);
   }
   else{
      if(loci=="DR"){
         loci <- "HLA.DRB1";
      }
      else if(loci == "DQ"){
         loci <- "HLA.DQB1";
      }
      cc <- get_allele_counts(loci=loci, dx="Control");
      ac <- get_allele_counts(loci=loci, dx="Affected");
   }

   padded_seq_hash <- get_padded_seqs(ac, cc);
   padded_seq_mat <- get_hash_values_as_matrix(padded_seq_hash);
   # Find the polymorphic positions in this dataset.
   posSet <- get_polys(padded_seq_mat); 
   if(!is.null(posV)){
      posSet <- intersect(posSet, posV);
   }
   if(is.null(maxN) || maxN > length(posSet)){
      maxN <- length(posSet);
   }
   for(n in 1:maxN){
      print(sprintf("n=%i",n));
      allPV <- c(allPV, testSetsofN(n, posSet,loci=loci,dataset=dataset, ac=ac, cc=cc, FDR=FDR));
   }
   return(allPV);
}


pocketList <- c();
pocketList[[1]] <- pocket1B;
pocketList[[2]] <- pocket4B;
pocketList[[3]] <- pocket6B;
pocketList[[4]] <- pocket7B;
pocketList[[5]] <- pocket9B;
pocketNum <- c(1,4,6,7,9);
eachPocket <- function(loci="DR", dataset="copd", silent=TRUE){
   allPV <- c();
   for(i in 1:length(pocketList)){
      print(sprintf("Pocket %i:", pocketNum[i]));
      ps <- test_pocket(posV=pocketList[[i]],loci=loci,dataset=dataset,
                        silent=silent);
      allPV <- c(allPV, ps);
   }
   fdr(allPV, .15);
   return(allPV);
}

testp <- function(){
   print("copd DQ");
   copdq <- eachPocket(loci="DQ", dataset="copd", silent=FALSE);
   print("copd DR");
   copdr <- eachPocket(loci="DR", dataset="copd", silent=FALSE);
   print("RA DQ");
   radq <- eachPocket(loci="DQ", dataset="RA");
   print("RA DR");
   radr <- eachPocket(loci="DR", dataset="RA");
}

# Get the matrix of sequences for the given loci and dataset.
# Rows are sequences, columns are aligned, sequences are 0-padded.
# see also get_padded_seqs
get_seq_mat <- function(loci="DR", dataset="copd"){
   print(sprintf("%s %s", dataset, loci));
   if(dataset == "copd"){
      cc <- get_allele_counts_copd(dx=0, loci=loci); # control
      ac <- get_allele_counts_copd(dx=4, loci=loci); # affected
   }
   else{
      if(loci == "DR"){
         loci <- "HLA.DRB1";
      }
      else if(loci == "DQ"){
         loci <- "HLA.DQB1";
      }
      cc <- get_allele_counts(loci=loci, dx="Control");
      ac <- get_allele_counts(loci=loci, dx="Affected");
   }
   mat <- get_hash_values_as_matrix(get_padded_seqs(ac, cc));
   return(mat);
}

# Calculate p-values and test for significance with fdr.
# counts is a 2-column matrix
# See two_lists_as_mat() to convert 2 count lists into a 2-column matrix.
# counts <- get_module_counts(...);
# Does not actually correct p-values for multiple comparisons.
p_correct <- function(counts, printAccepted=TRUE, orderByDiff=TRUE, FDR=.05,
                                                significantFigures=4,
                                                n_affected=50, n_control=50){
   p_values <- rep(0, nrow(counts));
   accepted <- rep(0, nrow(counts));
   for(i in 1:length(p_values)){
      p_values[i] <- signif(pvalue_ue2(counts[i,1], counts[i,2], 
                              n_affected=n_affected, n_control=n_control), 
                                                      significantFigures);
   }
   fdri <- fdr_dep_index(p_values, r=FDR);
   accepted[fdri] <- 1;
   mat <- cbind(counts, p_values, accepted);
   colnames(mat)[1] <- "Affected";
   colnames(mat)[2] <- "Control";
   if(orderByDiff && nrow(mat) > 1){
      sort_order <- sort((counts[,1]-counts[,2]), decreasing=TRUE, index.return=TRUE);
      mat <- mat[sort_order$ix,];
   }
   if(printAccepted){
      w <- which(mat[,4] == 1);
      print(mat[w,]);
   }
   return(mat);
}


# Risler amino acid distance matrix.
aadm <- matrix(c(
rep(0,1),92,50,12,40,39,71,13,21,22,29,23,61,10,17,4,7,6,78,49,
rep(0,2),98,93,95,98,99,96,94,92,95,94,99,90,93,88,91,90,100,83,
rep(0,3),30,63,65,88,56,53,60,67,36,86,41,57,38,55,54,89,66,
rep(0,4),40,47,71,18,21,32,41,19,58,2,8,9,14,15,81,50,
rep(0,5),65,83,30,53,31,59,46,83,37,46,42,48,34,77,4,
rep(0,6),86,56,58,59,64,51,85,51,52,37,49,52,87,61,
rep(0,7),75,80,78,84,63,96,68,64,66,77,73,97,74,
rep(0,8),29,3,32,33,70,20,20,14,16,1,73,45,
rep(0,9),38,44,31,72,13,3,19,26,25,82,43,
rep(0,10),9,36,76,27,24,23,26,5,76,43,
rep(0,11),54,84,25,28,39,35,24,87,60,
rep(0,12),79,16,24,7,28,27,82,57,
rep(0,13),69,62,62,68,69,96,85,
rep(0,14),5,10,12,17,80,42,
rep(0,15),6,8,18,75,35,
rep(0,16),2,11,74,45,
rep(0,17),15,80,47,
rep(0,18),72,48,
rep(0,19),70,
rep(0,20)), 
nrow=20, ncol=20, byrow=TRUE);

# Calculate the distance between the two strings using the distance matrix dm.
distance <- function(string1, string2, dm=aadm){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=''));
   s1 <- unlist(strsplit(string1, split=''));
   s2 <- unlist(strsplit(string2, split=''));
   len <- min(length(s1), length(s2));
   total_dist = 0;
   for(pos in 1:len){
      posi <- c(which(s1[pos]==AAvec), which(s2[pos]==AAvec));
      r_i <- min(posi)
      c_i <- max(posi);
      #print(sprintf("%s,%s  :  %i", s1[pos], s2[pos], dm[r_i,c_i]));
      total_dist <- total_dist + dm[r_i, c_i];
   }
   return(total_dist);
}


get_module_counts <- function(posV=pocket4B, dataset="copd", loci="DR"){
   # loci="HLA.DRB1"
   # loci="HLA.DQB1"
   posV <- sort(posV);
   print("From positions (posV):");
   print(posV);
   if(dataset == "copd"){
      cc <- get_allele_counts_copd(dx=0, loci=loci); # control
      ac <- get_allele_counts_copd(dx=4, loci=loci); # affected
   }
   else{
      if(loci == "DR"){
         loci <- "HLA.DRB1";
      }
      else if(loci == "DQ"){
         loci <- "HLA.DQB1";
      }

      cc <- get_allele_counts(loci=loci, dx="Control");
      ac <- get_allele_counts(loci=loci, dx="Affected");
   }
   padded_seq_hash <- get_padded_seqs(ac, cc);
   #padded_seq_mat <- get_hash_values_as_matrix(padded_seq_hash);

   pc <- get_position_counts(posV, cc, ac, padded_seq_hash, doPrint=FALSE);
   # pc is return(list(affected(hash), control(hash), pvaluesVec, allInfoList));
   module_counts <- two_lists_as_mat(pc[[1]], pc[[2]]);
   return(module_counts);
}

# Return a square, symmetric distance matrix.
pair_dists <- function(posV=pocket4B, dataset="copd", loci="DR"){
   module_counts <- get_module_counts(posV, dataset, loci);
   #for(row_1 in 1:(nrow(module_counts)-1)){
      #for(row_2 in (1+row_1):nrow(module_counts)){
   # Make a square matrix of zeros.
   my_dm <- matrix(0, nrow=nrow(module_counts), ncol=nrow(module_counts));
   for(row_1 in 1:nrow(module_counts)){
      for(row_2 in 1:nrow(module_counts)){
         #print(sprintf("row1:%i, row2:%i", row_1, row_2));
         my_dist <- distance(rownames(module_counts)[row_1], 
                              rownames(module_counts)[row_2]);
         my_dm[row_1, row_2] <- my_dist;
         print(sprintf("%s to %s: %i", rownames(module_counts)[row_1],
               rownames(module_counts)[row_2], my_dist));
      }
   }
   rownames(my_dm) <- rownames(module_counts);
   colnames(my_dm) <- rownames(module_counts);
   return(my_dm);
}

# Return the distance between the two clusters.
# l1 and l2 are lists, possibly singletons.
# Distance is average distance between all cluster members.
clust_dist <- function(l1, l2){
   total_dist <- 0;
   count <- 0;
   if(length(l1) == 0 || length(l2)==0){
      return(0);
   }
   for(i_1 in 1:length(l1)){
      for(i_2 in 1:length(l2)){
         count <- count + 1;
         total_dist <- total_dist + distance(l1[[i_1]], l2[[i_2]]);
         #print(sprintf("count:%i td:%.3f", count, total_dist));
      }
   }
   #print(sprintf("final count:%i td:%.3f", count, total_dist));
   return(total_dist / count);
}


# Return a vector of values from the list l1.
# Similar to names() function.
values <- function(l1){
   vec <- c();
   for(pos in 1:length(l1)){
      #print(l1[[pos]]);
      vec <- c(vec, paste(l1[[pos]], collapse=' '));
   }
   return(vec);
}


# cluster row names of counts
# counts is a 2-d matrix of module counts (from get_module_counts())
# thresh is the minimum distance to merge clusters.
# limit is the maximum number of mergers to perform.
cluster <- function(counts, thresh=15, limit=Inf, quiet=TRUE){
   # Make a list of rownames.
   # Each cluster is intially a singleton.
   clusters <- c();
   for(rn in rownames(counts)){
      clusters <- c(clusters, list(rn));
   }
   minDist <- Inf;
   #while(minDist < thresh){
   loop_count <- 0;
   while(TRUE){
      loop_count <- loop_count+1;
      # Find the two nearest clusters.
      for(ind_1 in 1:(length(clusters)-1)){
         for(ind_2 in (ind_1+1):length(clusters)){
            #print(sprintf("%i %i   : len:%i",ind_1, ind_2, length(clusters)));
            d <- clust_dist(clusters[[ind_1]], clusters[[ind_2]]);
            if(!quiet){
               print(sprintf("%i %s,  %i %s   : %.3f", ind_1, 
                  paste(clusters[[ind_1]], collapse=' '),
                     ind_2, paste(clusters[[ind_2]], collapse=' '), d));
            }
            if(d < minDist){
               minDist <- d;
               minDisti <- c(ind_1, ind_2);
            }
         }
      }
      #print(sprintf("minDist: %.3f", minDist));
      if(minDist < thresh){
         # Merge the two nearest clusters.
         counts <- merge_counts(counts, minDisti);
         clusters <- merge_clust(clusters, minDisti);
         minDist <- Inf;
      }
      else{
         break;
      }
      # Stop if the limit is reached, or if this is the only cluster.
      if(loop_count >= limit || length(clusters)==1){
         break;
      }
   }

   rownames(counts) <- values(clusters);
   return(counts);
}


# Merge counts.
# counts is a 2-d matrix.
# Returned matrix will have 1 less row than mat.
merge_counts <- function(mat, indices){
   indices <- sort(indices);
   #print(sprintf("merge_counts: %i, %i", indices[1], indices[2]));
   #print(mat[indices[1],]);
   #print(mat[indices[2],]);
   newMat <- matrix(0, nrow=(nrow(mat)-1), ncol=ncol(mat)); 
   rowCount <- 0;
   for(matRow in 1:nrow(mat)){
      #print(sprintf("matRow:%i",matRow));
      if(matRow != indices[2]){
         rowCount <- rowCount + 1;
         newMat[rowCount,] <- mat[matRow,];
      }
      else{
         newMat[indices[1],] <- newMat[indices[1],] + mat[matRow,];
      }
      #print(sprintf("to row:%i",rowCount));
      #print(newMat);
   }
   return(newMat);
}


# Merge the cluster elements 1 and 2.
# clusters is a list, elements is a vector of 2 indices to merge.
merge_clust <- function(clusters, elements){
   elements <- sort(elements);
   #print(sprintf("merge_clust %i, %i", elements[1], elements[2]));
   newClusts <- c();
   for(clusti in 1:length(clusters)){
      if(clusti == elements[1]){
         newClusts <- c(newClusts, list(c(clusters[[elements[1]]], clusters[[elements[2]]])));
      }
      else if(clusti != elements[2]){
         newClusts <- c(newClusts, list(clusters[[clusti]]));
      }
   }
   return(newClusts);
}

# Return allele names containing the given module at the given set of positions.
# Either a sequnce matrix or a loci+dataset must be provided.
which_has_module <- function(module, posV, seq_mat=c(), loci=c(), dataset=c()){
   # Make sure we have a sequence matrix.
   if(is.null(seq_mat)){
      seq_mat <- get_seq_mat(loci=loci, dataset=dataset);
   }
   # Break the module string in a character vector.
   mod_vec <- unlist(strsplit(module,split=c(),fixed=TRUE));
   # Break the position vector into numbers if it is characters.
   if(is.character(posV)){
      posV <- as.numeric(unlist(strsplit(posV, split=',', fixed=TRUE)));
   }
   #alleles <- c();
   alleles <- "";
   # Inspect each allele at the given position vector.
   for(rowi in 1:nrow(seq_mat)){
      #print(sprintf("look at row %i: %s",rowi,seq_mat[rowi,posV]));
      if(length(posV) > 1){
         if(identical(seq_mat[rowi,posV], mod_vec)){
            #print(sprintf("allele: %s at row: %i", rownames(seq_mat)[rowi], rowi));
            #print(rownames(seq_mat)[rowi]);
            alleles <- paste(alleles, rownames(seq_mat)[rowi], sep="   ", 
                                                                  collapse='');
         }
      }
      else{
         if(seq_mat[rowi,posV] == mod_vec){
            alleles <- paste(alleles, rownames(seq_mat)[rowi], sep="   ",
                                                              collapse='');
         }
      }
   }
   return(alleles);
}

# Print out all modules at the given position vector posV,
# with the allele(s) containing that module at that positon set.
# The module counts (mc) matrix must be named; see get_module_counts()
# If the module counts matrix is not provided it will be found using
# loci+dataset combination.
which_alleles <- function(posV, mc=c(), seq_mat=c(), loci=c(), dataset=c()){
   if(is.null(mc)){ 
      mc <- get_module_counts(posV=posV, dataset=dataset, loci=loci);
   }
   for(module in rownames(mc)){
      # Get the allele names as a string.
      str <- which_has_module(module=module, posV=posV, seq_mat=seq_mat, loci=loci, dataset=dataset);
      print(sprintf("%s: %s", module, str));
   }
}


# Add 2 lists.
# Lists are hashes keyed on names.  Values are numbers to be added.
# Return hash.
add2Lists <- function(l1, l2){
   allNames <- union(names(l1), names(l2));
   merged <- list();
   for(thisName in allNames){
      if(thisName %in% names(l1)){
         l1value <- l1[[thisName]];
      }
      else{
         l1value <- 0;
      }
      if(thisName %in% names(l2)){
         l2value <- l2[[thisName]];
      }
      else{
         l2value <- 0;
      }
      merged[[thisName]] <- l1value + l2value;
   }
   return(merged);
}


alleleCountsUnequal <- function(alleleStr, alleleCounts, expected_aff, expected_con){
   alleleNames <- unlist(strsplit(alleleStr, " +",)); # split on whitespace
   allRows <- c();
   for(thisName in alleleNames){
      aff <- alleleCounts[thisName, 1];
      con <- alleleCounts[thisName, 2];
      p <- signif(pvalue(aff, con), 4);
      p_uneq <- signif(pvalue_unequal(aff, con, expected_aff, expected_con), 4);
      thisRow <- cbind(aff, con, p, p_uneq);
      allRows <- rbind(allRows, thisRow);
      #print(sprintf(%s\t%i\t%i\t%e\t%e), quote=FALSE);
   }
   rownames(allRows) <- alleleNames;
   colnames(allRows) <- c("Affected", "Control", "   p null ", sprintf("   p %i:%i", expected_aff, expected_con));
   allRows;
}


# fixClusterCounts
# The function cluster is counting individuals multiple times.
# Use the orignial function to do the clustering, then use the group names
# of the returned clusters to do a real count of the patients in each group.
# countMatrix is the matrix of counts returned by the cluster function.
# Row names of countMatrix contain the clustered modules. 
# pm is the dataMat of a patientMatrix.  
# Look through pm and count the number of patients
# that have any of the modules for each cluster.
# posV is the vector of positions.
# seqHash is the hash of sequences for all alleles, keyed on allele name.
fixClusterCounts <- function(countMatrix, apm, cpm, posV, seqHash){

   # Build a list of clusters where each cluster is represented by a vector
   # of module names from the white-space-split rownames of countMatrix.
   clusterNames <- rownames(countMatrix);
   clustersList <- list();
   for(groupNumber in 1:length(clusterNames)){
      clustersList[[groupNumber]] <-
                                       # split on whitespace
                           unlist(strsplit(clusterNames[groupNumber], " +",));
   }
   affectedCounts <- rep(0, length(clustersList));
   controlCounts <- rep(0, length(clustersList));

   # Look at each row of the each patientMatrix and count patients with a 
   # module in each cluster set.
   alleleNames <- names(seqHash);
   for(row_i in 1:(max(nrow(apm),nrow(cpm)))){
      # Get the modules for this affected patient.
      module1a <- c(); # Module 1 affected.
      module2a <- c(); # Module 2 affected.
      # Make sure we don't run off the end of the affected matrix.
      if(row_i <= nrow(apm)){
         # Module 1 affected.
         if( (!is.null(apm[row_i, 1])) && (apm[row_i, 1] %in% alleleNames) ){
            allele1a <- apm[row_i, 1];
            module1a <-paste(seqHash[[allele1a]][posV], collapse='');
         }
         # Module 2 affected.
         if( (!is.null(apm[row_i, 2])) && (apm[row_i, 2] %in% alleleNames) ){
            allele2a <- apm[row_i, 2];
            module2a <-paste(seqHash[[allele2a]][posV], collapse='');
         }
      }

      # Get the modules for this control patient.
      module1c <- c(); # Module 1 control.
      module2c <- c(); # Module 2 control.
      # Make sure we don't run off the end of the control matrix.
      if(row_i <= nrow(cpm)){
         if( (!is.null(cpm[row_i, 1])) && (cpm[row_i, 1] %in% alleleNames) ){
            allele1c <- cpm[row_i, 1];
            module1c <-paste(seqHash[[allele1c]][posV], collapse='');
         }
         if( (!is.null(cpm[row_i, 2])) && (cpm[row_i, 2] %in% alleleNames) ){
            allele2c <- cpm[row_i, 2];
            module2c <-paste(seqHash[[allele2c]][posV], collapse='');
         }
      }

      # Check each cluster to see if this patient should be counted here.
      # Patients may be counted in more than one cluster.
      for(groupNumber in 1:length(clustersList)){
         # If either module of this patient is in this cluster group,
         # add one to this cluster's patient count.
         if((!is.null(module1a) && 
               (module1a %in% clustersList[[groupNumber]])) ||
            (!is.null(module2a) && 
               (module2a %in% clustersList[[groupNumber]])) ){
            affectedCounts[groupNumber] <- affectedCounts[groupNumber] + 1;
         }
         # Do the same thing for the control patientMatrix.
         if((!is.null(module1c) && 
               (module1c %in% clustersList[[groupNumber]])) ||
            (!is.null(module2c) && 
               (module2c %in% clustersList[[groupNumber]])) ){
            controlCounts[groupNumber] <- controlCounts[groupNumber] + 1;
         }
      }
   } # for row_i in 1:max(nrow(apm),nrow(cpm))


   # We have counts. Now format like the return from the cluster function.
   counts <- cbind(affectedCounts, controlCounts);
   rownames(counts) <- clusterNames;
   return(counts);
}

# Logging.
# Open the file name for appending and add the supplied message.
# If firsttime=TRUE, overwrite existing file, otherwise append.
logToFile <- function(logfilename, logmessage, timestamp=TRUE, firsttime=FALSE, echo=FALSE){
   # Should we add a timestamp?
   if(timestamp){
      logmessage <- sprintf("%s  %s", Sys.time(), logmessage);
   }
   if(firsttime==TRUE){
      # Over write existing.
      logfile <- file(logfilename, open="wt"); 
   }else{
      # Open the file for appending, in text mode.
      logfile <- file(logfilename, open="at"); 
   }
   writeLines(logmessage, con=logfile);
   if(echo){
      print(logmessage, quote=FALSE);
   }
   close(logfile);
}
