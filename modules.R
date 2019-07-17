#
# modules.R
#

# class modules
setClass("modules",
   representation(moduleFrame="data.frame")#,
   #prototype(affectedHash, controlHash, seqHash, positions)
);

# modules constructor
# 6/25/19 this function is not currently used.
modules <- function(
      affectedHash, controlHash, seqHash, positions, FDR, groupsOfN){
   stop("modules constructor is not used!");
   #print("modules constructor");
   startTime <- Sys.time();
   this <- new("modules");
   workingMat <- c();
   clusters <- c();
   # Get the set of all nchoosek combinations of positions.
   setList <- nofx(groupsOfN, positions);
   for(listI in 1:length(setList)){
      print(sprintf("combination %i, %s", listI, paste(setList[[listI]], collapse=' ')));

      #positionCounts <- get_position_counts(posV=positions, 
      positionCounts <- get_position_counts(posV=setList[[listI]], 
         controlCounts=controlHash, affectedCounts=affectedHash, 
            paddedSeqs=seqHash);
      # get_position_counts returns a list: 
      #     affected(hash), control(hash), pvaluesVec(vector)
      # see testSetsofN for an example of how to get all combinations
      matLen <- length(positionCounts[[3]]);
      isSignificant <- rep(0, matLen);
      sigIndices <- fdr_dep_index(positionCounts[[3]], r=FDR);
      isSignificant[sigIndices] <- 1;
      #workingMat <- rbind(workingMat, data.frame(cbind(
         #rep(paste(setList[[listI]], collapse=','), matLen),
         #two_lists_as_mat(positionCounts[[1]], positionCounts[[2]]), 
            #positionCounts[[3]], isSignificant), check.names=FALSE)); # todo add counts w alleleA

      countMat <- two_lists_as_mat(positionCounts[[1]], positionCounts[[2]]);
      print(sprintf("Cluster positions:  %s", paste(setList[[listI]], collapse=' ')));
      x<-cluster(countMat);
      thisCluster <- p_correct(x, printAccepted=FALSE, FDR=FDR);
      print(thisCluster);
      clusters <- rbind(clusters, thisCluster); 
      workingMat <- rbind(workingMat, cbind(rownames(countMat),
         rep(paste(setList[[listI]], collapse=','), matLen),
            countMat, 
               positionCounts[[3]], isSignificant)); # todo add counts w alleleA
   }
    
   colnames(workingMat) <- c("Module", "Positions", "Affected", "Control", "p-value", "Significant");
   #sort_order <- sort((as.numeric(workingMat$Affected) - 
                        #as.numeric(workingMat$Control)), 
                           #decreasing=TRUE, index.return=TRUE);
   sort_order <- sort((as.numeric(workingMat[,3]) - 
                        as.numeric(workingMat[,4])), 
                           decreasing=TRUE, index.return=TRUE);
   workingMat <- workingMat[sort_order$ix, ];
   matLen <- nrow(workingMat);
   isSignificant <- rep(0, matLen);
   sigIndices <- fdr_dep_index(as.numeric(workingMat[,5]), r=FDR);
   isSignificant[sigIndices] <- 1;
   #workingMat <- cbind(workingMat, isSignificant); 
   workingMat[, 6] <- isSignificant; 

   #this@moduleFrame <- data.frame(workingMat, check.names=FALSE);
   this@moduleFrame <- data.frame(workingMat);
   endTime <- Sys.time();
   print((endTime-startTime));
   print("All clusters:");
   print(clusters);

   this;

}


# printModules
printModules <- function(
      toScreen=TRUE, 
      toFile, # filepath for unclustered output or NULL
      doCluster=FALSE, # same as toScreen, but for clusters 
      clusterFile=c(), # filepath for cluster output or NULL 
      seqMat, # matrix of padded sequences
      seqHash, # precomputed hash of padded sequences, same as seqMat
      locus, # ex: DRB1
      affectedPatients, # patientMatrix 
      controlPatients,  # patientMatrix
      FDR=0.05, # false discovery rate 
      allPositions, # positions to look at
      groupsOfN,    # look at positions in groupsOfN, NULL means all at once
      printAllelesWithModules=FALSE, # list alleles containing each module
      awmToFile=FALSE){ # write AllelesWithModules to a file

   #print("printModules");
   if(is.null(allPositions)){
      stop("printModules got no positions. Check input file.", call.=FALSE);
   }
   startTime <- Sys.time();
   p_significantFigures <- 4; # number of signif. figures to save for p-values
   #this <- new("modules");
   # Initialization.
   allModules <- c(); 
   allClusters <- c();
   if(toScreen){
      # Make the screen print width wider if necessary.
      scrwid <- getOption("width");
      if(scrwid < 500){
         options(width=500);
      }
   }
   # Get the set of all nchoosek combinations of positions.
   setList <- nofx(groupsOfN, allPositions);
   loopCounter <- 0;
   for(posSet in setList){
      loopCounter <- loopCounter + 1;
      positionString <- paste(posSet, collapse=',');
      if(toScreen){
         print('', quote=FALSE);
         print(sprintf("Combination %i,  positions:  %s", loopCounter, 
                                    paste(posSet, collapse=' ')), quote=FALSE);
      }
      # Get counts at this set of positions for affected and control patients.
      Affected <- moduleCountHash(affectedPatients, locus, posV=posSet, seqMat);
      Control <- moduleCountHash(controlPatients, locus, posV=posSet, seqMat);
      # Combine the two hashes into a matrix.
      countMatrix <- two_lists_as_mat(Affected, Control);
      # Get p-values for individual modules and test for significance.
      countsWithPs <- p_correct(countMatrix, printAccepted=FALSE, FDR=FDR,
                                    significantFigures=p_significantFigures,
                                    n_affected=affectedPatients@patientCount,
                                    n_control=controlPatients@patientCount);
      allModules <- rbind(allModules, 
               cbind(rownames(countsWithPs), positionString, countsWithPs));
      #print(x);
      #countsWithPs <- data.frame(countsWithPs);
      #allModules <- rbind(allModules, countsWithPs);
      if(toScreen){
         print(countsWithPs);
         #countsWithPs <- data.frame(countsWithPs, check.names=FALSE, stringsAsFactors=FALSE);
      }
      if(doCluster){
         clusterSet <- cluster(countMatrix);
         clusterSet <- fixClusterCounts(clusterSet, 
                                 apm=affectedPatients@dataMat,
                                 cpm=controlPatients@dataMat,
                                 posV=posSet, seqHash=seqHash);

         clustersWithPs <- p_correct(clusterSet, printAccepted=FALSE, FDR=FDR,
                                       significantFigures=p_significantFigures,
                                    n_affected=affectedPatients@patientCount,
                                    n_control=controlPatients@patientCount);
         allClusters <- rbind(allClusters, #clustersWithPs);
               cbind(rownames(clustersWithPs), positionString, clustersWithPs));
         if(toScreen){
            print('', quote=FALSE);
            print(sprintf("Clusters  for positions:  %s", 
                                    paste(posSet, collapse=' ')), quote=FALSE);
            print(clustersWithPs);
         }
      }
   }
   # End of build up.

   # Multiple comparison test for significance for individual modules.
   allModules <- fdrAdjust(allModules, FDR=FDR,
                                       significantFigures=p_significantFigures,
                                       nAffected=affectedPatients@patientCount,
                                       nControl=controlPatients@patientCount);
   # Set column names.
   #if(is.null(dim(allModules))){
     #if(length(allModules) == 8){
         #names(allModules) <- c("Module", "Positions", "Affected", "Control", 
                                   #"p-value", "Accepted", "Accepted", "p-adj");
      #}
   #}
   #else{
      colnames(allModules) <- c("Module", "Positions", "Affected", "Control", 
                                   "p-value", "Accepted", "Accepted", "p-adj");
   #}

   # Do same for clusters if we have them.
   if(doCluster){
      allClusters <- fdrAdjust(allClusters, FDR=FDR,
                                       significantFigures=p_significantFigures,
                                       nAffected=affectedPatients@patientCount,
                                       nControl=controlPatients@patientCount);
      #if(is.null(dim(allClusters))){
         #if(length(allClusters) == 8){
            #names(allClusters) <-  c("Module", "Positions", "Affected", 
                        #"Control", "p-value", "Accepted", "Accepted", "p-adj");
         #}
      #}
      #else{
         colnames(allClusters) <- c("Module", "Positions", "Affected", 
                        "Control", "p-value", "Accepted", "Accepted", "p-adj");
      #}
   }

   # Write to file if necessary.
   if(!is.null(toFile)){
      outputFileHandle <- file(toFile, open="w");
      if(printAllelesWithModules || awmToFile){
         # Add alleles containing each module.
         print("write output, include allelesWithModules", q=F);
         #modulesWithAlleles <- cbind(allModules, ""); 
         modulesWithAlleles <- cbind(allModules, Alleles=""); 
         newcoli <- ncol(modulesWithAlleles);
         for(rowi in 1:nrow(allModules)){
            #an_allele <- which_has_module(module=allModules[rowi,1], posV=allModules[rowi,2], seq_mat=seqMat);
            #modulesWithAlleles[rowi, newcoli] <- an_allele;
            modulesWithAlleles[rowi, newcoli] <- which_has_module(module=allModules[rowi,1], posV=allModules[rowi,2], seq_mat=seqMat);
         }
         # Tab delimited table.
         write.table(modulesWithAlleles, file=outputFileHandle, quote=FALSE, row.names=FALSE, sep="\t");
      } else{
         print("write output w/o allelesWithModules", q=F);
         # Tab delimited table.
         write.table(allModules, file=outputFileHandle, quote=FALSE, row.names=FALSE, sep="\t");
      }
      close(outputFileHandle);
      if(doCluster && !is.null(clusterFile)){
         clusterFileHandle <- file(clusterFile, open="w");
         write.table(allClusters, file=clusterFileHandle, quote=FALSE, 
                                                   row.names=FALSE, sep="\t");
         close(clusterFileHandle);
         #allelesInFile(fileName=clusterFile, seqMat=seqMat, alpha=FDR);
      }

      # Create alleles file(s).
      #allelesInFile(fileName=toFile, seqMat=seqMat, alpha=FDR);
   }

   # Print everything to screen, if necessary.
   if(toScreen){
      print('', quote=FALSE);
      print("All modules:", quote=FALSE);
      print(data.frame(allModules, check.names=FALSE, stringsAsFactors=FALSE,
                        row.names=1:nrow(allModules)));
      #print(allClusters, quote=FALSE, row.names=FALSE);
      #write.table(allModules, file="", quote=FALSE, row.names=FALSE, sep="\t");
      if(doCluster){
         print('', quote=FALSE);
         print("All clusters:", quote=FALSE);
         allClusters <- data.frame(allClusters, check.names=FALSE,
                     stringsAsFactors=FALSE, row.names=1:nrow(allClusters));
         #write.table(allClusters, file="", quote=FALSE, row.names=FALSE, sep="\t");
         print(allClusters, quote=FALSE, row.names=FALSE);
      }
   }

   # Print each module followed by the list of alleles having it.
   if(printAllelesWithModules || awmToFile){
      if(printAllelesWithModules){
         print('', quote=FALSE);
         print("Alleles with the given module:", quote=FALSE);
      }
      # Open a connection to write alleles with module to a file.
      #if(awmToFile){
      #   awmFileCon <- openAWMfile(toFile); # awm written to main output file
      #}
      for(rowi in 1:nrow(allModules)){
         aLine <- sprintf("%s %s : %s", allModules[rowi,2], allModules[rowi,1], 
           which_has_module(module=allModules[rowi,1], 
           posV=allModules[rowi,2], seq_mat=seqMat));
         if(printAllelesWithModules){
            print(aLine, quote=FALSE);
         }
         #if(awmToFile){
         #   writeLines(aLine, con=awmFileCon); # awm written to main output file
         #}
      }
      #if(awmToFile){
      #   close(awmFileCon); # awm written to main output file
      #}
   }

   endTime <- Sys.time();
   print('', quote=FALSE);
   #print("printModules DONE", quote=FALSE);
   #timediff <- print((endTime-startTime));
   return((endTime-startTime));
   #return(list(allModules, allClusters));
   #return(allModules);
}

fdrAndSort <- function(dataMat, FDR=0.05, doSort=TRUE){
   # expected columns in dataMat
   modNameCol <- 1; # module name (character) ex: "QRAA"
   positionStrCol <- 2; # position string (character) ex: "70,71,72,73"
   affectedCol <- 3; # count of affected 
   controlCol <- 4; # count of control
   pCol <- 5; # p-values
   accCol <- 6; # Accepted from previous FDR
   newAccCol <- 7; # This is what we fill.
   if(doSort){
      sort_order <- sort((as.numeric(dataMat[,affectedCol]) - 
                          as.numeric(dataMat[,controlCol])), 
                                          decreasing=TRUE, index.return=TRUE);
      dataMat <- dataMat[sort_order$ix,];
   }
   dataMat <- cbind(dataMat, rep(0,nrow(dataMat)));
   fdri <- fdr_dep_index(as.numeric(dataMat[,pCol]), r=FDR);
   dataMat[fdri, newAccCol] <- 1;
   return(dataMat);
}

# Augment the given matrix of counts and pvalues with FDR adjusted pvalues.
fdrAdjust <- function(dataMat, FDR=0.05, doSort=TRUE, significantFigures=4,
               nAffected=0, nControl=0){
   # dataMat should have 6 columns when passed in. Return value has 8 columns.
   # expected columns in dataMat
   modNameCol <- 1; # module name (character) ex: "QRAA"
   positionStrCol <- 2; # position string (character) ex: "70,71,72,73"
   affectedCol <- 3; # count of affected 
   controlCol <- 4; # count of control
   pCol <- 5; # p-values
   accCol <- 6; # Accepted from previous FDR
   newAccCol <- 7; # 1 if this row passed cutoff after multiple comp. adjust. 
   pAdj <- 8; # Fill this with adjusted p-values.
   # This part is modified from the Benjamini & Yekutieli method of p.adjust.
   ps <- as.numeric(dataMat[,pCol]);
   n <- length(ps);
   i <- n:1;
   o <- order(ps, decreasing=TRUE);
   ro <- order(o);
   q <- sum(1/i);
   adjustedPs <- signif(pmin(1, cummin(q * n/i * ps[o]))[ro], 
                                                         significantFigures);
   dataMat <- cbind(dataMat, rep(0, n), adjustedPs);
   dataMat[(adjustedPs<=FDR), newAccCol] <- 1;
   # Break the matrix into 2 halves based on wether
   # affected or control counts is greater, then sort each half independently
   # on p-values, and recombine the halves.  
   if(doSort){
      # Get the rows that have more affected than control (susceptible).
      if(nAffected == 0 || nControl == 0){
         # This is the old way, which doesn't account for differences in 
         # affected vs control population sizes.
         suscep_i <- as.numeric(dataMat[,affectedCol]) > 
                                             as.numeric(dataMat[,controlCol]);
      }
      else{
         # Use population sizes to determine susceptibility/resistance.
         affRates <- as.numeric(dataMat[,affectedCol]) / nAffected;
         conRates <- as.numeric(dataMat[,controlCol]) / nControl;
         suscep_i <- affRates > conRates;
      }

      suscepMat <- dataMat[suscep_i,];
      resistMat <- dataMat[!suscep_i,];

      if(is.null(dim(suscepMat)) || is.null(dim(resistMat))){
         # Fall back to sorting just on affected-control difference.
         sort_order <- sort((as.numeric(dataMat[,affectedCol]) - 
                          as.numeric(dataMat[,controlCol])), 
                                          decreasing=TRUE, index.return=TRUE);
         dataMat <- dataMat[sort_order$ix,];
      }
      else{
         # Sort this group.
         suscep_o <- sort(as.numeric(suscepMat[, pAdj]), decreasing=FALSE, 
                                                   index.return=TRUE);
         # Sort the resistant group.
         resist_o <- sort(as.numeric(resistMat[, pAdj]), decreasing=TRUE, 
                                                   index.return=TRUE);
         dataMat <- rbind(suscepMat[suscep_o$ix,], resistMat[resist_o$ix,]);
      }

   }
   # Make sure what we return is a matrix.
   if(is.null(dim(dataMat))){
      dataMat <- matrix(dataMat, nrow=1);
   }
   return(dataMat);
}
  
# Like fdrAdjust, but takes output from countPatientsWithAlleles.
fdra <- function(dataMat, FDR=0.05, doSort=TRUE, significantFigures=3){
   # dataMat should have 2 columns when passed in. Return value has 6 columns.
   # expected columns in dataMat
   affectedCol <- 1; # count of affected 
   controlCol <- 2; # count of control
   pCol <- 3; # p-values
   accCol <- 4; # Accepted pre-FDR
   pAdj <- 5; # Fill this with adjusted p-values.
   newAccCol <- 6; # 1 if this row passed cutoff after multiple comp. adjust. 
   # This part is modified from the Benjamini & Yekutieli method of p.adjust.
   ps <- as.numeric(dataMat[,pCol]);
   n <- length(ps);
   i <- n:1;
   o <- order(ps, decreasing=TRUE);
   ro <- order(o);
   q <- sum(1/i);
   adjustedPs <- signif(pmin(1, cummin(q * n/i * ps[o]))[ro], 
                                                         significantFigures);
   dataMat <- cbind(dataMat, adjustedPs, rep(0, n));
   dataMat[(adjustedPs<=FDR), newAccCol] <- 1;
   # Break the matrix into 2 halves based on wether
   # affected or control counts is greater, then sort each half independently
   # on p-values, and recombine the halves.  
   if(doSort){
      # Get the rows that have more affected than control (susceptible).
      suscep_i <- as.numeric(dataMat[,affectedCol]) > 
                                          as.numeric(dataMat[,controlCol]);
      suscepMat <- dataMat[suscep_i,];
      resistMat <- dataMat[!suscep_i,];

      if(is.null(dim(suscepMat)) || is.null(dim(resistMat))){
         # Fall back to sorting just on affected-control difference.
         sort_order <- sort((as.numeric(dataMat[,affectedCol]) - 
                          as.numeric(dataMat[,controlCol])), 
                                          decreasing=TRUE, index.return=TRUE);
         dataMat <- dataMat[sort_order$ix,];
      }
      else{
         # Sort this group.
         suscep_o <- sort(as.numeric(suscepMat[, pAdj]), decreasing=FALSE, 
                                                   index.return=TRUE);
         # Sort the resistant group.
         resist_o <- sort(as.numeric(resistMat[, pAdj]), decreasing=TRUE, 
                                                   index.return=TRUE);
         dataMat <- rbind(suscepMat[suscep_o$ix,], resistMat[resist_o$ix,]);
      }

   }
   colnames(dataMat) <- c("Affected", "Control", "       p  ", "accepted", "   FDRadj p", " accepted");
   return(dataMat);
}



