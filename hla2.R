#
# hla2.R
# 

# Use abind from multiple dimensional matrices.
library("abind");


# Check if any of the elements of the vector argument are different.
# Return FALSE if all elements of vec are the same.
any_different <- function(vec){
   if(length(vec) < 2){
      return(FALSE);
   }
   for(index in 2:length(vec)){
      if(vec[index] != vec[1]){
         return(TRUE);
      }
   }
   return(FALSE);
}

# Return a vector of counts for each AA indexed as AAvec below.
count_residues <- function(vec, ignore_char="0"){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   aa_counts <- rep(0,length(AAvec));
   for(index in 1:length(vec)){
      if(vec[index] == ignore_char){
         next;
      }
      for(j_index in 1:length(AAvec)){
         if(vec[index] == AAvec[j_index]){
            aa_counts[j_index] <- aa_counts[j_index] + 1;
            break;
         }
      }
   }
   return(aa_counts);
}

# Return a list with 2 componets.
# Vectors of counts in each gold category, indexed as AAvec.
count_freqs_gold <- function(vec, gold0, gold4, ignore_char="0"){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   gold0count <- rep(0,length(AAvec));
   gold4count <- rep(0,length(AAvec));
   # Look at each position in vec.
   for(index in 1:length(vec)){
      if(vec[index] == ignore_char){
         next;
      }
      # Find the residue in AAvec that matches the current residue in vec.
      for(j_index in 1:length(AAvec)){
         if(vec[index] == AAvec[j_index]){
            #aa_counts[j_index] <- aa_counts[j_index] + 1;
            gold0count[j_index] <- gold0count[j_index] + gold0[index];
            gold4count[j_index] <- gold4count[j_index] + gold4[index];
            break;
         }
      }
   }
   return(list(gold0count, gold4count));
}


# Find any positions in data_matrix that are not all the same (by columns).
# Ignore ignore_char.
# lowerLim and upperLim are the bounds on the residue positions that
# we are interested in.  Residues beyond this bounds won't contact the
# peptide or the T cell receptor.
# lowerLim=8, upperLim=93 are good values for class II molecules,
# what about class I?
# If a limit is desired, both must be specified. 
get_polys <- function(data_matrix, ignore_char="0", silent=TRUE,
                  lowerLim=c(), upperLim=c()){
   polymorphic_pos <- c();
   for(col in 1:ncol(data_matrix)){
      non_zero <- data_matrix[,col] != ignore_char;
      if(any_different(data_matrix[non_zero,col])){
         if(silent == FALSE){
            print(col);
         }
         polymorphic_pos <- c(polymorphic_pos, col);
      }
   }
   if(!is.null(lowerLim) && !is.null(upperLim)){
      polymorphic_pos <- polymorphic_pos[polymorphic_pos>=lowerLim &
                           polymorphic_pos<=upperLim];
   }
   return(polymorphic_pos);
}


# Version of rbind that will pad vectors of different lengths
# with pad_char, rather than recycling values.
rbind_pad <- function(mat, row, pad_char="0", deparse.level=0){
   if(length(row) == 0){
      return(mat);
   }
   if(is.null(mat)){
      return(row);
   }
   if(is.matrix(row) && is.vector(mat)){
      # swap
      #return(rbind_pad(row,mat,pad_char));
      tmp <- mat;
      mat <- row;
      row <- tmp;
      if(ncol(mat) > length(row)){
         # pad row
         padded_row <- c(row,rep(pad_char,ncol(mat)-length(row)));
         return(rbind(padded_row,mat, deparse.level=deparse.level));
      }
      else if(ncol(mat) < length(row)){
         # pad mat
         padding <- matrix(pad_char,nrow=nrow(mat),ncol=(length(row)-ncol(mat)));
         mat <- cbind(mat,padding);
         return(rbind(row,mat,deparse.level=deparse.level));
      }
      # matrix and row have the same number of columns
      return(rbind(row,mat,deparse.level=deparse.level));
   }
   if(is.vector(mat) && is.vector(row)){
      # Two vectors.
      if(length(mat) > length(row)){
         # pad row
         padded_row <- c(row,rep(pad_char,length(mat)-length(row)));
         return(rbind(mat,padded_row,deparse.level=deparse.level));
      }
      else if(length(mat) < length(row)){
         # pad mat
         padded_mat <- c(mat,rep(pad_char,length(row)-length(mat)));
         return(rbind(padded_mat,row,deparse.level=deparse.level));
      }
      else{
         # equal length
         return(rbind(mat,row,deparse.level=deparse.level));
      }
   }
   if(is.matrix(mat) && is.matrix(row)){
      if(ncol(mat) > ncol(row)){
         # pad row matrix
         padding <- matrix(pad_char,nrow=nrow(row),ncol=(ncol(mat)-ncol(row)));
         row <- cbind(row,padding);
         return(rbind(mat,row,deparse.level=deparse.level));
      }
      else if(ncol(mat) < ncol(row)){
         # pad mat
         padding <- matrix(pad_char,nrow=nrow(mat),ncol=(ncol(row)-ncol(mat)));
         mat <- cbind(mat,padding);
         return(rbind(mat,row,deparse.level=deparse.level));
      }
      else{
         return(rbind(mat,row,deparse.level=deparse.level));
      }
   }
   # Only possiblity left: matrix and a vector
   if(ncol(mat) > length(row)){
      # pad row
      padded_row <- c(row,rep(pad_char,ncol(mat)-length(row)));
      return(rbind(mat,padded_row,deparse.level=deparse.level));
   }
   else if(ncol(mat) < length(row)){
      # pad mat
      padding <- matrix(pad_char,nrow=nrow(mat),ncol=(length(row)-ncol(mat)));
      mat <- cbind(mat,padding);
      return(rbind(mat,row,deparse.level=deparse.level));
   }
   # matrix and row have the same number of columns
   return(rbind(mat,row,deparse.level=deparse.level));
}

# which not
# Find and print position and content of matrix where column pos
# is not L1 or L2. L2 defaults to the pad_char.
wn <- function(m, pos, L1, L2="0",dm_list=c()){
   w <- intersect(which(m[,pos]!=L1),which(m[,pos]!=L2));
   if(length(w) == 0){
      print("none");
   }
   for(i in 1:length(w)){
      if(is.null(dm_list)){
         print(sprintf("%s at row %i",m[w[i],pos],w[i]));
      }
      else{
         allele <- names(dm_list)[w[i]];
         print(sprintf("%s at allele %s (row %i)",m[w[i],pos],allele,w[i]));
      }
   }
}

# Find positions which are L1 for the column pos.
which_are <- function(m, pos, L1, dm_list=c()){
   for(w in which(m[,pos] == L1)){
      if(is.null(dm_list)){
         print(sprintf("%s at row %i",m[w,pos], w));
      }
      else{
         allele <- names(dm_list)[w];
         print(sprintf("%s at allele %s (row %i)",m[w,pos], allele,w));
      }
   }
}


# Print out the residues in column pos that are not in the majority.
minority <- function(dm, pos, dm_list=c()){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   aa_count <- count_residues(dm[,pos]);
   majority_i <- which(aa_count == max(aa_count));
   majority_aa <- AAvec[majority_i];
   all_nz_indices <- which(aa_count > 0);
   if(length(all_nz_indices) > 1 && length(all_nz_indices) <= 6){
      for(nz_index in which(aa_count > 0)){
         if(all(AAvec[nz_index] != majority_aa)){
            #wn(dm, pos, majority_aa, dm_list=dm_list);
            which_are(dm, pos, AAvec[nz_index], dm_list=dm_list);
         }
      }
   }
   for(nz_index in which(aa_count > 0)){
      print(sprintf("%i %s",aa_count[nz_index],AAvec[nz_index]));
   }
}


# depricated
# This function is no longer used.  Still here for reference.
# Given a vector of polymorphic positions, create a matrix containg the 
# count of each AA at each polymorphic position for Gold 0 and Gold 4.
# Data is from get_data_matrix().
poly_freqs <- function(which_allele="DR", data_file="COPD_Cleaned_Up.txt"){

   if(which_allele == "DR"){
      allele_index <- 5;
      allele_file="DRB1.txt";
   }
   else if(which_allele == "DQ"){
      allele_index <- 7;
      allele_file="alignedDQ.txt";
   }
   else{
      print(sprintf("%s : unknown allele",which_allele));
      return(sprintf("%s : unknown allele",which_allele));
   }
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   gdm <- get_data_matrix(allele_file);
   data <- gdm[[1]]; # matrix of all aa
   data_allele_indexed <- gdm[[2]]; # list of all aa indexed by allele name
   poly_pos <- get_polys(data);
   gold0_list <- list(); 
   gold4_list <- list(); 
   for(pp in poly_pos){
      # initialize
      for(aa in AAvec){
         gold0_list[[as.character(pp)]][[aa]] <- 0;
         gold4_list[[as.character(pp)]][[aa]] <- 0;
      }
   }

   count <- 0;
   con <- file(data_file, open="r");
   header_line <- readLines(con,1);
   while(1){
      line <- readLines(con,1);
      if(length(line) == 0){
         break;
      }
      count <- count + 1;
      element <- unlist(strsplit(line,",",fixed=TRUE));
      element <- gsub("\"","",element,fixed=TRUE);
      gold <- as.numeric(element[2]);
      allele1 <- element[allele_index];
      allele2 <- element[(allele_index+1)];
      #print(sprintf("%f %s %s",gold,allele1,allele2));
      if(allele2 == "-"){
         allele2 <- allele1;
      }

      for(pp in poly_pos){
      #for(pp in c(3,9))
         ppc <- as.character(pp);
         if(gold == 0){
            if(pp <= length(data_allele_indexed[[allele1]])){
               letter <- data_allele_indexed[[allele1]][pp];
               if(letter != "0"){
                  gold0_list[[ppc]][[letter]] <- gold0_list[[ppc]][[letter]]+1;
               }
            }
            if(pp <= length(data_allele_indexed[[allele2]])){
               letter <- data_allele_indexed[[allele2]][pp];
               if(letter != "0"){
                  gold0_list[[ppc]][[letter]] <- gold0_list[[ppc]][[letter]]+1;
               }
            }
         }
         else if(gold == 4){
            if(pp <= length(data_allele_indexed[[allele1]])){
               letter <- data_allele_indexed[[allele1]][pp];
               if(letter != "0"){
                  gold4_list[[ppc]][[letter]] <- gold4_list[[ppc]][[letter]]+1;
               }
            }
            if(pp <= length(data_allele_indexed[[allele2]])){
               letter <- data_allele_indexed[[allele2]][pp];
               if(letter != "0"){
                  gold4_list[[ppc]][[letter]] <- gold4_list[[ppc]][[letter]]+1;
               }
            }
         }
         else{
            print(sprintf("skipping gold: %f",gold));
         }
      }

   }
   close(con);
   result <- matrix("0",nrow=length(poly_pos),ncol=length(AAvec));
   row_index <- 0;
   ratios <- c();
   count <- 0;
   pvals <- c();
   for(pp in poly_pos){
      row_index <- row_index + 1;
      ppc <- as.character(pp);
      col_index <- 0;
      for(aa in AAvec){
         col_index <- col_index + 1;
         g0 <- gold0_list[[ppc]][[aa]];
         g4 <- gold4_list[[ppc]][[aa]];
         tot <- g0 + g4;
         #print(sprintf("%i,%i: %i + %i : %i",row_index,col_index,g0,g4,tot));
         if(tot > 0){
            count <- count + 1;
            ratios <- c(ratios,g4/tot);
            pval <- pvalue(g4,g0);
            pvals <- c(pvals,pval);
            result[row_index,col_index] <- 
                  sprintf("%i %i %s: %.2f (%i:%i) p=%.3f",
            count,pp,aa,(g4/tot),g4,g0,pval);
            print(result[row_index,col_index]);
         }
      }
   }

   plot(pvals);
   #return(gold0_list);
   #return(result);
   return(pvals);
}


# P-value: probability of seeing numbers as or more extreme than these, given
# the null hypothesis that the positive and negative classes are equally likely.
# pos: number of postive examples
# neg: number of netative examples
# This assumes a two-sided test.
# Example: pvalue(9,1) = Pr(9,1)+Pr(10,0)+Pr(1,9)+Pr(0,10)
pvalue <- function(pos, neg){
   # Make it commutative.
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   trials=pos+neg;
   # Faster to precompute this, then multiply after the summation.
   one_trial <- .5 ^ trials; 
   total <- 0;
   for(i in pos:trials){
      #total <- total + choose(trials, i) * one_trial;
      total <- total + choose(trials, i);
   }
   total <- total * one_trial;
   if(!is.finite(total)){
      return(pvalue_log(pos, neg));
   }
   total <- min(1, (2 * total)); # make it two-sided.
   return(total);
}

# P-value
# Add logs instead of mulplying probabilities to prevent overflow.
pvalue_log <- function(pos, neg){
   # Make it commutative.
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   trials=pos+neg;
   one_trial <- trials * log(.5);
   total <- 0;
   for(i in pos:trials){
      x <- one_trial + lchoose(trials, i);
      total <- total + exp(x);
   }
   total <- min(1, (2 * total)); # make it two-sided.
   return(total);
}


#  Make pvalue that accounts for different number of affected and 
#  control individuals.  
#  pvalue unequal 2
#  affected = number of affected individuals with this module
#  control = number of control individuals with this module
#  n_affected = sample size, total number of affected individuals
#  n_control = sample size, total number of control individuals
pvalue_ue2 <- function(affected, control, n_affected, n_control){
   # Build the contingency table.
   c_table <- matrix(data=c(affected, control, 
         (n_affected-affected), (n_control-control)), 
         nrow=2, ncol=2, byrow=TRUE);
   # Decide which test to use.
   if(any(c_table < 10)){
      # Some values too small for chisq test, use fisher's exact.
      test.result <- fisher.test(c_table, 
         alternative = "two.sided", conf.int = FALSE);
   }
   else{
      # Use chi square test.
      test.result <- chisq.test(c_table, correct=FALSE);
   }
   return(test.result$p.value);
}

pvalue_unequal <- function(pos, neg, expected_pos, expected_neg){
  # Make it commutative.
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
      exp_tmp <- expected_pos;
      expected_pos <- expected_neg;
      expected_neg <- exp_tmp;
   }
   trials=pos+neg;
   p_pos <- expected_pos / (expected_pos + expected_neg);
   p_neg <- expected_neg / (expected_pos + expected_neg);

   total <- 0;
   for(i in pos:trials){
      total <- total + choose(trials, i) * (p_pos^i) * (p_neg^(trials-i));
   }
   if(!is.finite(total)){
      return(pvalue_unequal_log(pos, neg, expected_pos, expected_neg));
   }
   total <- min(1, (2 * total)); # make it two-sided.
   return(total);
}

pvalue_unequal_log <- function(pos, neg, expected_pos, expected_neg){
   # Make it commutative.
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
      exp_tmp <- expected_pos;
      expected_pos <- expected_neg;
      expected_neg <- exp_tmp;
   }
   trials=pos+neg;
   p_pos <- expected_pos / (expected_pos + expected_neg);
   p_neg <- expected_neg / (expected_pos + expected_neg);
   log_p_pos <- log(p_pos);
   log_p_neg <- log(p_neg);

   #one_trial <- trials * log(.5);
   total <- 0;
   for(i in pos:trials){
      x <- i*log_p_pos + (trials-i)*log_p_neg + lchoose(trials, i);
      total <- total + exp(x);
   }
   total <- min(1, (2 * total)); # make it two-sided.
   return(total);
}

# Another way of calculating the p-value.  
# Result is half the value of using pvalue().
# This one is a little slower.
# Note: it is also one-sided.
apv <- function(pos, neg){
   # Make it commutative.
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   trials=pos+neg;
   # Faster to precompute this, then multiply after the summation.
   #one_trial <- .5 ^ trials; 
   total <- 0;
   for(i in pos:trials){
      #total <- total + choose(trials, i) * one_trial;
      total <- total + choose(trials, i);
   }
   #total <- total * one_trial;
   less <- 0;
   for(i in 0:(pos-1)){
      less <- less + choose(trials, i);
   }
   #less <- less * one_trial;
   return(total/(total+less));
}


# P-value where we assume the probability of drawing a sample from this
# category (either pos or neg) is equal to the frequency of this category 
# in the  dataset.
pvf <- function(pos, neg, pf=c()){
   # pf: Pos Frequency or the prior probability of seeing a pos or neg example.
   if(is.null(pf)){
      pf <- (pos+neg)/199;
   }
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   trials=pos+neg;
   total <- 0;
   #one_trial <- .5 ^ trials;
   for(i in pos:trials){
      #total <- total + choose(trials, i) * (pf^i) * ((1-pf)^(trials-i));
      total <- total + choose(trials, i) * (pf^i) * (pf^(trials-i));
   }
   return(total);
}


# Compute the probability of seeing results as or more extreme than (pos,neg)
# from all alleles assuming the null hypothesis that the probabilities of each
# allele are the same in both categories.  Probability of an allele is the
# frequency of that allele in the data set. 
# This can be sped up if necessary: precompute (freq_allele_i^sum_posneg)
# One-sided.
pallele <- function(pos, neg, cat1, cat2){
   # cat1 and cat2 : vectors of counts in each category
   #        ie gold0 and gold4
   #           cat1[1] = count of allele 1 in category 1
   #           cat2[9] = count of allele 9 in category 2

   # Make it commutative.
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
      tmp <- cat1;
      cat1 <- cat2;
      cat2 <- tmp;
   }

   total <- 0; # result
   sum_posneg <- pos + neg;
   sum_cat <- sum(cat1,cat2); # scalar count of all instances
   both_cat <- cat1 + cat2; # vector sum of events in both categories
   # Calculate the probability of seeing results as or more extreme than
   # this from each allele.
   for(allele_i in 1:length(both_cat)){ # for every allele
      # Frequency of allele_i (both categories).
      freq_allele_i <- both_cat[allele_i] / sum_cat; 
      for(x_i in pos:sum_posneg){ # x_i is the index for x or more of cat1
         #total <- total + choose(sum_posneg,x_i) * (freq_allele_i^x_i) *
                          #((1-freq_allele_i)^(sum_posneg-x_i));
         total <- total +
            choose(sum_cat,pos) * choose((sum_cat-pos),neg) * 
               (freq_allele_i ^ sum_posneg) * 
                  ((1-freq_allele_i)^(sum_cat-sum_posneg));
      }
   }
   return(total);
}

# Measure of how extreme the given results are.
# Sums the probabilities of seeing these values or values that are more extreme
# and the probabilities of seeing less extreme values.
# Value returned is the ratio of P(asOrMore_ex)/(P(asOrMore_ex)+P(less_ex))
# In other words: ratio of sum of the probabilities of seeing this or something
# more extreme to the sum of the probabilities of all possible outcomes.
# Similar to a p-value, but provides more meaningful relative values.
# pos: # of positive examples from our category of interest
# neg: # of negative examples " "
# outof: total # of samples in all categories
# f_pn: probability of seeing a sample from this category among all others,
#       and should be left as the default: (pos+net)/outof
as_or_more_x <- function(pos, neg, outof, f_pn=c()){
   # outof is the total number of samples (ie 199 for DR)
   if(neg > pos){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   if(pos == 0){
      # If pos=0 and neg=0
      return(-Inf);
   }
   if(is.null(f_pn)){
      # Frequency of samples from this category in compared to all others.
      f_pn <- (pos+neg)/outof;
   }

   # Start with prob of exactly pos and neg of outof.
   #total <- choose(outof,pos)*choose((outof-pos),neg)*
               #(f_pn^(pos+neg)) * ((1-f_pn)^(outof-pos-neg));
   asOrMoreEx <- 0;
   lessEx <- 0;
   as_x <- neg/pos;
   #countAsOrMore <- 0;
   #countLess <- 0;
   # as or more extreme
   for(p in 1:outof){
      for(n in 0:outof){
         if((p+n) > outof){
            # Not possible.
            next;
         }
         if( (n/p) > as_x){
            # Less extreme.
            #countLess <- countLess + 1;
            lessEx <- lessEx + 
                     choose(outof,p) * choose((outof-p),n) *
                        (f_pn^(p+n)) * ((1-f_pn)^(outof-p-n));
         }
         else{
            # As or more extreme.
            #countAsOrMore <- countAsOrMore + 1;
            asOrMoreEx <- asOrMoreEx + 
                     choose(outof,p) * choose((outof-p),n) *
                        (f_pn^(p+n)) * ((1-f_pn)^(outof-p-n));
         }

      }
   }
   #print(sprintf("more:%f    less: %f     tot:%f   m/t: %f    m/l:%f",
      #asOrMoreEx, lessEx, (asOrMoreEx+lessEx), asOrMoreEx/(asOrMoreEx+lessEx), 
      #asOrMoreEx/lessEx));
   #return(c( (asOrMoreEx/(asOrMoreEx+lessEx)), asOrMoreEx));
   #return(c( (asOrMoreEx/(asOrMoreEx+lessEx)), 
    #           (countAsOrMore/(countAsOrMore+countLess))));
   # This is the one used:
   return(asOrMoreEx/(asOrMoreEx+lessEx));
   #return(c((asOrMoreEx/lessEx),(asOrMoreEx/(asOrMoreEx+lessEx))));
}

x_rep <- function(pos, neg, outof, reps=1000){
   if(pos < neg){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   both <- pos + neg;
   freq_either <- both/outof;
   as_or_more <- 0;
   for(i in 1:reps){
      ps <- 0;
      ns <- 0;
      for(j in 1:outof){
         if(runif(1) <= freq_either){
            if(runif(1) >= .5){
               ps <- ps + 1;
            }
            else{
               ns <- ns + 1;
            }
         }
      }
      if(ps >= pos && ns <= neg){
         as_or_more <- as_or_more + 1;
      }
   }
   print(sprintf("%i of %i as or more ex:  %f", as_or_more, reps, 
      as_or_more/reps));
}

# Another variation on the p-value, keeping the number of p+n the same,
# but including the probability of choosing in this category.
# Not really what is intended.
p3 <- function(pos, neg, outof){
   if(pos < neg){
      tmp <- pos;
      pos <- neg;
      neg <- tmp;
   }
   pn <- pos + neg;
   fr <- pn/outof;
   first <- (.5 ^ pn) * choose(outof,pn) * (fr^pn) * ((1-fr)^(outof-pn));
   tot <- 0;
   for(x in pos:pn){
      tot <- tot + choose(pn,x);
   }
   return(tot*first);
}


p2 <- function(pos, neg, loops=1000000){
   x <- c();
   for(i in 1:loops){
      r <- runif((pos+neg));
      x[i] <- sum(r > .5);
   }
   result <- sum(x >= pos)/loops;
   return(result);
}

check_p <- function(filename="DRB1.txt"){
   data <- get_data_matrix(filename);
   dm <- data[[1]];
   dm_list <- data[[2]];
   polys <- get_polys(dm);
   for(i in polys){
      x <- readline("next");
      minority(dm,i,dm_list);
      print(sprintf("from col %i",i));
   }
}


# Read the padded alleles file.
# The file is assumed to have no header.
# File format is csv: allele_name,count_in_gold0,count_in_gold4,padded_seq_data
read_copd <- function(whichFile="DR"){
   if(whichFile == "DR"){
      filename <- "../data/DRallelesInCopdPadded.csv";
   }
   else{
      filename <- "../data/DQallelesInCopdPadded.csv";
   }
   con <- file(filename, open="r");
   count <- 0;
   seq_matrix <- c();
   gold0 <- c();
   gold4 <- c();
   allele_names <- c();
   while(1){
      line <- readLines(con,1);
      if(length(line) == 0){
         break;
      }
      count <- count + 1;
      element <- unlist(strsplit(line,",",fixed=TRUE)); # seperate csv
      # element 4 is the padded sequence
      #allele_name <- gsub("\"","",element[2],fixed=TRUE); # drop quotes
      seq_data <- element[4];
      aas <- unlist(strsplit(seq_data,split=""));
      seq_matrix <- rbind_pad(seq_matrix, aas);
      allele_names <- c(allele_names, element[1]);
      g0 <- as.numeric(element[2]);
      gold0 <- c(gold0, g0);
      g4 <- as.numeric(element[3]);
      gold4 <- c(gold4, g4);
      pv <- pvalue(g0,g4);
      #print(sprintf("%s   pvalue: %f       odds %2i / %2i      %f",element[1],pv,g0,g4, g0/g4));
      #data <- paste(aas,sep="",collapse="");
      #print(data);
   }
   close(con);

   #return(seq_matrix);
   # seq_matrix: sequence data matrix, offset adjusted and padded
   # allele_names, gold0, gold4 are vectors whose positions correspond to
   # rows of seq_matrix
   return(list(allele_names, gold0, gold4, seq_matrix));
}


# Calculate frequencies and p values at the allele level.
allele_freqs <- function(whichFile="DR"){
   all_data <- read_copd(whichFile);
   allele_names <- all_data[[1]];
   gold0 <- all_data[[2]];
   gold4 <- all_data[[3]];
   total_sample_count <- sum(gold0, gold4);
   seq_mat <- all_data[[4]]; # sequence data matrix, offset adjusted and padded
   nrowSM <- nrow(seq_mat);
   poly_pos <- get_polys(seq_mat); 
   output_file_base <- "../data/alleles";
   output_file_name <- paste(output_file_base, whichFile, ".csv", sep='');
   output_file_handle <- file(output_file_name, open="w");
   for(row in 1:nrow(seq_mat)){
      g0 <- gold0[row];
      g4 <- gold4[row];
      px <- as_or_more_x(g0, g4, total_sample_count);
      pv <- pvalue(g0,g4);
      poly_pos <- get_polys(seq_mat); 
      this_string <- sprintf("%s,%s,%f,%f,%i,%i,%f",allele_names[row],
         paste(seq_mat[row,poly_pos],collapse=""),px,pv,g4,g0,g4/g0);
      #print(this_string);
      writeLines(this_string, con=output_file_handle);
   }
   close(output_file_handle);
}


MI <- function(abvec){
   aa <- abvec[1];
   ab <- abvec[2];
   ba <- abvec[3];
   bb <- abvec[4];
   a1 <- aa + ab;
   b1 <- ba + bb;
   a2 <- aa + ba;
   b2 <- ab + bb;
   xybase <- 2;
   total <- aa + ab + ba + bb;
   hx <- 0; hy <- 0; hxy <- 0;
   if(a1 > 0){
      hx <- hx - a1 * logb((a1/total), base=2);
   }
   if(a2 > 0){
      hy <- hy - a2 * logb((a2/total), base=2);
   }
   if(b1 > 0){
      hx <- hx - b1 * logb((b1/total), base=2);
   }
   if(b2 > 0){
      hy <- hy - b2 * logb((b2/total), base=2);
   }
   if(aa > 0){
      hxy <- hxy - aa * logb((aa/total), base=xybase);
   }
   if(ab > 0){
      hxy <- hxy - ab * logb((ab/total), base=xybase);
   }
   if(ba > 0){
      hxy <- hxy - ba * logb((ba/total), base=xybase);
   }
   if(bb > 0){
      hxy <- hxy - bb * logb((bb/total), base=xybase);
   }

   mi <- hx + hy - hxy;
   print(sprintf("MI: %f,    scaled: %f", mi, mi/total));
   #return(mi);
}

singles <- function(whichFile="DR"){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   lenAA <- length(AAvec);

   all_data <- read_copd(whichFile);
   allele_names <- all_data[[1]];
   gold0 <- all_data[[2]];
   gold4 <- all_data[[3]];
   total_sample_count <- sum(gold0, gold4);
   seq_mat <- all_data[[4]]; # sequence data matrix, offset adjusted and padded
   nrowSM <- nrow(seq_mat);
   poly_pos <- get_polys(seq_mat); 
   output_file_base <- "../data/singles_v2";
   output_file_name <- paste(output_file_base, whichFile, ".csv", sep='');
   output_file_handle <- file(output_file_name, open="w");
   for(posi in 1:length(poly_pos)){
      pos <- poly_pos[posi];
      G0countAA <- rep(0,lenAA);
      G4countAA <- rep(0,lenAA);
      for(seq_mat_row in 1:nrowSM){
         for(aaindex in 1:lenAA){
            if(seq_mat[seq_mat_row,pos] == AAvec[aaindex]){
               G0countAA[aaindex] <- G0countAA[aaindex] + gold0[seq_mat_row];
               G4countAA[aaindex] <- G4countAA[aaindex] + gold4[seq_mat_row];
            }
         }
      }
      for(aaindex in 1:lenAA){
         if((G0countAA[aaindex] > 0) ||
            (G4countAA[aaindex] > 0) ){
            g0 <- G0countAA[aaindex];
            g4 <- G4countAA[aaindex];
            pv <- pvalue(g0, g4);
            px <- as_or_more_x(g0,g4,total_sample_count);
            this_string <- sprintf("%i,%i,%s,%f,%f,%i,%i,%f", 
               posi, pos, AAvec[aaindex],px, pv, g4,g0,g4/g0);
            #print(this_string);
            writeLines(this_string, con=output_file_handle);
         }
      }
   }
   close(output_file_handle);
}




pairs <- function(whichFile="DR"){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   lenAA <- length(AAvec);

   all_data <- read_copd(whichFile);
   # read_copd returns list(allele_names, gold0, gold4, seq_matrix)
   # seq_matrix: sequence data matrix, offset adjusted and padded
   # allele_names, gold0, gold4 are vectors whose positions correspond to
   # rows of seq_matrix
   allele_names <- all_data[[1]];
   gold0 <- all_data[[2]];
   gold4 <- all_data[[3]];
   total_sample_count <- sum(gold0, gold4);
   seq_mat <- all_data[[4]]; # sequence data matrix, offset adjusted and padded
   nrowSM <- nrow(seq_mat);

   print("pas");
   for(index in 1:length(allele_names)){
      pa <- pallele(gold0[index],gold4[index],gold0,gold4);
      pv <- pvalue(gold0[index], gold4[index]);
      print(sprintf("%s  pa: %f   pvalue: %f   odds %2i / %2i  %f",
         allele_names[index], pa, pv,gold0[index],gold4[index], 
         gold0[index]/gold4[index]));
   }

   # Get the vector of polymorphic positions.  Note that this drops any
   # beyond our range of interest by default.  See get_polys().
   poly_pos <- get_polys(seq_mat); 


   # testing............................
   #poly_pos <- c(9,10,11,12);

   #px_vals <- c();
   #p_vals <- c();

   output_file_base <- "../data/pairs_v2";
   output_file_name <- paste(output_file_base, whichFile, ".csv", sep='');
   output_file_handle <- file(output_file_name, open="w");
   fourDimG0 <- list();
   fourDimG4 <- list();
   pos1_index <- 0;
   pos2_index <- 0;
   n_poly_pos <- length(poly_pos);
   for(pos1_index in 1:(n_poly_pos-1)){
      #pos1_index <- pos1_index + 1;
      pos1 <- poly_pos[pos1_index];
      all_pair_freqsG0 <- list();
      all_pair_freqsG4 <- list();
      for(pos2_index in (pos1_index+1):n_poly_pos){
         pos2 <- poly_pos[pos2_index];
         #print(sprintf("make matrix for %i:%i,  %i:%i", pos1_index, pos1, pos2_index, pos2));
         pair_freq_matG0 <- matrix(data=0,nrow=lenAA,ncol=lenAA);
         pair_freq_matG4 <- matrix(data=0,nrow=lenAA,ncol=lenAA);
         for(aarow in 1:lenAA){
            for(aacol in 1:length(AAvec)){
               for(sm_row in 1:nrowSM){ # seq mat index 1..nrow seq mat
                  if((seq_mat[sm_row,pos1] == AAvec[aarow]) &&
                     (seq_mat[sm_row,pos2] == AAvec[aacol])){
                        pair_freq_matG0[aarow,aacol] <- 
                           pair_freq_matG0[aarow,aacol] + gold0[sm_row];
                        pair_freq_matG4[aarow,aacol] <- 
                           pair_freq_matG4[aarow,aacol] + gold4[sm_row];
                  }
               }
            }
         }
         # check some values
         for(aarow in 1:lenAA){
            for(aacol in 1:length(AAvec)){
               if(pair_freq_matG0[aarow,aacol] > 0 ||
                  pair_freq_matG4[aarow,aacol] > 0){
                  g0 <- pair_freq_matG0[aarow,aacol];
                  g4 <- pair_freq_matG4[aarow,aacol];
                  pv <- pvalue(g0,g4);
                  px <- as_or_more_x(g0, g4, total_sample_count);
                  #px_vals <- c(px_vals, px);
                  #p_vals <- c(p_vals, pv);
                  this_string <- sprintf("%i,%i,%i,%i,%s%s,%f,%f,%i,%i,%f", 
                     pos1_index, pos2_index, pos1, pos2, AAvec[aarow], 
                     AAvec[aacol], px, pv, g4,g0,g4/g0);
                  #print(this_string);
                  writeLines(this_string, con=output_file_handle);

               }
            }
         }

         all_pair_freqsG0 <- c(all_pair_freqsG0, list(pair_freq_matG0));
         all_pair_freqsG4 <- c(all_pair_freqsG4, list(pair_freq_matG4));
      }
      fourDimG0 <- c(fourDimG0,list(all_pair_freqsG0));
      fourDimG4 <- c(fourDimG4,list(all_pair_freqsG4));

   }
   close(output_file_handle);
   #return(list(fourDimG0, fourDimG4));
   #return(list(p_vals,px_vals));
   # Each matrix is 20x20, with counts for residue 1=row, residue 2=col
   # fourDimG0 and G4 are 2-d lists
   # To access the matrix for counts of the second and 4th polymoric pos: 
   # fourDimG0[[2]][[2]];
   #   fourDimG0[[firstpos]][[secondpos-firstpos]]
   # keep in mind these are indexes, so the actual polypos is polypos[ppindex]
}


triples <- function(whichFile="DR"){
   AAstring <- "ACDEFGHIKLMNPQRSTVWY";
   AAvec <- unlist(strsplit(AAstring,split=""));
   lenAA <- length(AAvec);

   all_data <- read_copd(whichFile);
   # read_copd returns list(allele_names, gold0, gold4, seq_matrix)
   # seq_matrix: sequence data matrix, offset adjusted and padded
   # allele_names, gold0, gold4 are vectors whose positions correspond to
   # rows of seq_matrix
   allele_names <- all_data[[1]];
   gold0 <- all_data[[2]];
   gold4 <- all_data[[3]];
   total_sample_count <- sum(gold0,gold4);
   seq_mat <- all_data[[4]]; # sequence data matrix, offset adjusted and padded
   nrowSM <- nrow(seq_mat);

   # Get the vector of polymorphic positions.  Note that this drops any
   # beyond our range of interest by default.  See get_polys().
   poly_pos <- get_polys(seq_mat); 


   # testing............................
   #poly_pos <- c(9,10,11,12,13,14);
   #poly_pos <- c(9,10,11,12);


   #bestpx <- Inf;
   #beststr <- c();

   output_file_base <- "../data/triples_v2";
   output_file_name <- paste(output_file_base, whichFile, ".csv", sep='');
   output_file_handle <- file(output_file_name, open="w");

   #fourDimG0 <- list();
   #fourDimG4 <- list();
   pos1_index <- 0;
   pos2_index <- 0;
   pos3_index <- 0;
   n_poly_pos <- length(poly_pos);
   for(pos1_index in 1:(n_poly_pos-2)){
      #pos1_index <- pos1_index + 1;
      pos1 <- poly_pos[pos1_index];
      #all_pair_freqsG0 <- list();
      #all_pair_freqsG4 <- list();
      for(pos2_index in (pos1_index+1):(n_poly_pos-1)){
         pos2 <- poly_pos[pos2_index];
         #G02d <- matrix(data=0,nrow=lenAA,ncol=lenAA);
         #G42d <- matrix(data=0,nrow=lenAA,ncol=lenAA);
         for(pos3_index in (pos2_index+1):n_poly_pos){
            pos3 <- poly_pos[pos3_index];

            #print(sprintf("check %i,%i,%i : %i,%i,%i", 
             #  pos1_index, pos2_index, pos3_index, pos1, pos2, pos3));

           G0threeD <- c();
           G4threeD <- c();
           # [pos1,pos2,pos3]

            for(aadepth in 1:lenAA){
               G0threeD <- abind(G0threeD, matrix(0,nrow=lenAA,ncol=lenAA),
                  along=3);
               G4threeD <- abind(G4threeD, matrix(0,nrow=lenAA,ncol=lenAA),
                  along=3);
               for(aacol in 1:lenAA){
                  for(aarow in 1:lenAA){
                     for(sm_row in 1:nrowSM){ # seq mat index 1..nrow seq mat
                        if((seq_mat[sm_row,pos1] == AAvec[aarow]) &&
                           (seq_mat[sm_row,pos2] == AAvec[aacol]) &&
                           (seq_mat[sm_row,pos3] == AAvec[aadepth])){
                           G0threeD[aarow,aacol,aadepth] <-
                              G0threeD[aarow,aacol,aadepth] + gold0[sm_row];
                           G4threeD[aarow,aacol,aadepth] <-
                              G4threeD[aarow,aacol,aadepth] + gold4[sm_row];
                        }
                     }
                  }
               }
            }

            # Now check it.
            for(aarow in 1:lenAA){
               for(aacol in 1:lenAA){
                  for(aadepth in 1:lenAA){
                     if((G0threeD[aarow,aacol,aadepth] > 0) ||
                        (G4threeD[aarow,aacol,aadepth] > 0) ){
                        g0 <- G0threeD[aarow,aacol,aadepth]; 
                        g4 <- G4threeD[aarow,aacol,aadepth];
                        pv <- pvalue(g0,g4);
                        px <- as_or_more_x(g0,g4,total_sample_count);
                        #this_string <- sprintf("polyposi:%i,%i,%i : pos:%i,%i,%i : %s%s%s : px:%f, pv:%f, cg4/cg0: %i / %i   ratio: %f", pos1_index, pos2_index, pos3_index, pos1, pos2, pos3, AAvec[aarow], AAvec[aacol], AAvec[aadepth], px, pv, g4,g0,g4/g0);
                        this_string <- sprintf("%i,%i,%i,%i,%i,%i,%s%s%s,%f,%f,%i,%i,%f", pos1_index, pos2_index, pos3_index, pos1, pos2, pos3, AAvec[aarow], AAvec[aacol], AAvec[aadepth], px, pv, g4,g0,g4/g0);
                        #print(this_string);
                        writeLines(this_string, con=output_file_handle);
                        #if(px < bestpx){
                           #bestpx <- px;
                           #beststr <- this_string;
                        #}
                     }
                  }
               }
            }


         }
      }
   }
   close(output_file_handle);
   # output file format csv:
   # 3 polymorphic position indices, 3 positions, letters, px, pv, 
   # count of these 3 letter at these positions in gold 4, count in gold 0, 
   # ratio countg4/countg0 which = odds ratio when totals in each category are =
   #print(sprintf("best: %s",beststr));
}



# which file: 1=singles, 2=doubles, 3=triples
read_freqs <- function(whichFile=1, drdq="DR"){
   if(whichFile==1){
      if(drdq=="DR"){
         file_base <- "../data/singles";
      }
      else{
         file_base <- "../data/singles_v2";
      }
      # Indices for singles.csv:
      px_index <- 4;
      pv_index <- 5;
      g4_index <- 6;
      g0_index <- 7;
      ratio40 <- 8;
   }
   else if(whichFile==2){
      if(drdq=="DR"){
         file_base <- "../data/pairs";
      }
      else{
         file_base <- "../data/pairs_v2";
      }

      # Indices for pairs.csv:
      px_index <- 6;
      pv_index <- 7;
      g4_index <- 8;
      g0_index <- 9;
      ratio40 <- 10;
   }
   else if(whichFile==3){
      if(drdq=="DR"){
         file_base <- "../data/triples";
      }
      else{
         file_base <- "../data/triples_v2";
      }
      # Indices for threeD.csv:
      px_index <- 8;
      pv_index <- 9;
      g4_index <- 10;
      g0_index <- 11;
      ratio40 <- 12;
   }
   else{
      warning("unknown file type");
      return;
   }
   fn <- paste(file_base, drdq, ".csv", sep='');
   print(sprintf("Read file %s",fn));
   input <- file(fn, open="r");
  
   # Get the residue position distance matrix.
   dm <- get_residue_distance_matrix(drdq);

   px <- c();
   pv <- c();
   all_dists <- c();
   count <- 0;
   all_lines <- c();
   g0 <- c();
   g4 <- c();
   ratios <- c();
   while(1){
      line <- readLines(input,1);
      if(length(line) == 0){
         break;
      }
      count <- count + 1;
      element <- unlist(strsplit(line,",",fixed=TRUE));
      if(whichFile>=2){
         # Get residue distances.
         if(whichFile==2){
            # pairs
            pos1 <- as.numeric(element[3]);
            pos2 <- as.numeric(element[4]);
            d <- dm[pos1,pos2];
         }
         else{
            # triples
            pos1 <- as.numeric(element[4]);
            pos2 <- as.numeric(element[5]);
            pos3 <- as.numeric(element[6]);
            d <- max(dm[pos1,pos2], dm[pos1,pos3], dm[pos2,pos3]);
         }
         all_dists <- c(all_dists,d);
         line <- paste(line,",",d,sep='');
      }
      all_lines <- c(all_lines,line);
      #element <- gsub("\"","",element,fixed=TRUE);
      px <- c(px, as.numeric(element[px_index]));
      pv <- c(pv, as.numeric(element[pv_index]));
      g0 <- c(g0, as.numeric(element[g0_index]));
      g4 <- c(g4, as.numeric(element[g4_index]));
      ratios <- c(ratios, as.numeric(element[ratio40]));
   }
   
   close(input);
   print(sprintf("min px: %f", min(px)));
   print(sprintf("min pv: %f", min(pv)));
   ordered_i <- order(pv); # order by p-value
   #ordered_i <- order(ratios); # order by g4:g0 ratio
   ordered_lines <- all_lines[ordered_i];
   ret_list <- (list(px[ordered_i],pv[ordered_i],ratios[ordered_i],g0[ordered_i],g4[ordered_i],ordered_lines));
   if(whichFile>=2){
      ordered_dists <- all_dists[ordered_i];
      ret_list <- c(ret_list,
               list(ordered_dists));
   }
   nz_hist(ret_list);
   return(ret_list);
}

nz_hist <- function(xlist,b=20,silent=FALSE, quiet=TRUE){
   if(silent){
      quiet <- TRUE;
   }
   rat <- xlist[[3]];
   g0 <- xlist[[4]];
   g4 <- xlist[[5]];
   mg4 <- max(g4[g0==0]);
   mg0 <- max(g0[g4==0]);
   mig4 <- g4 == mg4 & g0==0;
   mig0 <- g0 == mg0 & g4==0;
   if(!quiet){
      print(xlist[[6]][mig4]);
      print(xlist[[6]][mig0]);
   }
   gt0 <- rat > 0;
   fin <- is.finite(rat);
   rat <- rat[gt0 & fin];
   maxratnz <- max(rat);
   minratnz <- min(rat);
   maxi <- xlist[[3]]==maxratnz;
   mini <- xlist[[3]]==minratnz;

   if(!quiet){
      print(xlist[[6]][maxi]);
      print(xlist[[6]][mini]);
   }
   if(!silent){
      print(xlist[[6]][xlist[[1]]==min(xlist[[1]])]);
      print(xlist[[6]][xlist[[2]]==min(xlist[[2]])]);
      print(sprintf("non-0 ratios   min:%f,   max:%f",min(rat), max(rat)));
   }
   hist(rat,breaks=b);
}


plot_junk <- function(x, xh=c(), xl=1){
   if(is.null(xh)){
      xh <- length(x[[1]]);
   }
   plotrange <- c(1,100);
   px <- x[[1]][xl:xh];
   pr <- x[[2]][xl:xh];
   pv <- x[[3]][xl:xh];
   px <- px/max(px);
   pr <- pr/max(pr);
   pv <- pv/max(pv);
   xmax <- length(px);
   ymax <- max(px,pv,pr);
   ymin <- min(px,pv,pr);

   t <- "b";
   plot(px, xlim=c(0,xmax), ylim=c(ymin,ymax),type=t);
   par(new=TRUE);
   plot(pr, xlim=c(0,xmax), ylim=c(ymin,ymax), col=2, pch=2,type=t);
   par(new=TRUE);
   plot(pv, xlim=c(0,xmax), ylim=c(ymin,ymax), col=4, pch=3,type=t);

}

# Read the PDB file and get distances between residues.
get_residue_distance_matrix <- function(drdq="DR"){
   fn <- paste("../data/",drdq,"pos.csv",sep='');
   posFile <- file(fn, open="r");
   count <- 0;
   max_residues <- 1000; 
   dist_mat <- matrix(0, nrow=max_residues, ncol=max_residues); 
   xs <- rep(0,max_residues);
   ys <- rep(0,max_residues);
   zs <- rep(0,max_residues);
   while(1){
      line <- readLines(posFile,1);
      if(length(line) == 0){
         break;
      }
      count <- count + 1;
      element <- unlist(strsplit(line,",",fixed=TRUE));
      #element <- gsub("\"","",element,fixed=TRUE);
      posi <- as.numeric(element[1]); # residue position in chain
      if(posi > max_residues){
         break;
      }

      xs[posi] <- as.numeric(element[3]); # x
      ys[posi] <- as.numeric(element[4]); # y
      zs[posi] <- as.numeric(element[5]); # z
   }
   close(posFile);

   # Resize.
   xs <- xs[1:count];
   ys <- ys[1:count];
   zs <- zs[1:count];

   # Make a symmetric distance matrix.
   dist_mat <- dist(cbind(xs,ys,zs),method="euclidean",diag=TRUE,upper=TRUE);

   return(as.matrix(dist_mat));
}
 
# False discovery rate for p values.
# This version assumes only positive dependence.
fdr <- function(pv, r=.1){
   pv <- sort(pv);
   k <- length(pv);
   for(i in k:1){
      #print(sprintf("i:%i, pv:%g",i,pv[i]));
      distar <- (i/k) * r;
      #print(sprintf("i:%i, p:%g, di*:%f",i,pv[i], distar));
      if(pv[i] <= distar){
         print(sprintf("accept %g at %i and better (lowest %i)",pv[i],i,i));
         return(TRUE);
      }
   }
   return(FALSE);
}

# FDR index
# This will return the indices from the argument pv of the accepted p values.
# This version assumes only positive dependence.
fdr_index <- function(pv, r=.1){
   sort_return <- sort(pv, index.return=TRUE);
   sorted_pv <- sort_return$x;
   k <- length(sorted_pv);
   for(i in k:1){
      #print(sprintf("i:%i, pv:%g",i,sorted_pv[i]));
      distar <- (i/k) * r;
      #print(sprintf("i:%i, p:%g, di*:%f",i,sorted_pv[i], distar));
      if(sorted_pv[i] <= distar){
         print(sprintf("accept %g at %i and better (lowest %i)",
               sorted_pv[i],i,i));
         return(sort_return$ix[1:i]);
      }
   }
   return(c());
}


# False discovery rate.
# This will return the indices from the argument pv of the accepted p values.
# This version makes no assumptions about dependence.
fdr_dep_index <- function(pv, r=.1){
   sort_return <- sort(pv, index.return=TRUE);
   P_values <- sort_return$x;
   k <- length(P_values);
   fdr_values <- rep(0, k);
   is_sig <- rep(0, k);
   # Calculate critical values 
   for(i in k:1){
      fdr_values[i] <- min(r, r * (k / ((k+1-i)^2)));
      if(P_values[i] <= fdr_values[i]){
         #is_sig[i:1] <- 1;
         #break;
         return(sort_return$ix[1:i]);
      }

   }
   return(c()); 
}

# False discovery rate.
# This version makes no assumptions about dependence.
# Returns a 3-column matrix [pvalues, d*[i], is_sigificant].
fdr_dep <- function(pv, r=.1){
   P_values <- sort(pv);
   k <- length(P_values);
   fdr_values <- c();
   is_sig <- rep(0, k);
   # Calculate critical values 
   for(i in 1:k){
      fdr_values[i] <- min(r, r * (k/((k+1-i)^2)));
      if(P_values[i] <= fdr_values[i]){
         is_sig[i] <- 1;
      }

   }
   res <- matrix(c(P_values, fdr_values, is_sig),byrow=FALSE, nrow=k, ncol=3);
   return(res); 
}



#####################

# Possibly useful FDR packages:
# p.adjust {stats}
# Threshold.FDR {AnalyzeFMRI}
# package `qvalue'











