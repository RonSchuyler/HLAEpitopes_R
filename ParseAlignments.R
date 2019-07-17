#
# ParseAlignments.R
#
#
# Parse the archive of HLA allele alignments (Alignments_Rel_3360.zip) to replace loading from AlleleImport2.txt.
#
# Original parser of AlleleImport.txt and AlleleImport2.txt padded with 0s,
# so we will do that here.
#
# See get_padded_seqs_from_alignments(affectedAlleles, controlAlleles, zipfile)
#     as a replacement for get_padded_seqs()
#
# zipfile from: 
#   curl -O ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Alignments_Rel_3360.zip
#   ln -s Alignments_Rel_3360.zip Alignments.zip

# filenames <- unzip(alignments_zip_file, list=TRUE);
# prot_files <- grep(pattern="_prot.txt$", x=filenames, value=TRUE);
# Skip this one:  "alignments/ClassI_prot.txt" it is redunant with A, B, C, and less complete (as of version 3360).


# parse_prot_line()
# Parse a line from the _prot.txt alignment file.
# Return a single-element list containing sequence data from this line, named with alleleName:
# ret[[alleleName]] "MILNKALLLGALA"
# OR NULL if aline is comment, blank. 
# offset is characters from the start of the line, or 0 to ignore.
parse_prot_line <- function(aline=" DQA1*01:01:01:01         MIL NKALLLGALA ", offset=0){
   if(length(grep("^#", aline)) > 0){
      # # comment line
      return(c());
   }
   if(length(grep("^Please", aline)) > 0){
      # Please see http://hla.alleles.org/terms.html for terms of use.
      return(c());
   }

   spl <- unlist(strsplit(aline, split=" "));
   alleleName <- spl[2];

   if(is.na(alleleName) || alleleName == ""){
      # empty line, or spaces only
      # also catches alignment position pointers:                   |                                |
      return(c());
   }
   if(alleleName == "Prot"){
      # Prot              -30                              1
      return(c());
   }

   if(offset == 0){
      sequence <- paste(spl[3:length(spl)], collapse="");
   } else{
      aline_vec <- unlist(strsplit(aline, split=""));
      sequence <- paste(aline_vec[offset:length(aline_vec)], collapse="");
      sequence <- gsub(pattern=" ", replacement="", x=sequence, fixed=TRUE);
   }
   ret <- list(sequence);
   names(ret) <- alleleName;
   return(ret);
}

checkLength <- function( allAllelesList, ref_allele_name){
   maxLen <- 0;
   for(key in names(allAllelesList)){
      value <- allAllelesList[[key]];
      print(sprintf("%s :: %i", key, nchar(value)), q=F);
      if(nchar(value) > maxLen){
         maxLen <- nchar(value);
      }
   }
}

# 1. 
# fillFromRef()
# (must be done prior to rmDots(), or alignment will be off)
# replace '-' with amino acid from reference allele
# Values for allAllelesList on input and output are character strings: "MILNK..."
# allAllelesList keys are allele names
fillFromRef <- function(allAllelesList, ref_allele_name){
   ref_seq <- allAllelesList[[ref_allele_name]];
   # maxLen <- nchar(ref_seq); unused
   ref_seq_chars <- unlist(strsplit(ref_seq, split=""));
   for(alleleName in names(allAllelesList)){
      if(alleleName == ref_allele_name){
         next;
      }
      aseq <- allAllelesList[[alleleName]];

      # Fill from ref: replace '-' with amino acid from reference allele.
      zz <- unlist(strsplit(aseq, split=""));
      wdash <- which(zz == "-");
      if(length(wdash) > 0){
         zz[wdash] <- ref_seq_chars[wdash];
         aseq <- paste(zz, collapse="");
      }
      allAllelesList[[alleleName]] <- aseq;
   }
   return(allAllelesList);
}

# 2.
# rmDotsX()
# . -> 0 (if collapseDots==FALSE)   OR  cut it out (if collapseDots==TRUE)
# * -> 0
# X -> 0
# Values for allAllelesList on input are character strings: "MILNK..."
# Values for allAllelesList on output are character vectors:  "M" "I" "L" "N" "K"
# allAllelesList keys are allele names
rmDotsX_and0pad <- function(allAllelesList, collapseDots=TRUE){
   maxLen <- 0;
   # 1st pass: may change string lengths, so just track longest
   for(alleleName in names(allAllelesList)){
      aseq <- allAllelesList[[alleleName]];
      zeroseq <- gsub(pattern="*", replacement="0", x=aseq, fixed=TRUE); # replace '*'
      zeroseq <- gsub(pattern="X", replacement="0", x=zeroseq, fixed=TRUE); # remove stop codon 'X'
      if(collapseDots){
         zeroseq <- gsub(pattern=".", replacement="", x=zeroseq, fixed=TRUE);
      }else{
         zeroseq <- gsub(pattern=".", replacement="0", x=zeroseq, fixed=TRUE);
      }
      if(nchar(zeroseq) > maxLen){
         maxLen <- nchar(zeroseq);
      }
      allAllelesList[[alleleName]] <- zeroseq;
   } # 1st pass

   # 2nd pass: 
   #    pad the ends with 0's
   #    expand character strings to vectors
   for(alleleName in names(allAllelesList)){
      aseq <- allAllelesList[[alleleName]];
      if(nchar(aseq) < maxLen){
         # pad the end
         aseq <- sprintf("%s%s", aseq, paste(rep("0", (maxLen-nchar(aseq))), collapse=""));
      }
      # Split string into individual characters
      allAllelesList[[alleleName]] <- unlist(strsplit(aseq, split=""));
   } # 2nd pass
   return(allAllelesList);
}


# get_padded_seqs_from_alignments()
# Get a list (hash) of aligned 0-padded sequences keyed by alleleName,
# given lists of allele names: affectedAlleles, controlAlleles (countHashes from get_allele_counts)
# Actual counts are not used, the hash names are just used as a list of allele names to get.
# zipfile from: 
#   curl -O ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Alignments_Rel_3360.zip
#   ln -s Alignments_Rel_3360.zip Alignments.zip
get_padded_seqs_from_alignments <- function(affectedAlleles, controlAlleles, collapseDots=TRUE, zipfile="../Alignments.zip"){
   allAllelesList <- list();
   alleleNames <- unique( c(names(affectedAlleles), names(controlAlleles)) );

   locusToGet <- unique(sub("\\*.*$", replacement="", x=alleleNames)); # should only be one, but allow for multiple?

   #protfile="alignments/DQA1_prot.txt"; zipfile="../Alignments_Rel_3360.zip"
   #protfile="alignments/DRB_prot.txt"; zipfile="../Alignments_Rel_3360.zip"
   # In the alignments files, DRB1,3,4,5 are all in DRB_prot.txt
   if(locusToGet %in% c("DRB1","DRB3","DRB4","DRB5")){
      protfile <- sprintf("alignments/%s_prot.txt", "DRB");
   } else{
      protfile <- sprintf("alignments/%s_prot.txt", locusToGet);
   }
   pf <- unz(description=zipfile, filename=protfile);
   alignment <- readLines(pf);
   close(pf);

   protline <- unlist(strsplit(alignment[8], split=""));
   offset1 <- length(protline);
   # make sure protline endes in " 1" for expected offsets
   if( !(protline[offset1]==1) || !(protline[offset1-1]==" ") ){
      stop(sprintf("bad protline:%s", alignment[8]));
   }

   #refline <- unlist(strsplit(alignment[10], split=""));
   #startPos_line <- length(protline);
   #startPos_line <- nchar(alignment[8]);
   #strsplit(refline, split="")
   #grep(pattern, x, ignore.case = FALSE, perl = FALSE, value = FALSE, fixed = FALSE, useBytes = FALSE, invert = FALSE)
  
   space_split_ref <- unlist(strsplit(alignment[10], split=" "));
   ref_allele_name <- space_split_ref[2];
   #ref_line_numbers <- grep(ref_allele_name, alignment, fixed=TRUE); # unused

   # First 9 lines are header information. 
   # Start scanning allele sequences at line 10, which is the reference allele.
   for(line_number in 10:length(alignment)){
      aline <- alignment[line_number];
      if( length(grep("Prot", aline, fixed=TRUE)) > 0){ 
         # if this is another Prot line reset offset.
         offset1 <- 0;
      }

      # Skip this line if it doesn't contain locusToGet and is not the referece allele.
      # eg skip DRB1,3,4 if locusToGet is DRB5
      # but allow DRB1*01:01:01 (allele from 1st line (line 10) for DRB_prot.txt. It is the reference for DRB3,4,5 also.
      if( length(grep(locusToGet, aline, fixed=TRUE)) == 0 &&
          length(grep(ref_allele_name, aline, fixed=TRUE)) == 0 ){
         next;
      }

      newProt <- parse_prot_line(aline, offset=offset1);
      if(is.null(newProt)){
         next;
      }
      newName <- names(newProt);
      newSeq <- newProt[[newName]];
      if( length(grep(locusToGet, newName, fixed=TRUE)) == 0 &&
         newName != ref_allele_name){ # allow DRB1*01:01:01 for DRB_prot.txt. It is the reference for DRB3,4,5 also.
         # Don't know what this is, but it doesn't contain the locus.
         next;
      }
      if(newName %in% names(allAllelesList)){
         # Add this new bit of sequence onto the end.
         allAllelesList[[newName]] <- sprintf("%s%s", allAllelesList[[newName]], newSeq);
      }else{
         # insert new sequence
         allAllelesList[[newName]] <- newSeq;
      }
   }
   allAllelesList <- fillFromRef(allAllelesList, ref_allele_name);
   allAllelesList <- rmDotsX_and0pad(allAllelesList, collapseDots=collapseDots);
   return(allAllelesList);
}




# just comments and notes on parsing one protein alignment file. 
dont_call_parse_prot <- function(protfile="alignments/DPA1_prot.txt", zipfile="../Alignments_Rel_3360.zip"){

   # for testing:
   #protfile="alignments/DPA1_prot.txt"; zipfile="../Alignments_Rel_3360.zip"
   protfile="alignments/DQA1_prot.txt"; zipfile="../Alignments_Rel_3360.zip"

   pf <- unz(description=zipfile, filename=protfile);
   #  filename: a filename within a zip file.
   #  description:
   #     'unz' reads (only) single files within zip files, in binary mode.
   #     The description is the full path to the zip file, with '.zip'
   #     extension if required.

   # For file-like connections on Windows,
   #     translation of line endings (between LF and CRLF) is done in text
   #     mode only (but text read operations on connections such as
   #     'readLines', 'scan' and 'source' work for any form of line ending).
   #pa <- read.table(pf);
   alignment <- readLines(pf);
   for(line_number in 1:length(alignment)){
      aline <- alignment[line_number];
      # skip header info
      #  [1] "# file: DQA1_prot.txt"
      # [2] "# date: 2019-04-17"
      # [3] "# version: IPD-IMGT/HLA 3.36.0"
      # [4] "# origin: http://hla.alleles.org/wmda/DQA1_prot.txt"
      # [5] "# repository: https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/alignments/DQA1_prot.txt"
      # [6] "# author: WHO, Steven G. E. Marsh (steven.marsh.ac.uk)"
      # blank line:
      # [7] "                   "
      # get offset from line:  " Prot              -40                                         1"
      # skip lines with only pipes '|'
      # strip whitespace
      # first data line is reference. Reference is repeated 
      # Prot              -40                                         1
      #                   |                                           |
      # DPA1*01:03:01:01           M RPEDRMFHIR AVILRALSLA FLLSLRGAGA IKADHVSTYA AFVQTHRPTG EFMFEFDEDE MFYVDLDKKE TVWHLEEFGQ AFSFEAQGGL
      # DPA1*01:03:01:02           - ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
      # ...
      # blank line
      # Prot              61
      #                   |
      # DPA1*01:03:01:01  ANIAILNNNL NTLIQRSNHT QATNDPPEVT VFPKEPVELG QPNTLICHID KFFPPVLNVT WLCNGELVTE GVAESLFLPR TDYSFHKFHY LTFVPSAEDF
      # DPA1*01:03:01:02  ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
      #
      # Prot              161
      #                   |
      # DPA1*01:03:01:01  YDCRVEHWGL DQPLLKHWEA QEPIQMPETT ETVLCALGLV LGLVGIIVGT VLIIKSLRSG HDPRAQGTL
      # DPA1*01:03:01:02  ---------- ---------- ---------- ---------- ---------- ---------- ---------

      # Example: DQA1_prot.txt  
      #     NOTE: X is stop codon, lines may be of variable length, with or w/o X, may end padded with *, or blank entirely.
      # 
      # Prot              -30                              1
      #                   |                                |
      # DQA1*01:01:01:01         MIL NKALLLGALA LTTVMSPCGG EDIVADHVAS CGVNLYQFYG PSGQYTHEFD GDEEFYVDLE RKETAWRWPE FSKFGGFDPQ GALRNMAVAK
      # DQA1*01:01:01:02         --- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ---------- ----------
      # DQA1*01:26               *** ********** ********** *****----- ---------- ---------- ---------- W--------- ---------- ----------
      # DQA1*02:01:01:01         --- ----M----- ---------- ---------- Y------S-- ----F----- ---------- ----V-KL-L -HRLRF.--- F--T-I--L-
      # DQA1*05:15N              *** ********** ********** *****----- Y------S-- ---------- ---Q-----G ----GVCLFS DNLD.LTRNL HX
      #
      # Prot              171
      #                   |
      # DQA1*01:01:01:01  WGLDQPLLKH WEPEIPAPMS ELTETVVCAL GLSVGLVGIV VGTVFIIQGL RSVGASRHQG PL
      # DQA1*05:15N
      # DQA1*06:02        ----E----- -********* ********** ********** ********** ********** **



      # Example: DPB1_prot.txt  NOTE: gap in reference
      #
      # Prot              -30                              1
      #                   |                                |
      # DPB1*01:01:01:01   MMVLQVSAA PRTVALTALL MVLLTSVVQG RATPENYVYQ GRQECYAF.. .....NGTQR FLERYIYNRE EYARFDSDVG EFRAVTELGR PAAEYWNSQK
      # DPB1*935:01Q       ********* ********** ********** *****----- L------LRQ ECYAF----- ---------- -FV------- ---------- -DED------ 
      # DPB1*936:01Q       ********* ********** ********** *****--LF- --------.. .....----- --D..----- -FV------- ---------- -DE------- 

   } # for each line in the alignment file
}







