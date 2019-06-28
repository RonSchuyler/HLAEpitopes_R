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
parse_prot_line <- function(aline=" DQA1*01:01:01:01         MIL NKALLLGALA "){
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

   sequence <- paste(spl[3:length(spl)], collapse="");
   ret <- list(sequence);
   names(ret) <- alleleName;
   return(ret);
}


# padAllAlleles_andFillFromRef()
# need to change:
# replace '-' with amino acid from reference allele
# * -> 0
# . -> 0  OR cut it out
# Check for truncated alleles:
# allAllelesList[["DQA1*05:15N"]]   "****************************-----Y------S---------------Q-----G----GVCLFSDNLD.LTRNLHX"
# X -> 0
# and pad the end with 0 so all have same length
padAllAlleles_andFillFromRef <- function(allAllelesList, ref_allele_name, collapseDots=TRUE){
   ref_seq <- allAllelesList[[ref_allele_name]];
   maxLen <- nchar(ref_seq);
   ref_seq_chars <- unlist(strsplit(ref_seq, split=""));
   for(alleleName in names(allAllelesList)){
      if(alleleName == ref_allele_name){
         next;
      }
      aseq <- allAllelesList[[alleleName]];
      zeroseq <- gsub(pattern="*", replacement="0", x=aseq, fixed=TRUE);
      if(collapseDots){
         zeroseq <- gsub(pattern=".", replacement="", x=zeroseq, fixed=TRUE);
      }else{
         zeroseq <- gsub(pattern=".", replacement="0", x=zeroseq, fixed=TRUE);
      }
      zeroseq <- gsub(pattern="X", replacement="0", x=zeroseq, fixed=TRUE); # remove stop codon 'X'
      if(nchar(zeroseq) < maxLen){
         # pad the end
         zeroseq <- sprintf("%s%s", zeroseq, paste(rep("0", (maxLen-nchar(zeroseq))), collapse=""));
      }

      # Fill from ref: replace '-' with amino acid from reference allele.
      zz <- unlist(strsplit(zeroseq, split=""));
      wdash <- which(zz == "-");
      if(length(wdash) > 0){
         zz[wdash] <- ref_seq_chars[wdash];
         zeroseq <- paste(zz, collapse="");
      }
       
      allAllelesList[[alleleName]] <- zeroseq;
   }
   return(allAllelesList);
}


# get_padded_seqs_from_alignments()
# Get a list (hash) of aligned 0-padded sequences keyed by alleleName,
# given lists of allele names: affectedAlleles, controlAlleles (countHashes from get_allele_counts)
# Actual counts are not used, the hash names are just used as a list of allele names to get.
# zipfile from: 
#   curl -O ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/Alignments_Rel_3360.zip
#   ln -s Alignments_Rel_3360.zip Alignments.zip
get_padded_seqs_from_alignments <- function(affectedAlleles, controlAlleles, zipfile="../Alignments.zip"){
   allAllelesList <- list();
   alleleNames <- unique( c(names(affectedAlleles), names(controlAlleles)) );

   locusToGet <- unique(sub("\\*.*$", replacement="", x=alleleNames)); # should only be one, but allow for multiple?

   #protfile="alignments/DQA1_prot.txt"; zipfile="../Alignments_Rel_3360.zip"
   protfile <- sprintf("alignments/%s_prot.txt", locusToGet);
   pf <- unz(description=zipfile, filename=protfile);
   alignment <- readLines(pf);

   protline <- unlist(strsplit(alignment[8], split=""))
   refline <- unlist(strsplit(alignment[10], split=""))
   startPos_line <- length(protline);
   startPos_line <- nchar(alignment[8]);
   #strsplit(refline, split="")
   #grep(pattern, x, ignore.case = FALSE, perl = FALSE, value = FALSE, fixed = FALSE, useBytes = FALSE, invert = FALSE)
  
   space_split_ref <- unlist(strsplit(alignment[10], split=" "));
   ref_allele_name <- space_split_ref[2];
   ref_line_numbers <- grep(ref_allele_name, alignment, fixed=TRUE);

   for(line_number in 1:length(alignment)){
      aline <- alignment[line_number];
      newProt <- parse_prot_line(aline);
      newName <- names(newProt);
      newSeq <- newProt[[newName]];
      if(is.null(newProt)){
         next;
      }
      if( length(grep(locusToGet, newName, fixed=TRUE)) == 0 ){
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
   allAllelesList <- padAllAlleles_andFillFromRef(allAllelesList, ref_allele_name);
   return(allAllelesList);
}




# just comments and notes on arsing one protein alignment file. 
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







