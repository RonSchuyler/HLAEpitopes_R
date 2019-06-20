

alignments_zip_file <- "Alignments_Rel_3360.zip";
filenames <- unzip(alignments_zip_file, list=TRUE);

filenames

prot_files <- grep(pattern="_prot.txt$", x=filenames, value=TRUE);


# SKip this one:  "alignments/ClassI_prot.txt" it is redunant with A, B, C, and less complete (as of version 3360).

# Parse one protein alignment file. 
parse_prot <- function(protfile="alignments/DPA1_prot.txt", zipfile="Alignments_Rel_3360.zip"){

protfile="alignments/DPA1_prot.txt"; zipfile="Alignments_Rel_3360.zip"


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
   while(readLines(pf){
      # skip header info
      # blank line
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




   }
}






