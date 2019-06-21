library(utils);
library(gWidgets);
#library(gWidgetsRGtk2);       # Comment this out for my local mac version.
#options(guiToolkit="RGtk2");  # Comment this out for my local mac version.
options(guiToolkit="tcltk"); # UNcomment this out for my local mac version.

EpitopesVersion <- "1.5";
lociOptions <- c("A","B","C","DPA1","DPB1","DQA1","DQB1","DRB1","DRB3","DRB4","DRB5");

source("main.R");

# Top level window
tlw <- gwindow(sprintf("HLA Epitopes v%s", EpitopesVersion), visible=TRUE);
#tlw <- gwindow(sprintf("HLA Epitopes v%s", EpitopesVersion), visible=TRUE, width=150, height=300);
group <- ggroup(horizontal=FALSE, container=tlw, expand=FALSE);

locusGroup <- ggroup(horizontal=FALSE, container=group, expand=FALSE);
locusLabel <- glabel("Locus", container=locusGroup, expand=FALSE);
#locusLabel <- glabel("Locus", container=locusGroup, expand=FALSE, font.attr=list(style="bold"))

addSpring(locusGroup);
locusList <- gdroplist(lociOptions, container=locusGroup, expand=FALSE);
svalue(locusList)<-"DRB1";
addSpace(locusGroup, 15, horizontal=FALSE);
#addSpace(group, 50);

positionsLabel <- glabel("Positions to Include: (Leave blank for all)", container=group);
positionsLabel2 <- glabel("or make patterns like this: 70:Q,71:!R,73:A", container=group);
positionsBox <- gedit("8:93", container=group);
addSpace(group, 10, horizontal=FALSE);
groupsLabel <- glabel("All groups of: (Leave blank for all)", container=group);
groupsBox <- gedit("2", container=group);

addSpace(group, 10, horizontal=FALSE);
fileLabel <- glabel("Input File (.csv) from data directory", container=group);
dataDirectoryContents <- dir("../data", pattern=".csv$");
fileList <- gdroplist(dataDirectoryContents, container=group);

addSpace(group, 10, horizontal=FALSE);
allelesWithModulesCheckbox <- gcheckbox("List alleles containing each module", container=group, checked=TRUE);
#toScreen <- gcheckbox("Print to screen", container=group, checked=FALSE);
clusterCheckbox <- gcheckbox("Cluster", container=group);

addSpace(group, 10, horizontal=FALSE);
sigLabel <- glabel("Significance threshold", container=group);
sigBox <- gedit("0.05", container=group);

#gstatusbar
#addSpace
#gpanedgroup
#glayout

addSpace(group, 15, horizontal=FALSE);
#sbw <- gstatusbar("sbw", container=tlw);
#status0 <- glabel("Ready0.", container=group);
gobutton <- gbutton("Go",container=group, handler = function(h,...){
#   svalue(status0) <- c("Running...",sprintf("Output directory: %s", normalizePath("../output")));
   svalue(status) <- c("Running...",sprintf("Output directory: %s", normalizePath("../output")));
   gofunction();
});

#gobutton <- gbutton("Go",container=group, handler=NULL);
#addhandlerclicked(gobutton, 
#   handler = function(h,...) {
#             svalue(h$obj) <-"click me again"
#             svalue(h$action) <- "Button has been clicked"
#             svalue(sbw) <- "new value sbw";
#             gofunction();
#             }, 
#   action = sbw);


status <- glabel("Ready.", container=group);
gofunction <- function(){
   #print(sprintf("size:%s",size(locusList)));
   #print(sprintf("size:%s",size(positionsBox)));
   #print(sprintf("groups of:%s",svalue(groupsBox)), quote=FALSE);
   #print(sprintf("listAlleles:%s",svalue(allelesWithModulesCheckbox)));
   #print(sprintf("cluster:%s",svalue(clusterCheckbox)));
   #print(sprintf("sig:%s",svalue(sigBox)));
   #print(sprintf("group size:%s", size(group)));

   svalue(status) <- c("Running");

   All <- "All";
   #print(sprintf("Locus:%s",svalue(locusList)),quote=FALSE);
   positionsTextbox <- svalue(positionsBox);
   #print(sprintf("Positions:%s",positionsTextbox),quote=FALSE);
   #print(sprintf("Positions:%s",svalue(positionsBox)));
   #print(sprintf("file:%s",svalue(fileList)));
   # Test for possible values of positionsTextbox that should mean "all"
   # "all" is a function, so test for that specifically
   # Other values to test for: All, null, Null, NULL
   if(is.null(positionsTextbox) || length(positionsTextbox)==0 ||
      is.function(positionsTextbox) || positionsTextbox == "All" || positionsTextbox == "all" ||
      positionsTextbox==""){
      positionsTextbox <- c();
      #print("Positions: all polymorphic");
   # Check for any uppercase letters for amino acids.
   }else if( length(grep("[[:alpha:]]+", positionsTextbox)) > 0){
      # Turn this into a 2 element list pattern.
      # This will tell the Analysis to do a pattern analysis.
      #cat("Using pattern:", positionsTextbox, "\n");
      positionsTextbox <- makePattern(positionsTextbox);
   }else{
      # Convert the string to an integer vector
      positionsTextbox <- eval(parse(text=positionsTextbox));
      #cat("Positions:", positionsTextbox, "\n");
   }
   dataSet <- sprintf("../data/%s", svalue(fileList));
   #cat("Data File:", dataSet,"\n");
   
   # If positionsTextbox is a list, then it is a pattern
   # so groupsOfN is ignored.
   groupsOfN <- svalue(groupsBox);
   if(!is.list(positionsTextbox)){
      # Test for possible values of groupsOfN that should mean "all"
      # "all" is a function, so test for that specifically
      # Other values to test for: All, null, Null, NULL
      if(is.null(groupsOfN) || length(groupsOfN)==0 ||
         is.function(groupsOfN) || groupsOfN == "All" || groupsOfN == "all" ||
         is.na(as.numeric(groupsOfN))){
         groupsOfN <- c();
         #cat("Groups of: Just this set.", "\n");
      }else{
         #groupsOfN <- as.numeric(groupsOfN);
         #cat("Groups of: ", groupsOfN, "\n");
      }
   }

   # Create the Analysis object.
   a <- Analysis(locus=svalue(locusList), dataFile=dataSet,
      positions=positionsTextbox, FDR=as.numeric(svalue(sigBox)), groupsOfN=groupsOfN,
      doCluster=svalue(clusterCheckbox),
      #outputToScreen=FALSE, 
      #outputToScreen=svalue(toScreen), 
      outputToScreen=FALSE,
      outputToFile=TRUE,
      #outputToFile=toFile,
      printAllelesWithModules=FALSE,
      #printAllelesWithModules=(svalue(toScreen) && svalue(allelesWithModulesCheckbox)),
      awmToFile=svalue(allelesWithModulesCheckbox),
      a.version=EpitopesVersion);
   #print(a@moduleSet@moduleFrame);

   #gmessage("Finished")
   svalue(status) <- a@uiStatusString;
}

#addSpace(group, 1, horizontal=FALSE);
#status <- glabel("", container=group, expand=FALSE);
#svalue(status) <- c("Running...",sprintf("Sending output to %s", normalizePath("../output")));

