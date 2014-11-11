# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_NAV_AUM_melt.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Melt NAV & AUM Eurekahedge Data
#
###############################################################################

###############################################################################
cat("SECTION: INITIAL SETUP", "\n")
###############################################################################

# Clear workspace
rm(list = ls(all = TRUE))
rm(list = ls(all.names = TRUE))

# Limit History to not exceed 500 lines
Sys.setenv(R_HISTSIZE = 500)

repo <- c("http://cran.us.r-project.org")
options(repos = structure(repo))
options(install.packages.check.source = FALSE)

# String as factors is False -- used for read.csv
options(StringsAsFactors = FALSE)

# Default maxprint option
options(max.print = 500)
# options(max.print=99999)

# Memory limit
#memory.limit(size = 8183)

#Remove scientific notation if digits less than 100
options("scipen"=100)

#Uknown Strings
#unknowns_strings <- c("",".",NA,"na","n/a","n\a","NA","N/A","N\\A","<NA>","null","NULL",NULL,"nan","NaN",NaN,
#                      NA_integer_,"NA_integer_",NA_complex_,"NA_complex_",NA_character_,
#                      "NA_character_",NA_real_,"NA_real_")
unknowns_strings <- c(" ","\n","",".","n/a","na","NA",NA,"<NA>","null","NULL",NULL,"nan","NaN",NaN,Inf,
                      NA_integer_,"NA_integer_",NA_complex_,"NA_complex_",
                      NA_character_,"NA_character_",NA_real_,"NA_real_")

# Set location (1=HOME,2=WORK,3=LAPTOP,4=CORALSEA FROM HOME,5=CORALSEA FROM WORK,6=CORALSEA FROM LAPTOP)
Location <- 1

if (Location == 1) {
  
  #input_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\", mustWork=TRUE)
  input_directory <- normalizePath("F:/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\", mustWork=TRUE)
  output_directory <- normalizePath("F:/Import_Data/Data/Eurekahedge",winslash="\\", mustWork=TRUE)
  #function_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)    
  function_directory <- normalizePath("F:/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)  
  
} else if (Location == 2) {
  
  input_directory <- normalizePath("C:/Users/bdaughdr/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\", mustWork=TRUE)
  output_directory <- normalizePath("C:/Import_Data/Data/Eurekahedge",winslash="\\", mustWork=TRUE)
  function_directory <- normalizePath("C:/Users/bdaughdr/Dropbox/Research_Methods/R",winslash="\\", mustWork=TRUE)   
  
} else if (Location == 3) {
  
  input_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\", mustWork=TRUE)
  output_directory <- normalizePath("C:/Import_Data/Data/Eurekahedge",winslash="\\", mustWork=TRUE)
  function_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)
  
} else if (Location == 4) {
  
  input_directory <- normalizePath("H:/Research/Hedge_Fund_Databases/Data", winslash = "\\", mustWork = TRUE)
  #output_directory <- normalizePath("C:/Users/bdaughdr/Documents/Import_Data/Data/Eurekahedge",winslash="\\", mustWork=TRUE)
  output_directory <- normalizePath("H:/Research/Import_Data/Data/Eurekahedge", winslash = "\\", mustWork = TRUE)
  #function_directory <- normalizePath("//tsclient/C/Users/S.Brad/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)
  function_directory <- normalizePath("//tsclient/F/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)
  
} else if (Location == 5) {
  
  input_directory <- normalizePath("H:/Research/Hedge_Fund_Databases/Data", winslash = "\\", mustWork = TRUE)
  #output_directory <- normalizePath("C:/Users/bdaughdr/Documents/Import_Data/Data/Eurekahedge",winslash="\\", mustWork=TRUE)
  output_directory <- normalizePath("H:/Research/Import_Data/Data/Eurekahedge", winslash = "\\", mustWork = TRUE)
  function_directory <- normalizePath("//tsclient/C/Users/bdaughdr/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)
  
} else if (Location == 6) {
  
  input_directory <- normalizePath("H:/Research/Hedge_Fund_Databases/Data", winslash = "\\", mustWork = TRUE)
  #output_directory <- normalizePath("C:/Users/bdaughdr/Documents/Import_Data/Data/Eurekahedge",winslash="\\", mustWork=TRUE)
  output_directory <- normalizePath("H:/Research/Import_Data/Data/Eurekahedge", winslash = "\\", mustWork = TRUE)
  #function_directory <- normalizePath("//tsclient/C/Users/S.Brad/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)
  function_directory <- normalizePath("//tsclient/F/Dropbox/Research_Methods/R", winslash = "\\", mustWork = TRUE)
  
} else {
  
  cat("ERROR ASSIGNING DIRECTORIES", "\n")
  
}
rm(Location)


###############################################################################
cat("SECTION: FUNCTIONS", "\n")
###############################################################################

source(file=paste(function_directory,"functions_db.R",sep="\\"),echo=FALSE)
source(file=paste(function_directory,"functions_statistics.R",sep="\\"),echo=FALSE)
source(file=paste(function_directory,"functions_text_analysis.R",sep="\\"),echo=FALSE)
#source(file=paste(function_directory,"functions_text_parse.R",sep="\\"),echo=FALSE)
source(file=paste(function_directory,"functions_utilities.R",sep="\\"),echo=FALSE)


###############################################################################
# LIBRARIES;
cat("SECTION: LIBRARIES", "\n")
###############################################################################

#Load External Packages

external_packages <- c("compare","cwhmisc","data.table","fastmatch","foreign","formatR","gdata","gtools",
                       "Hmisc","koRpus","lubridate","mitools","pbapply","plyr","R.oo","reshape2","rJava","RWeka","RWekajars",
                       "sqldf","stringr","tcltk","tm","zoo")
invisible(unlist(sapply(external_packages,load_external_packages, repo_str=repo, simplify=FALSE, USE.NAMES=FALSE)))
installed_packages <- list_installed_packages(external_packages)

rm(external_packages,installed_packages,repo)


###############################################################################
cat("SECTION: IMPORT NAV & AUM", "\n")
###############################################################################

NAV_AUM_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
NAV_AUM_yr <- c("2014","2013","2012","2011","2010","2009","2008","2007")

NAV_AUM_files <- c("EurekahedgeHF_EXCEL_2014Apr_NAV_AUM","EurekahedgeHF_EXCEL_2013Feb_NAV_AUM",
                   "EurekahedgeHF_EXCEL_2012Feb_NAV_AUM","EurekahedgeHF_EXCEL_2011Feb_NAV_AUM",
                   "EurekahedgeHF_EXCEL_2010Feb_NAV_AUM","EurekahedgeHF_EXCEL_2009Feb_NAV_AUM",
                   "EurekahedgeHF_EXCEL_2008Feb_NAV_AUM","EurekahedgeHF_EXCEL_2007Feb_NAV_AUM")

NAV_AUM_input <- data.frame(matrix(NA, ncol=3, nrow=length(NAV_AUM_files), dimnames=list(c(), c("pull","yr","file"))), 
                            stringsAsFactors=FALSE)

NAV_AUM_input[,"pull"] <- NAV_AUM_pull
NAV_AUM_input[,"yr"] <- NAV_AUM_yr
NAV_AUM_input[,"file"] <- NAV_AUM_files

rm2(NAV_AUM_pull,NAV_AUM_yr,NAV_AUM_files)


###############################################################################
cat("SECTION: IMPORT NAV & AUM", "\n")
###############################################################################

NAV_AUM_id_cols <- c("Ret_AUM","Fund.ID","Fund.Name")

a_ply(.data=NAV_AUM_input, .margins=1, .fun = function(x,directory_in,unknowns,id_cols){
  
  # x <- NAV_AUM_input[1,]
  # x <- NAV_AUM_input[7,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  # id_cols <- NAV_AUM_id_cols
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(directory_in,"\\",x[,"yr"],"\\",x[,"file"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      stringsAsFactors=FALSE)

  input <- input[rowSums(is.na(input[,1:ncol(input)]))<ncol(input),]
  row.names(input) <- seq(nrow(input))
  
  input[,"pull"] <- x[,"pull"]
  
  #input <- input[,colSums(is.na(input))<nrow(input)]
  
  date_cols <- colnames(input)[!(colnames(input) %in% c("pull",id_cols))]
  date_cols <- gsub("\\."," ",date_cols)
  date_cols <- gsub("-"," ",date_cols)
  date_cols <- gsub("  "," ",date_cols)
  date_cols <- gsub("  "," ",date_cols)
  date_cols <- gsub("  "," ",date_cols)
  
  #Convert to common format
  date_cols <- as.yearmon(date_cols,format="%b %y")
  date_cols <- as.character(date_cols)
  
  colnames(input) <- c("pull",id_cols,date_cols)
  
  for(i in which(sapply(input,class)=="character"))
  {
    input[[i]] = trim(input[[i]])
  }
  rm(i)
  for (i in 1:ncol(input))
  {
    input[,i] <- unknownToNA(input[,i], unknown=unknowns,force=TRUE)
    input[,i] <- ifelse(is.na(input[,i]),NA,input[,i])
  } 
  rm(i)
  
  assign(x[,"file"], input, envir = .GlobalEnv)

  gc()
  
}, directory_in=output_directory, unknowns=unknowns_strings, id_cols=NAV_AUM_id_cols, .expand = TRUE, .progress = "text")


###############################################################################
cat("SECTION: MELT NAV & AUM", "\n")
###############################################################################

#Check to see if melt folder exists.  If not, create it.
melt_folder_path <- paste(output_directory, "NAV_AUM_melt", sep = "//", collapse = "//")  
create_directory(melt_folder_path,remove=1)

a_ply(.data=NAV_AUM_input, .margins=1, .fun = function(x,directory_out,unknowns,id_cols){
  
  # x <- NAV_AUM_input[1,]
  # directory_out <- output_directory
  # unknowns <- unknowns_strings
  # id_cols <- NAV_AUM_id_cols
  
  NAV_AUM_concatentate <- get(x[,"file"])
  
  rm(list=c(x[,"file"]), envir = .GlobalEnv)
  
  # Seperate Names, NAV and AUM
  
  #Names_concatentate <- NAV_AUM_concatentate[,c("pull","Fund.ID","Fund.Name")]
  #Names_concatentate_u <- unique(Names_concatentate)
  #row.names(Names_concatentate_u) <- seq(nrow(Names_concatentate_u))
  
  AUM_concatentate <- NAV_AUM_concatentate[NAV_AUM_concatentate[,"Ret_AUM"]=="AUM",!colnames(NAV_AUM_concatentate) %in% c("")]
  AUM_concatentate <- AUM_concatentate[order(AUM_concatentate[,"pull"],AUM_concatentate[,"Ret_AUM"],AUM_concatentate[,"Fund.ID"]),]
  row.names(AUM_concatentate) <- seq(nrow(AUM_concatentate))
  
  NAV_concatentate <- NAV_AUM_concatentate[NAV_AUM_concatentate[,"Ret_AUM"]=="Return",!colnames(NAV_AUM_concatentate) %in% c("")]
  NAV_concatentate <- NAV_concatentate[order(NAV_concatentate[,"pull"],NAV_concatentate[,"Ret_AUM"],NAV_concatentate[,"Fund.ID"]),]
  row.names(NAV_concatentate) <- seq(nrow(NAV_concatentate))
  
  rm(NAV_AUM_concatentate)
  
  
  # Melt AUM
  
  AUM_melt <- melt(AUM_concatentate[,!colnames(AUM_concatentate) %in% c("Ret_AUM")], id=c("pull","Fund.ID","Fund.Name"), na.rm=FALSE)
  rm(AUM_concatentate)
  
  colnames(AUM_melt)[match("variable",names(AUM_melt))] <- "date"
  colnames(AUM_melt)[match("value",names(AUM_melt))] <- "AUM"
  
  AUM_melt[,"AUM"] <- ifelse(is.na(AUM_melt[,"AUM"]), NA, AUM_melt[,"AUM"])
  
  AUM_melt[,"AUM"] <- gsub(",","", AUM_melt[,"AUM"])
  AUM_melt[,"AUM"] <- gsub(",","", AUM_melt[,"AUM"])
  AUM_melt[,"AUM"] <- gsub(",","", AUM_melt[,"AUM"])
  AUM_melt[,"AUM"] <- as.numeric(AUM_melt[,"AUM"])
  AUM_melt[,"AUM"] <- as.integer(AUM_melt[,"AUM"])
  
  # Melt NAV
  
  NAV_melt <- melt(NAV_concatentate[,!colnames(NAV_concatentate) %in% c("Ret_AUM")], id=c("pull","Fund.ID","Fund.Name"), na.rm=FALSE)

  rm(NAV_concatentate)
  
  colnames(NAV_melt)[match("variable",names(NAV_melt))] <- "date"
  colnames(NAV_melt)[match("value",names(NAV_melt))] <- "Monthly_Ret"
  
  NAV_melt[,"Monthly_Ret"] <- ifelse(is.na(NAV_melt[,"Monthly_Ret"]), NA, NAV_melt[,"Monthly_Ret"])

  NAV_melt[,"Monthly_Ret"] <- as.numeric(NAV_melt[,"Monthly_Ret"])
 
  # Merge
  
  NAV_AUM_merge <- merge(NAV_melt, AUM_melt, 
                         by.x=c("pull","Fund.ID","Fund.Name","date"), 
                         by.y=c("pull","Fund.ID","Fund.Name","date"), 
                         all.x=TRUE, all.y=TRUE, sort=FALSE,suffixes=c(".x",".y"))
  
  rm(NAV_melt,AUM_melt)
  
  #assign(paste(x[,"file"],"melt",sep="_"), NAV_AUM_merge, envir = .GlobalEnv)
  write.csv(NAV_AUM_merge,file=paste(directory_out,"\\",x[,"file"],".csv",sep=""),na="",quote=TRUE,row.names=FALSE)
  
  gc()
  
}, directory_out=melt_folder_path, unknowns=unknowns_strings, id_cols=NAV_AUM_id_cols, .expand = TRUE, .progress = "text")

rm2(NAV_AUM_input,NAV_AUM_id_cols)
