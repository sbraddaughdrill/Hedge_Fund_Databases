# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Identifiers.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Combine Identifiers from Eurekahedge Data
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
cat("SECTION: IMPORT FUND DETAILS", "\n")
###############################################################################

Identifiers_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_Identifier_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)

# Identifiers_input0 <- data.frame(pull=NA,file_name=NA,read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
#                                       stringsAsFactors=FALSE)
# Identifiers_input0[,"file_name"] <- Identifiers_input0[,"file_clean"]
# Identifiers_input0[,"file_name"] <- gsub("\\\\","/",Identifiers_input0[,"file_name"])
# Identifiers_input0[,"file_name"] <- gsub("//","/",Identifiers_input0[,"file_name"])
# Identifiers_input0[,"file_name"] <- gsub("//","/",Identifiers_input0[,"file_name"])
# 
# Identifiers_input0[,"file_name"] <- encodeString(Identifiers_input0[,"file_name"])
# 
# #Identifiers_input0[,"pull"] <- regexpr("/[^/]*$", Identifiers_input0[,"file_name"])
# Identifiers_input0[,"pull"] <- sapply(gregexpr("\\/", Identifiers_input0[,"file_name"]), tail, 1)
# Identifiers_input0[,"file_name"] <- substr(Identifiers_input0[,"file_name"],Identifiers_input0[,"pull"]+1,nchar(Identifiers_input0[,"file_name"]))
# 
# Identifiers_input0[,"pull"] <- Identifiers_input0[,"file_name"]
# Identifiers_input0[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=Identifiers_input0[,"pull"])
# 
# Identifiers_input <- Identifiers_input0[grep("Identifiers",Identifiers_input0[,"file_name"]),]
# row.names(Identifiers_input) <- seq(nrow(Identifiers_input))
# rm2(Identifiers_input0)


# Identifiers_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
# Identifiers_yr <- c(2014,2013,2012,2011,2010,2009,2008,2007)
# 
# Identifiers_files <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Identifiers","EurekahedgeHF_EXCEL_2013Feb_Data_Identifiers",
#                             "EurekahedgeHF_EXCEL_2012Feb_Data_Identifiers","EurekahedgeHF_EXCEL_2011Feb_Data_Identifiers",
#                             "EurekahedgeHF_EXCEL_2010Feb_Data_Identifiers","EurekahedgeHF_EXCEL_2009Feb_Data_Identifiers",
#                             "EurekahedgeHF_EXCEL_2008Feb_Data_Identifiers","EurekahedgeHF_EXCEL_2007Feb_Data_Identifiers")
# 
# Identifiers_input <- data.frame(matrix(NA, ncol=3, nrow=length(Identifiers_files), dimnames=list(c(), c("pull","yr","file"))), 
#                                      stringsAsFactors=FALSE)
# 
# Identifiers_input[,"pull"] <- Identifiers_pull
# Identifiers_input[,"yr"] <- Identifiers_yr
# Identifiers_input[,"file"] <- Identifiers_files

Identifiers_concatenate0 <- alply(.data=Identifiers_input, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Identifiers_input[1,]
  # x <- Identifiers_input[5,]
  # x <- Identifiers_input[6,]
  # x <- Identifiers_input[7,]
  
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      #read.csv(file=paste(directory_in,"\\",x[,"yr"],"\\",x[,"file"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      stringsAsFactors=FALSE)
  
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
  
  input <- input[rowSums(is.na(input[,1:ncol(input)]))<ncol(input),]
  row.names(input) <- seq(nrow(input))
  
  input[,"pull"] <- x[,"pull"]
  
  input <- input[,colSums(is.na(input))<nrow(input)]
  
  input  <- input[order(input[,"Fund_ID"],input[,"Fund_Name"]),]
  row.names(input) <- seq(nrow(input))
  
  input_cols <- colnames(input)
  input_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=input_cols)
  input_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=input_cols)
  input_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=input_cols)
  input_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=input_cols)
  input_cols <- gsub(pattern="\\_{2,}", replacement="\\_", x=input_cols)
  input_cols <- gsub(pattern="\\_{2,}", replacement="\\_", x=input_cols)
  input_cols <- gsub(pattern="\\_{2,}", replacement="\\_", x=input_cols)
  input_cols <- gsub(pattern="\\_{2,}", replacement="\\_", x=input_cols)
  input_cols <- gsub(pattern="[[:punct:]]?$", replacement="", x=input_cols)
  
  
  #bad_ids  <- sort(unique(c(input[!is.na(input[,c("X")]),"Fund_ID"],input[!is.na(input[,c("X.1")]),"Fund_ID"])))
  #bad_rows  <- input[input[,c("Fund_ID")] %in% bad_ids,]
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  #Fix Date Added
  #   input[,"Date_Added"] <- gsub(pattern="/", replacement="-", x=input[,"Date_Added"])
  #   Date_Added_month <- gsub(pattern="([0-9]|[[:punct:]])", replacement="", x=input[,"Date_Added"])
  #   Date_Added_yr <- gsub(pattern="([^0-9])", replacement="", x=input[,"Date_Added"])
  #   Date_Added_yr <- as.integer(Date_Added_yr)
  #   Date_Added_yr_prefix <- Date_Added_yr
  #   Date_Added_yr_prefix <- ifelse(Date_Added_yr_prefix>=50,19,20)
  #   Date_Added_yr_suffix <- formatC(Date_Added_yr, width=2, format="d", flag="0") 
  #   
  #   input[,"Date_Added"] <- ifelse(is.na(input[,"Date_Added"]),NA,paste(Date_Added_month," ",Date_Added_yr_prefix,Date_Added_yr_suffix,sep=""))
  #   
  #   #input[,"Date_Added"] <- paste(Date_Added_month," ","20",Date_Added_yr,sep="")
  #   input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b %Y")
  #   input[,"Date_Added"] <- as.Date(input[,"Date_Added"],format="%m-%d-%Y")
  #   input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  #   
  #   rm(Date_Added_month,Date_Added_yr,Date_Added_yr_prefix,Date_Added_yr_suffix)
  
  input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
  #input[,"Date_Added"] <- as.Date(input[,"Date_Added"],format="%m-%d-%Y")
  input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  
  
  #Fix Dead Date
  #   input[,"Dead_Date"] <- gsub(pattern="/", replacement="-", x=input[,"Dead_Date"])
  #   Dead_Date_month <- gsub(pattern="([0-9]|[[:punct:]])", replacement="", x=input[,"Dead_Date"])
  #   Dead_Date_yr <- gsub(pattern="([^0-9])", replacement="", x=input[,"Dead_Date"])
  #   Dead_Date_yr <- as.integer(Dead_Date_yr)
  #   Dead_Date_yr_prefix <- Dead_Date_yr
  #   Dead_Date_yr_prefix <- ifelse(Dead_Date_yr_prefix>=50,19,20)
  #   Dead_Date_yr_suffix <- formatC(Dead_Date_yr, width=2, format="d", flag="0") 
  #   
  #   input[,"Dead_Date"] <- ifelse(is.na(input[,"Dead_Date"]),NA,paste(Dead_Date_month," ",Dead_Date_yr_prefix,Dead_Date_yr_suffix,sep=""))
  #   
  #   #input[,"Dead_Date"] <- paste(Dead_Date_month," ","20",Dead_Date_yr,sep="")
  #   input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b %Y")
  #   input[,"Dead_Date"] <- as.Date(input[,"Dead_Date"],format="%m-%d-%Y")
  #   input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  #   
  #   rm(Dead_Date_month,Dead_Date_yr,Dead_Date_yr_prefix,Dead_Date_yr_suffix)
  
  input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
  #input[,"Dead_Date"] <- as.Date(input[,"Dead_Date"],format="%m-%d-%Y")
  input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Identifiers_concatenate <- rbind.fill(Identifiers_concatenate0)

#Get colnames for all input files

Identifiers_common_colnames <- ldply(.data=Identifiers_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Identifiers_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                               reshape(Identifiers_common_colnames[,!(colnames(Identifiers_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                               Totals=NA,stringsAsFactors=FALSE)
Identifiers_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Identifiers_common_colnames),1) 
Identifiers_common_colnames_num_cols <- colnames(Identifiers_common_colnames)[!(colnames(Identifiers_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Identifiers_common_colnames[,"Totals"] <- rowMeans(Identifiers_common_colnames[,Identifiers_common_colnames_num_cols],na.rm=TRUE)
Identifiers_common_colnames <- Identifiers_common_colnames[order(Identifiers_common_colnames[,"Totals"],Identifiers_common_colnames[,"Order_All_Org"]),]
Identifiers_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Identifiers_common_colnames),1)
Identifiers_common_colnames[,"Totals"] <- rowSums(!is.na(Identifiers_common_colnames[,Identifiers_common_colnames_num_cols]))

Identifiers_good_cols <- Identifiers_common_colnames[Identifiers_common_colnames[,"Totals"]==length(Identifiers_common_colnames_num_cols),"cols"]
Identifiers_bad_cols <- Identifiers_common_colnames[Identifiers_common_colnames[,"Totals"]!=length(Identifiers_common_colnames_num_cols),"cols"]

#Identifiers_common_colnames <- rbind(Identifiers_common_colnames[Identifiers_common_colnames[,"Totals"]==length(Identifiers_common_colnames_num_cols),],
#                                 Identifiers_common_colnames[Identifiers_common_colnames[,"Totals"]!=length(Identifiers_common_colnames_num_cols),])
Identifiers_common_colnames <- rbind(Identifiers_common_colnames[Identifiers_common_colnames[,"cols"] %in% Identifiers_good_cols,],
                                          Identifiers_common_colnames[Identifiers_common_colnames[,"cols"] %in% Identifiers_bad_cols,])
Identifiers_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Identifiers_common_colnames),1)
colnames(Identifiers_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Identifiers_common_colnames))
row.names(Identifiers_common_colnames) <- seq(nrow(Identifiers_common_colnames))

#rm2(Identifiers_pull,Identifiers_yr,Identifiers_files)
rm2(Identifiers_concatenate0)

for(i in which(sapply(Identifiers_concatenate,class)=="character"))
{
  Identifiers_concatenate[[i]] = trim(Identifiers_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Identifiers_concatenate))
{
  Identifiers_concatenate[,i] <- unknownToNA(Identifiers_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Identifiers_concatenate[,i] <- ifelse(is.na(Identifiers_concatenate[,i]),NA, Identifiers_concatenate[,i])
} 
rm2(i)

Identifiers_concatenate  <- Identifiers_concatenate[order(Identifiers_concatenate[,"pull"],Identifiers_concatenate[,"Fund_ID"],Identifiers_concatenate[,"Fund_Name"]),]
row.names(Identifiers_concatenate) <- seq(nrow(Identifiers_concatenate))

rm2(Identifiers_input)


#Reorder Columns

# Identifiers_concatenate_all_cols <- colnames(Identifiers_concatenate)
# 
# Identifiers_concatenate_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
# Identifiers_concatenate_nonid_cols <- Identifiers_concatenate_all_cols[!(Identifiers_concatenate_all_cols %in% c(Identifiers_concatenate_id_cols))]
# 
# Identifiers_concatenate_primebroker_col <- Identifiers_concatenate_nonid_cols[grep("Prime.Broker", Identifiers_concatenate_nonid_cols)]
# Identifiers_concatenate_nonprimebroker_col <- Identifiers_concatenate_nonid_cols[!(Identifiers_concatenate_nonid_cols %in% c(Identifiers_concatenate_primebroker_col))]
# 
# Identifiers_concatenate_legaladvisor_col <- Identifiers_concatenate_nonprimebroker_col[grep("Legal.Advisor", Identifiers_concatenate_nonprimebroker_col)]
# Identifiers_concatenate_nonlegaladvisor_col <- Identifiers_concatenate_nonprimebroker_col[!(Identifiers_concatenate_nonprimebroker_col %in% c(Identifiers_concatenate_legaladvisor_col))]
# 
# 
# Identifiers_concatenate <- Identifiers_concatenate[,c(Identifiers_concatenate_id_cols,Identifiers_concatenate_nonlegaladvisor_col,
#                                                                 Identifiers_concatenate_primebroker_col,Identifiers_concatenate_legaladvisor_col)]
#
# rm2(Identifiers_concatenate_all_cols)
# rm2(Identifiers_concatenate_id_cols,Identifiers_concatenate_nonid_cols)
# rm2(Identifiers_concatenate_primebroker_col,Identifiers_concatenate_nonprimebroker_col)
# rm2(Identifiers_concatenate_legaladvisor_col,Identifiers_concatenate_nonlegaladvisor_col)


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)

write.csv(Identifiers_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Identifiers_cols",".csv",sep=""),row.names=FALSE)

rm2(Identifiers_common_colnames,common_col_folder_path)
rm2(Identifiers_common_colnames_num_cols,Identifiers_good_cols,Identifiers_bad_cols)


#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

write.csv(Identifiers_concatenate, file=paste(final_folder_path,"//","EurekahedgeHF_Identifiers",".csv",sep=""),row.names=FALSE)

rm2(Identifiers_concatenate,final_folder_path)
