# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Fund_Details.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Combine Fund Details from Eurekahedge Data
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

Details_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_Fund_Details_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)

# Details_input0 <- data.frame(pull=NA,file_name=NA,read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
#                            stringsAsFactors=FALSE)
# Details_input0[,"file_name"] <- Details_input0[,"file_clean"]
# Details_input0[,"file_name"] <- gsub("\\\\","/",Details_input0[,"file_name"])
# Details_input0[,"file_name"] <- gsub("//","/",Details_input0[,"file_name"])
# Details_input0[,"file_name"] <- gsub("//","/",Details_input0[,"file_name"])
# 
# Details_input0[,"file_name"] <- encodeString(Details_input0[,"file_name"])
# 
# #Details_input0[,"pull"] <- regexpr("/[^/]*$", Details_input0[,"file_name"])
# Details_input0[,"pull"] <- sapply(gregexpr("\\/", Details_input0[,"file_name"]), tail, 1)
# Details_input0[,"file_name"] <- substr(Details_input0[,"file_name"],Details_input0[,"pull"]+1,nchar(Details_input0[,"file_name"]))
# 
# Details_input0[,"pull"] <- Details_input0[,"file_name"]
# Details_input0[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=Details_input0[,"pull"])
# 
# Details_input <- Details_input0[grep("Fund_Details_Fund_Details",Details_input0[,"file_name"]),]
# row.names(Details_input) <- seq(nrow(Details_input))
# rm2(Details_input0)

#Details_yr <- unique(Details_input[,"yr"])

# Details_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
# Details_yr <- c(2014,2013,2012,2011,2010,2009,2008,2007)
# 
# Details_files <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Fund_Details","EurekahedgeHF_EXCEL_2013Feb_Data_Fund_Details",
#                    "EurekahedgeHF_EXCEL_2012Feb_Data_Fund_Details","EurekahedgeHF_EXCEL_2011Feb_Data_Fund_Details",
#                    "EurekahedgeHF_EXCEL_2010Feb_Data_Fund_Details","EurekahedgeHF_EXCEL_2009Feb_Data_Fund_Details",
#                    "EurekahedgeHF_EXCEL_2008Feb_Data_Fund_Details","EurekahedgeHF_EXCEL_2007Feb_Data_Fund_Details")
# 
# Details_input <- data.frame(matrix(NA, ncol=3, nrow=length(Details_files), dimnames=list(c(), c("pull","yr","file"))), 
#                             stringsAsFactors=FALSE)
# 
# Details_input[,"pull"] <- Details_pull
# Details_input[,"yr"] <- Details_yr
# Details_input[,"file"] <- Details_files


Details_concatenate0 <- alply(.data=Details_input, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Details_input[1,]
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
  
  input_cols <- gsub(pattern="Assets", replacement="Asset", x=input_cols)
  input_cols <- gsub(pattern="assets", replacement="assets", x=input_cols)
  
  input_cols <- gsub(pattern="Firm.s", replacement="Firms", x=input_cols)
  input_cols <- gsub(pattern="firm.s", replacement="firms", x=input_cols)
  input_cols <- gsub(pattern="Firm_s", replacement="Firms", x=input_cols)
  input_cols <- gsub(pattern="firm_s", replacement="firms", x=input_cols)
  
  input_cols <- gsub(pattern="Main_Investment.Strategy", replacement="Primary_Investment_Strategy_combcol", x=input_cols)
  input_cols <- gsub(pattern="Investment_Style", replacement="Primary_Investment_Strategy_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="Investment_Geography", replacement="Geography_combcol", x=input_cols)
  input_cols <- gsub(pattern="Geographical_Mandate", replacement="Geography_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="Currency", replacement="Currency_combcol", x=input_cols)
  input_cols <- gsub(pattern="Base_Currency", replacement="Currency_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="Accounting_Method_for_Performance_Fees", replacement="Accounting_Method_combcol", x=input_cols)
  input_cols <- gsub(pattern="Equalisation_Share_class", replacement="Accounting_Method_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="UCITS_Compliant", replacement="UCITS_combcol", x=input_cols)
  input_cols <- gsub(pattern="UCITS", replacement="UCITS_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  #Fix Date Added
input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
#input[,"Date_Added"] <- as.Date(input[,"Date_Added"],format="%m-%d-%Y")
input[,"Date_Added"] <- as.character(input[,"Date_Added"])


  #Fix Dead Date
input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
#input[,"Dead_Date"] <- as.Date(input[,"Dead_Date"],format="%m-%d-%Y")
input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])


  #Fix Inception Date
input[,"Inception_Date"] <- as.yearmon(input[,"Inception_Date"],format="%b-%y")
#input[,"Inception_Date"] <- as.Date(input[,"Inception_Date"],format="%m-%d-%Y")
input[,"Inception_Date"] <- as.character(input[,"Inception_Date"])

  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Details_concatenate <- rbind.fill(Details_concatenate0)

#Get colnames for all input files

Details_common_colnames <- ldply(.data=Details_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Details_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                      reshape(Details_common_colnames[,!(colnames(Details_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                      Totals=NA,stringsAsFactors=FALSE)
Details_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Details_common_colnames),1) 
Details_common_colnames_num_cols <- colnames(Details_common_colnames)[!(colnames(Details_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Details_common_colnames[,"Totals"] <- rowMeans(Details_common_colnames[,Details_common_colnames_num_cols],na.rm=TRUE)
Details_common_colnames <- Details_common_colnames[order(Details_common_colnames[,"Totals"],Details_common_colnames[,"Order_All_Org"]),]
Details_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Details_common_colnames),1)
Details_common_colnames[,"Totals"] <- rowSums(!is.na(Details_common_colnames[,Details_common_colnames_num_cols]))

Details_good_cols <- Details_common_colnames[Details_common_colnames[,"Totals"]==length(Details_common_colnames_num_cols),"cols"]
Details_bad_cols <- Details_common_colnames[Details_common_colnames[,"Totals"]!=length(Details_common_colnames_num_cols),"cols"]

#Details_common_colnames <- rbind(Details_common_colnames[Details_common_colnames[,"Totals"]==length(Details_common_colnames_num_cols),],
#                                 Details_common_colnames[Details_common_colnames[,"Totals"]!=length(Details_common_colnames_num_cols),])
Details_common_colnames <- rbind(Details_common_colnames[Details_common_colnames[,"cols"] %in% Details_good_cols,],
                                 Details_common_colnames[Details_common_colnames[,"cols"] %in% Details_bad_cols,])
Details_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Details_common_colnames),1)
colnames(Details_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Details_common_colnames))
row.names(Details_common_colnames) <- seq(nrow(Details_common_colnames))

#rm2(Details_pull,Details_yr,Details_files)
rm2(Details_concatenate0)

for(i in which(sapply(Details_concatenate,class)=="character"))
{
  Details_concatenate[[i]] = trim(Details_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Details_concatenate))
{
  Details_concatenate[,i] <- unknownToNA(Details_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Details_concatenate[,i] <- ifelse(is.na(Details_concatenate[,i]),NA, Details_concatenate[,i])
} 
rm2(i)

Details_concatenate  <- Details_concatenate[order(Details_concatenate[,"pull"],Details_concatenate[,"Fund_ID"],Details_concatenate[,"Fund_Name"]),]
row.names(Details_concatenate) <- seq(nrow(Details_concatenate))

rm2(Details_input)


#Reorder Columns

Details_concatenate_all_cols <- colnames(Details_concatenate)

Details_concatenate_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
Details_concatenate_nonid_cols <- Details_concatenate_all_cols[!(Details_concatenate_all_cols %in% c(Details_concatenate_id_cols))]

Details_concatenate_exposure_col <- Details_concatenate_nonid_cols[grep("Exposure", Details_concatenate_nonid_cols)]
Details_concatenate_nonexposure_col <- Details_concatenate_nonid_cols[!(Details_concatenate_nonid_cols %in% c(Details_concatenate_exposure_col))]

Details_concatenate_strategy_col <- Details_concatenate_nonexposure_col[grep("Strategy", Details_concatenate_nonexposure_col)]
Details_concatenate_nonstrategy_col <- Details_concatenate_nonexposure_col[!(Details_concatenate_nonexposure_col %in% c(Details_concatenate_strategy_col))]

#Details_concatenate_size_col <- Details_concatenate_nonstrategy_col[grep("Investment_Size", Details_concatenate_nonstrategy_col)]
Details_concatenate_size_col <- Details_concatenate_nonstrategy_col[grep("Size|Capacity|Asset", Details_concatenate_nonstrategy_col)]
Details_concatenate_nonsize_col <- Details_concatenate_nonstrategy_col[!(Details_concatenate_nonstrategy_col %in% c(Details_concatenate_size_col))]

Details_concatenate_currency_col <- Details_concatenate_nonsize_col[grep("Investment_Currency", Details_concatenate_nonsize_col)]
Details_concatenate_noncurrency_col <- Details_concatenate_nonsize_col[!(Details_concatenate_nonsize_col %in% c(Details_concatenate_currency_col))]

Details_concatenate_annualstats_col <- Details_concatenate_noncurrency_col[grep("Annualized", Details_concatenate_noncurrency_col)]
Details_concatenate_nonannualstats_col <- Details_concatenate_noncurrency_col[!(Details_concatenate_noncurrency_col %in% c(Details_concatenate_annualstats_col))]

Details_concatenate <- Details_concatenate[,c(Details_concatenate_id_cols,Details_concatenate_strategy_col,Details_concatenate_size_col,
                                              Details_concatenate_currency_col,Details_concatenate_exposure_col,Details_concatenate_annualstats_col,
                                              Details_concatenate_nonannualstats_col)]

rm2(Details_concatenate_all_cols)
rm2(Details_concatenate_id_cols,Details_concatenate_nonid_cols)
rm2(Details_concatenate_exposure_col,Details_concatenate_nonexposure_col)
rm2(Details_concatenate_strategy_col,Details_concatenate_nonstrategy_col)
rm2(Details_concatenate_size_col,Details_concatenate_nonsize_col)
rm2(Details_concatenate_currency_col,Details_concatenate_noncurrency_col)
rm2(Details_concatenate_annualstats_col,Details_concatenate_nonannualstats_col)


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)

write.csv(Details_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Fund_Details_cols",".csv",sep=""),row.names=FALSE)

rm2(Details_common_colnames,common_col_folder_path)
rm2(Details_common_colnames_num_cols,Details_good_cols,Details_bad_cols)


#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

write.csv(Details_concatenate, file=paste(final_folder_path,"//","EurekahedgeHF_Fund_Details",".csv",sep=""),row.names=FALSE)

rm2(Details_concatenate,final_folder_path)
