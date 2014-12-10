# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Stats.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Combine Stats from Eurekahedge Data
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
cat("SECTION: IMPORT STATS", "\n")
###############################################################################

Stats_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_Stats_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)

# Stats_input0 <- data.frame(pull=NA,file_name=NA,read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
#                            stringsAsFactors=FALSE)
# Stats_input0[,"file_name"] <- Stats_input0[,"file_clean"]
# Stats_input0[,"file_name"] <- gsub("\\\\","/",Stats_input0[,"file_name"])
# Stats_input0[,"file_name"] <- gsub("//","/",Stats_input0[,"file_name"])
# Stats_input0[,"file_name"] <- gsub("//","/",Stats_input0[,"file_name"])
# 
# Stats_input0[,"file_name"] <- encodeString(Stats_input0[,"file_name"])
# 
# #Stats_input0[,"pull"] <- regexpr("/[^/]*$", Stats_input0[,"file_name"])
# Stats_input0[,"pull"] <- sapply(gregexpr("\\/", Stats_input0[,"file_name"]), tail, 1)
# Stats_input0[,"file_name"] <- substr(Stats_input0[,"file_name"],Stats_input0[,"pull"]+1,nchar(Stats_input0[,"file_name"]))
# 
# Stats_input0[,"pull"] <- Stats_input0[,"file_name"]
# Stats_input0[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=Stats_input0[,"pull"])
# 
# Stats_input <- Stats_input0[grep("Statistics",Stats_input0[,"file_name"]),]
# row.names(Stats_input) <- seq(nrow(Stats_input))
# rm2(Stats_input0)

Stats_yr <- unique(Stats_input[,"yr"])

# Stats_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
# Stats_yr <- c(2014,2013,2012,2011,2010,2009,2008,2007)
# 
# Stats_files <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Stats","EurekahedgeHF_EXCEL_2013Feb_Data_Stats",
#                  "EurekahedgeHF_EXCEL_2012Feb_Data_Stats","EurekahedgeHF_EXCEL_2011Feb_Data_Stats",
#                  "EurekahedgeHF_EXCEL_2010Feb_Data_Stats","EurekahedgeHF_EXCEL_2009Feb_Data_Stats",
#                  "EurekahedgeHF_EXCEL_2008Feb_Data_Stats","EurekahedgeHF_EXCEL_2007Feb_Data_Stats")
# 
# Stats_input <- data.frame(matrix(NA, ncol=3, nrow=length(Stats_files), dimnames=list(c(), c("pull","yr","file"))), 
#                           stringsAsFactors=FALSE)
# 
# Stats_input[,"pull"] <- Stats_pull
# Stats_input[,"yr"] <- Stats_yr
# Stats_input[,"file"] <- Stats_files

Stats_concatenate0 <- alply(.data=Stats_input, .margins=1, .fun = function(x,directory_in,unknowns,years){
  
  # x <- Stats_input[1,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  # years <- Stats_yr
  
  years_expand <- seq(min(years)-2,max(years),by=1)
  years_trim <- substr(years_expand, 3, 4)
  
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
  
  input_cols <- gsub(pattern="Returns", replacement="Return", x=input_cols)
  input_cols <- gsub(pattern="returns", replacement="return", x=input_cols)
  
  #Common Names -  Yearly Returns
  for (i in 1:length(years_expand))
  {
    # i <- 1
    
    input_cols <- gsub(pattern=paste("X",years_expand[i],"_","Return",sep=""), replacement=paste("Return","_",years_expand[i],"_yearlyret",sep=""), x=input_cols)
  } 
  rm(i)
  
  #Common Names -  Monthly Returns
  for (i in 1:length(years_expand))
  {
    # i <- 1
    
    input_cols <- gsub(pattern=paste("Jan","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","01","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Feb","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","02","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Mar","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","03","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Apr","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","04","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("May","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","05","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Jun","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","06","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Jul","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","07","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Aug","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","08","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Sep","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","09","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Oct","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","10","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Nov","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","11","_monthlyret",sep=""), x=input_cols)
    input_cols <- gsub(pattern=paste("Dec","_",years_trim[i],"_","Return",sep=""), 
                       replacement=paste("Return","_",years_expand[i],"_","12","_monthlyret",sep=""), x=input_cols)
    
  } 
  rm(i)
  
  rm(years_expand,years_trim)
  
  #Common Names -  VaR
  input_cols <- gsub(pattern="VaR_90", replacement="VaR_90pct", x=input_cols)
  input_cols <- gsub(pattern="VaR_95", replacement="VaR_95pct", x=input_cols)
  input_cols <- gsub(pattern="VaR_99", replacement="VaR_99pct", x=input_cols)
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

#rm(Date_Added_month,Date_Added_yr,Date_Added_yr_prefix,Date_Added_yr_suffix)

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
 
#rm(Dead_Date_month,Dead_Date_yr,Dead_Date_yr_prefix,Dead_Date_yr_suffix)

   input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
   #input[,"Dead_Date"] <- as.Date(input[,"Dead_Date"],format="%m-%d-%Y")
   input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])

  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, years=Stats_yr, .expand = TRUE, .progress = "text")

Stats_concatenate <- rbind.fill(Stats_concatenate0)

#Get colnames for all input files

Stats_common_colnames <- ldply(.data=Stats_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Stats_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                    reshape(Stats_common_colnames[,!(colnames(Stats_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                    Totals=NA,stringsAsFactors=FALSE)
Stats_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Stats_common_colnames),1) 
Stats_common_colnames_num_cols <- colnames(Stats_common_colnames)[!(colnames(Stats_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Stats_common_colnames[,"Totals"] <- rowMeans(Stats_common_colnames[,Stats_common_colnames_num_cols],na.rm=TRUE)
Stats_common_colnames <- Stats_common_colnames[order(Stats_common_colnames[,"Totals"],Stats_common_colnames[,"Order_All_Org"]),]
Stats_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Stats_common_colnames),1)
Stats_common_colnames[,"Totals"] <- rowSums(!is.na(Stats_common_colnames[,Stats_common_colnames_num_cols]))

Stats_good_cols <- Stats_common_colnames[Stats_common_colnames[,"Totals"]==length(Stats_common_colnames_num_cols),"cols"]
Stats_bad_cols <- Stats_common_colnames[Stats_common_colnames[,"Totals"]!=length(Stats_common_colnames_num_cols),"cols"]

#Stats_common_colnames <- rbind(Stats_common_colnames[Stats_common_colnames[,"Totals"]==length(Stats_common_colnames_num_cols),],
#                                 Stats_common_colnames[Stats_common_colnames[,"Totals"]!=length(Stats_common_colnames_num_cols),])
Stats_common_colnames <- rbind(Stats_common_colnames[Stats_common_colnames[,"cols"] %in% Stats_good_cols,],
                               Stats_common_colnames[Stats_common_colnames[,"cols"] %in% Stats_bad_cols,])
Stats_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Stats_common_colnames),1)
colnames(Stats_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Stats_common_colnames))
row.names(Stats_common_colnames) <- seq(nrow(Stats_common_colnames))

#rm2(Stats_pull,Stats_yr,Stats_files)
rm2(Stats_concatenate0)

for(i in which(sapply(Stats_concatenate,class)=="character"))
{
  Stats_concatenate[[i]] = trim(Stats_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Stats_concatenate))
{
  Stats_concatenate[,i] <- unknownToNA(Stats_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Stats_concatenate[,i] <- ifelse(is.na(Stats_concatenate[,i]),NA, Stats_concatenate[,i])
} 
rm2(i)

Stats_concatenate  <- Stats_concatenate[order(Stats_concatenate[,"pull"],Stats_concatenate[,"Fund_ID"],Stats_concatenate[,"Fund_Name"]),]
row.names(Stats_concatenate) <- seq(nrow(Stats_concatenate))

rm2(Stats_input)


#Reorder Columns

# XStats_concatenate_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
# XStats_concatenate_yearlyret_cols <- sort(colnames(Stats_concatenate)[grep("yearlyret", colnames(Stats_concatenate))])
# XStats_concatenate_monthlyret_cols <- sort(colnames(Stats_concatenate)[grep("monthlyret", colnames(Stats_concatenate))])
# XStats_concatenate_other_cols <- colnames(Stats_concatenate)[!(colnames(Stats_concatenate) %in% c(XStats_concatenate_id_cols,
#                                                                                                   XStats_concatenate_yearlyret_cols,
#                                                                                                   XStats_concatenate_monthlyret_cols))]
# 
# XStats_concatenate <- Stats_concatenate[,c(XStats_concatenate_id_cols,XStats_concatenate_yearlyret_cols,
#                                           XStats_concatenate_monthlyret_cols,XStats_concatenate_other_cols)]
# rm2(XStats_concatenate_id_cols,XStats_concatenate_yearlyret_cols,XStats_concatenate_monthlyret_cols,XStats_concatenate_other_cols)


Stats_concatenate_all_cols <- colnames(Stats_concatenate)

Stats_concatenate_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
Stats_concatenate_nonid_cols <- Stats_concatenate_all_cols[!(Stats_concatenate_all_cols %in% c(Stats_concatenate_id_cols))]

Stats_concatenate_yearlyret_cols <- Stats_concatenate_nonid_cols[grep("yearlyret", Stats_concatenate_nonid_cols)]
Stats_concatenate_nonyearlyret_cols <- Stats_concatenate_nonid_cols[!(Stats_concatenate_nonid_cols %in% c(Stats_concatenate_yearlyret_cols))]

Stats_concatenate_monthlyret_cols <- Stats_concatenate_nonyearlyret_cols[grep("monthlyret", Stats_concatenate_nonyearlyret_cols)]
Stats_concatenate_nonmonthlyret_cols <- Stats_concatenate_nonyearlyret_cols[!(Stats_concatenate_nonyearlyret_cols %in% c(Stats_concatenate_monthlyret_cols))]

Stats_concatenate <- Stats_concatenate[,c(Stats_concatenate_id_cols,
                                          sort(Stats_concatenate_yearlyret_cols),sort(Stats_concatenate_monthlyret_cols),
                                          Stats_concatenate_nonmonthlyret_cols)]

rm2(Stats_concatenate_all_cols)
rm2(Stats_concatenate_id_cols,Stats_concatenate_nonid_cols)
rm2(Stats_concatenate_yearlyret_cols,Stats_concatenate_nonyearlyret_cols)
rm2(Stats_concatenate_monthlyret_cols,Stats_concatenate_nonmonthlyret_cols)


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)

write.csv(Stats_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Stats_cols",".csv",sep=""),row.names=FALSE)

rm2(Stats_common_colnames,common_col_folder_path)
rm2(Stats_common_colnames_num_cols,Stats_good_cols,Stats_bad_cols)


#Check to see if nav_aum_ret folder exists.  If not, create it.
nav_aum_ret_folder_path <- paste(output_directory, "NAV_AUM_Ret", sep = "//", collapse = "//")  
create_directory(nav_aum_ret_folder_path,remove=1)

write.csv(Stats_concatenate, file=paste(nav_aum_ret_folder_path,"//","EurekahedgeHF_Stats_with_Rets",".csv",sep=""),row.names=FALSE)

rm2(Stats_concatenate,nav_aum_ret_folder_path)

