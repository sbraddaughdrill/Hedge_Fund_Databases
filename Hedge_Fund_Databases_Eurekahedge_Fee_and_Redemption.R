# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Fee_and_Redemption.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Combine Fee and Redemption from Eurekahedge Data
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

Fee_Redemption_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_Fee_Redemption_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)

# Fee_Redemption_input0 <- data.frame(pull=NA,file_name=NA,read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
#                                     stringsAsFactors=FALSE)
# Fee_Redemption_input0[,"file_name"] <- Fee_Redemption_input0[,"file_clean"]
# Fee_Redemption_input0[,"file_name"] <- gsub("\\\\","/",Fee_Redemption_input0[,"file_name"])
# Fee_Redemption_input0[,"file_name"] <- gsub("//","/",Fee_Redemption_input0[,"file_name"])
# Fee_Redemption_input0[,"file_name"] <- gsub("//","/",Fee_Redemption_input0[,"file_name"])
# 
# Fee_Redemption_input0[,"file_name"] <- encodeString(Fee_Redemption_input0[,"file_name"])
# 
# #Fee_Redemption_input0[,"pull"] <- regexpr("/[^/]*$", Fee_Redemption_input0[,"file_name"])
# Fee_Redemption_input0[,"pull"] <- sapply(gregexpr("\\/", Fee_Redemption_input0[,"file_name"]), tail, 1)
# Fee_Redemption_input0[,"file_name"] <- substr(Fee_Redemption_input0[,"file_name"],Fee_Redemption_input0[,"pull"]+1,nchar(Fee_Redemption_input0[,"file_name"]))
# 
# Fee_Redemption_input0[,"pull"] <- Fee_Redemption_input0[,"file_name"]
# Fee_Redemption_input0[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=Fee_Redemption_input0[,"pull"])
# 
# Fee_Redemption_input <- Fee_Redemption_input0[grep("Fee_and_Redemption_Structure",Fee_Redemption_input0[,"file_name"]),]
# row.names(Fee_Redemption_input) <- seq(nrow(Fee_Redemption_input))
# rm2(Fee_Redemption_input0)

# Fee_Redemption_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
# Fee_Redemption_yr <- c(2014,2013,2012,2011,2010,2009,2008,2007)
# 
# Fee_Redemption_files <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Fee_and_Redemption","EurekahedgeHF_EXCEL_2013Feb_Data_Fee_and_Redemption",
#                           "EurekahedgeHF_EXCEL_2012Feb_Data_Fee_and_Redemption","EurekahedgeHF_EXCEL_2011Feb_Data_Fee_and_Redemption",
#                           "EurekahedgeHF_EXCEL_2010Feb_Data_Fee_and_Redemption","EurekahedgeHF_EXCEL_2009Feb_Data_Fee_and_Redemption",
#                           "EurekahedgeHF_EXCEL_2008Feb_Data_Fee_and_Redemption","EurekahedgeHF_EXCEL_2007Feb_Data_Fee_and_Redemption")
# 
# Fee_Redemption_input <- data.frame(matrix(NA, ncol=3, nrow=length(Fee_Redemption_files), dimnames=list(c(), c("pull","yr","file"))), 
#                                    stringsAsFactors=FALSE)
# 
# Fee_Redemption_input[,"pull"] <- Fee_Redemption_pull
# Fee_Redemption_input[,"yr"] <- Fee_Redemption_yr
# Fee_Redemption_input[,"file"] <- Fee_Redemption_files

Fee_Redemption_concatenate0 <- alply(.data=Fee_Redemption_input, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Fee_Redemption_input[1,]
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

Fee_Redemption_concatenate <- rbind.fill(Fee_Redemption_concatenate0)

#Get colnames for all input files

Fee_Redemption_common_colnames <- ldply(.data=Fee_Redemption_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Fee_Redemption_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                      reshape(Fee_Redemption_common_colnames[,!(colnames(Fee_Redemption_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                      Totals=NA,stringsAsFactors=FALSE)
Fee_Redemption_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Fee_Redemption_common_colnames),1) 
Fee_Redemption_common_colnames_num_cols <- colnames(Fee_Redemption_common_colnames)[!(colnames(Fee_Redemption_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Fee_Redemption_common_colnames[,"Totals"] <- rowMeans(Fee_Redemption_common_colnames[,Fee_Redemption_common_colnames_num_cols],na.rm=TRUE)
Fee_Redemption_common_colnames <- Fee_Redemption_common_colnames[order(Fee_Redemption_common_colnames[,"Totals"],Fee_Redemption_common_colnames[,"Order_All_Org"]),]
Fee_Redemption_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Fee_Redemption_common_colnames),1)
Fee_Redemption_common_colnames[,"Totals"] <- rowSums(!is.na(Fee_Redemption_common_colnames[,Fee_Redemption_common_colnames_num_cols]))

Fee_Redemption_good_cols <- Fee_Redemption_common_colnames[Fee_Redemption_common_colnames[,"Totals"]==length(Fee_Redemption_common_colnames_num_cols),"cols"]
Fee_Redemption_bad_cols <- Fee_Redemption_common_colnames[Fee_Redemption_common_colnames[,"Totals"]!=length(Fee_Redemption_common_colnames_num_cols),"cols"]

#Fee_Redemption_common_colnames <- rbind(Fee_Redemption_common_colnames[Fee_Redemption_common_colnames[,"Totals"]==length(Fee_Redemption_common_colnames_num_cols),],
#                                 Fee_Redemption_common_colnames[Fee_Redemption_common_colnames[,"Totals"]!=length(Fee_Redemption_common_colnames_num_cols),])
Fee_Redemption_common_colnames <- rbind(Fee_Redemption_common_colnames[Fee_Redemption_common_colnames[,"cols"] %in% Fee_Redemption_good_cols,],
                                        Fee_Redemption_common_colnames[Fee_Redemption_common_colnames[,"cols"] %in% Fee_Redemption_bad_cols,])
Fee_Redemption_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Fee_Redemption_common_colnames),1)
colnames(Fee_Redemption_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Fee_Redemption_common_colnames))
row.names(Fee_Redemption_common_colnames) <- seq(nrow(Fee_Redemption_common_colnames))

#rm2(Fee_Redemption_pull,Fee_Redemption_yr,Fee_Redemption_files)
rm2(Fee_Redemption_concatenate0)

for(i in which(sapply(Fee_Redemption_concatenate,class)=="character"))
{
  Fee_Redemption_concatenate[[i]] = trim(Fee_Redemption_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Fee_Redemption_concatenate))
{
  Fee_Redemption_concatenate[,i] <- unknownToNA(Fee_Redemption_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Fee_Redemption_concatenate[,i] <- ifelse(is.na(Fee_Redemption_concatenate[,i]),NA, Fee_Redemption_concatenate[,i])
} 
rm2(i)

Fee_Redemption_concatenate  <- Fee_Redemption_concatenate[order(Fee_Redemption_concatenate[,"pull"],Fee_Redemption_concatenate[,"Fund_ID"],Fee_Redemption_concatenate[,"Fund_Name"]),]
row.names(Fee_Redemption_concatenate) <- seq(nrow(Fee_Redemption_concatenate))

rm2(Fee_Redemption_input)


#Reorder Columns

Fee_Redemption_concatenate_all_cols <- colnames(Fee_Redemption_concatenate)

Fee_Redemption_concatenate_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
Fee_Redemption_concatenate_nonid_cols <- Fee_Redemption_concatenate_all_cols[!(Fee_Redemption_concatenate_all_cols %in% c(Fee_Redemption_concatenate_id_cols))]

Fee_Redemption_concatenate_single_col <- Fee_Redemption_concatenate_nonid_cols[grep("Key_Man_Clause", Fee_Redemption_concatenate_nonid_cols)]
Fee_Redemption_concatenate_nonsingle_col <- Fee_Redemption_concatenate_nonid_cols[!(Fee_Redemption_concatenate_nonid_cols %in% c(Fee_Redemption_concatenate_single_col))]

Fee_Redemption_concatenate <- Fee_Redemption_concatenate[,c(Fee_Redemption_concatenate_id_cols,
                                                            Fee_Redemption_concatenate_nonsingle_col,Fee_Redemption_concatenate_single_col)]

rm2(Fee_Redemption_concatenate_all_cols)
rm2(Fee_Redemption_concatenate_id_cols,Fee_Redemption_concatenate_nonid_cols)
rm2(Fee_Redemption_concatenate_single_col,Fee_Redemption_concatenate_nonsingle_col)


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)

write.csv(Fee_Redemption_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Fee_and_Redemption_cols",".csv",sep=""),row.names=FALSE)

rm2(Fee_Redemption_common_colnames,common_col_folder_path)
rm2(Fee_Redemption_common_colnames_num_cols,Fee_Redemption_good_cols,Fee_Redemption_bad_cols)


#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

write.csv(Fee_Redemption_concatenate, file=paste(final_folder_path,"//","EurekahedgeHF_Fee_and_Redemption",".csv",sep=""),row.names=FALSE)

rm2(Fee_Redemption_concatenate,final_folder_path)
