# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Service_Providers.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Combine Service Providers from Eurekahedge Data
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

Service_Providers_input0 <- data.frame(pull=NA,file_name=NA,read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                             stringsAsFactors=FALSE)
Service_Providers_input0[,"file_name"] <- Service_Providers_input0[,"file_clean"]
Service_Providers_input0[,"file_name"] <- gsub("\\\\","/",Service_Providers_input0[,"file_name"])
Service_Providers_input0[,"file_name"] <- gsub("//","/",Service_Providers_input0[,"file_name"])
Service_Providers_input0[,"file_name"] <- gsub("//","/",Service_Providers_input0[,"file_name"])

Service_Providers_input0[,"file_name"] <- encodeString(Service_Providers_input0[,"file_name"])

#Service_Providers_input0[,"pull"] <- regexpr("/[^/]*$", Service_Providers_input0[,"file_name"])
Service_Providers_input0[,"pull"] <- sapply(gregexpr("\\/", Service_Providers_input0[,"file_name"]), tail, 1)
Service_Providers_input0[,"file_name"] <- substr(Service_Providers_input0[,"file_name"],Service_Providers_input0[,"pull"]+1,nchar(Service_Providers_input0[,"file_name"]))

Service_Providers_input0[,"pull"] <- Service_Providers_input0[,"file_name"]
Service_Providers_input0[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=Service_Providers_input0[,"pull"])

Service_Providers_input <- Service_Providers_input0[grep("Service_Provider",Service_Providers_input0[,"file_name"]),]
row.names(Service_Providers_input) <- seq(nrow(Service_Providers_input))
rm2(Service_Providers_input0)


# Service_Providers_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
# Service_Providers_yr <- c(2014,2013,2012,2011,2010,2009,2008,2007)
# 
# Service_Providers_files <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Service_Providers","EurekahedgeHF_EXCEL_2013Feb_Data_Service_Providers",
#                              "EurekahedgeHF_EXCEL_2012Feb_Data_Service_Providers","EurekahedgeHF_EXCEL_2011Feb_Data_Service_Providers",
#                              "EurekahedgeHF_EXCEL_2010Feb_Data_Service_Providers","EurekahedgeHF_EXCEL_2009Feb_Data_Service_Providers",
#                              "EurekahedgeHF_EXCEL_2008Feb_Data_Service_Providers","EurekahedgeHF_EXCEL_2007Feb_Data_Service_Providers")
# 
# Service_Providers_input <- data.frame(matrix(NA, ncol=3, nrow=length(Service_Providers_files), dimnames=list(c(), c("pull","yr","file"))), 
#                                       stringsAsFactors=FALSE)
# 
# Service_Providers_input[,"pull"] <- Service_Providers_pull
# Service_Providers_input[,"yr"] <- Service_Providers_yr
# Service_Providers_input[,"file"] <- Service_Providers_files


Service_Providers_concatenate0 <- alply(.data=Service_Providers_input, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Service_Providers_input[1,]
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
  
  input_cols <- gsub(pattern="_Broker_Broker", replacement="_Broker", x=input_cols)
  
  input_cols <- gsub(pattern="Secondary_Prime_Broker", replacement="temp_col_secondary", x=input_cols)
  input_cols <- gsub(pattern="Synthetic_Prime_Broker", replacement="temp_col_synthetic", x=input_cols)
  
  input_cols <- gsub(pattern="Executive_Prime_Broker", replacement="temp_col_main", x=input_cols)
  input_cols <- gsub(pattern="Principal_Prime_Broker", replacement="temp_col_main", x=input_cols)
  input_cols <- gsub(pattern="Prime_Broker", replacement="temp_col_main", x=input_cols)
  
  input_cols <- gsub(pattern="temp_col_main", replacement="Principal_Prime_Broker_combcol", x=input_cols)
  input_cols <- gsub(pattern="temp_col_secondary", replacement="Secondary_Prime_Broker", x=input_cols)
  input_cols <- gsub(pattern="temp_col_synthetic", replacement="Synthetic_Prime_Broker", x=input_cols)
  
  input_cols <- gsub(pattern="offshore", replacement="Offshore", x=input_cols)
  input_cols <- gsub(pattern="onshore", replacement="Onshore", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
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

Service_Providers_concatenate <- rbind.fill(Service_Providers_concatenate0)

#Get colnames for all input files

Service_Providers_common_colnames <- ldply(.data=Service_Providers_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Service_Providers_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                      reshape(Service_Providers_common_colnames[,!(colnames(Service_Providers_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                      Totals=NA,stringsAsFactors=FALSE)
Service_Providers_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Service_Providers_common_colnames),1) 
Service_Providers_common_colnames_num_cols <- colnames(Service_Providers_common_colnames)[!(colnames(Service_Providers_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Service_Providers_common_colnames[,"Totals"] <- rowMeans(Service_Providers_common_colnames[,Service_Providers_common_colnames_num_cols],na.rm=TRUE)
Service_Providers_common_colnames <- Service_Providers_common_colnames[order(Service_Providers_common_colnames[,"Totals"],Service_Providers_common_colnames[,"Order_All_Org"]),]
Service_Providers_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Service_Providers_common_colnames),1)
Service_Providers_common_colnames[,"Totals"] <- rowSums(!is.na(Service_Providers_common_colnames[,Service_Providers_common_colnames_num_cols]))

Service_Providers_good_cols <- Service_Providers_common_colnames[Service_Providers_common_colnames[,"Totals"]==length(Service_Providers_common_colnames_num_cols),"cols"]
Service_Providers_bad_cols <- Service_Providers_common_colnames[Service_Providers_common_colnames[,"Totals"]!=length(Service_Providers_common_colnames_num_cols),"cols"]

#Service_Providers_common_colnames <- rbind(Service_Providers_common_colnames[Service_Providers_common_colnames[,"Totals"]==length(Service_Providers_common_colnames_num_cols),],
#                                 Service_Providers_common_colnames[Service_Providers_common_colnames[,"Totals"]!=length(Service_Providers_common_colnames_num_cols),])
Service_Providers_common_colnames <- rbind(Service_Providers_common_colnames[Service_Providers_common_colnames[,"cols"] %in% Service_Providers_good_cols,],
                                           Service_Providers_common_colnames[Service_Providers_common_colnames[,"cols"] %in% Service_Providers_bad_cols,])
Service_Providers_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Service_Providers_common_colnames),1)
colnames(Service_Providers_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Service_Providers_common_colnames))
row.names(Service_Providers_common_colnames) <- seq(nrow(Service_Providers_common_colnames))

#rm2(Service_Providers_pull,Service_Providers_yr,Service_Providers_files)
rm2(Service_Providers_concatenate0)

for(i in which(sapply(Service_Providers_concatenate,class)=="character"))
{
  Service_Providers_concatenate[[i]] = trim(Service_Providers_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Service_Providers_concatenate))
{
  Service_Providers_concatenate[,i] <- unknownToNA(Service_Providers_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Service_Providers_concatenate[,i] <- ifelse(is.na(Service_Providers_concatenate[,i]),NA, Service_Providers_concatenate[,i])
} 
rm2(i)

Service_Providers_concatenate  <- Service_Providers_concatenate[order(Service_Providers_concatenate[,"pull"],Service_Providers_concatenate[,"Fund_ID"],Service_Providers_concatenate[,"Fund_Name"]),]
row.names(Service_Providers_concatenate) <- seq(nrow(Service_Providers_concatenate))

rm2(Service_Providers_input)


#Reorder Columns

Service_Providers_concatenate_all_cols <- colnames(Service_Providers_concatenate)

Service_Providers_concatenate_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
Service_Providers_concatenate_nonid_cols <- Service_Providers_concatenate_all_cols[!(Service_Providers_concatenate_all_cols %in% c(Service_Providers_concatenate_id_cols))]

Service_Providers_concatenate_primebroker_col <- Service_Providers_concatenate_nonid_cols[grep("Prime_Broker", Service_Providers_concatenate_nonid_cols)]
Service_Providers_concatenate_nonprimebroker_col <- Service_Providers_concatenate_nonid_cols[!(Service_Providers_concatenate_nonid_cols %in% c(Service_Providers_concatenate_primebroker_col))]

Service_Providers_concatenate_legaladvisor_col <- Service_Providers_concatenate_nonprimebroker_col[grep("Legal_Advisor", Service_Providers_concatenate_nonprimebroker_col)]
Service_Providers_concatenate_nonlegaladvisor_col <- Service_Providers_concatenate_nonprimebroker_col[!(Service_Providers_concatenate_nonprimebroker_col %in% c(Service_Providers_concatenate_legaladvisor_col))]


Service_Providers_concatenate <- Service_Providers_concatenate[,c(Service_Providers_concatenate_id_cols,Service_Providers_concatenate_nonlegaladvisor_col,
                                                                  Service_Providers_concatenate_primebroker_col,Service_Providers_concatenate_legaladvisor_col)]

rm2(Service_Providers_concatenate_all_cols)
rm2(Service_Providers_concatenate_id_cols,Service_Providers_concatenate_nonid_cols)
rm2(Service_Providers_concatenate_primebroker_col,Service_Providers_concatenate_nonprimebroker_col)
rm2(Service_Providers_concatenate_legaladvisor_col,Service_Providers_concatenate_nonlegaladvisor_col)


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)

write.csv(Service_Providers_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Service_Providers_cols",".csv",sep=""),row.names=FALSE)

rm2(Service_Providers_common_colnames,common_col_folder_path)
rm2(Service_Providers_common_colnames_num_cols,Service_Providers_good_cols,Service_Providers_bad_cols)


#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

write.csv(Service_Providers_concatenate, file=paste(final_folder_path,"//","EurekahedgeHF_Service_Providers",".csv",sep=""),row.names=FALSE)

rm2(Service_Providers_concatenate,final_folder_path)
