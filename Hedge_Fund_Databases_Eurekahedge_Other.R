# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Other.R
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
                       "Hmisc","koRpus","lubridate","mitools","pbapply","plyr","R.oo","reshape","reshape2","rJava","RWeka","RWekajars",
                       "sqldf","stringr","tcltk","tm","zoo")
invisible(unlist(sapply(external_packages,load_external_packages, repo_str=repo, simplify=FALSE, USE.NAMES=FALSE)))
installed_packages <- list_installed_packages(external_packages)

rm(external_packages,installed_packages,repo)


###############################################################################
cat("SECTION: GET FILE TYPES", "\n")
###############################################################################

Other_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Fund_Details_Other_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)

Other_manager_details <- Other_input[grep("Manager_Details",Other_input[,"pull"],invert=FALSE),]
Other_nonmanager_details <- Other_input[grep("Manager_Details",Other_input[,"pull"],invert=TRUE),]

Other_service_provider <- Other_nonmanager_details[grep("Service_Provider",Other_nonmanager_details[,"pull"],invert=FALSE),]
Other_nonservice_provider <- Other_nonmanager_details[grep("Service_Provider",Other_nonmanager_details[,"pull"],invert=TRUE),]

Other_industry_focus <- Other_nonservice_provider[grep("Industry_Focus",Other_nonservice_provider[,"pull"],invert=FALSE),]
Other_nonindustry_focus <- Other_nonservice_provider[grep("Industry_Focus",Other_nonservice_provider[,"pull"],invert=TRUE),]

Other_country_focus <- Other_nonindustry_focus[grep("Country_Focus",Other_nonindustry_focus[,"pull"],invert=FALSE),]
Other_noncountry_focus <- Other_nonindustry_focus[grep("Country_Focus",Other_nonindustry_focus[,"pull"],invert=TRUE),]

#rm2(Other_input)
rm2(Other_nonmanager_details,Other_nonservice_provider,Other_nonindustry_focus)


###############################################################################
cat("SECTION: IMPORT MANAGER DETAILS", "\n")
###############################################################################

Other_manager_details_concatenate0 <- alply(.data=Other_manager_details, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Other_manager_details[1,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input_cols <- gsub(pattern="Countries_Focus", replacement="temp_col_focus_location", x=input_cols)
  input_cols <- gsub(pattern="Country_Focus", replacement="temp_col_focus_location", x=input_cols)
  
  input_cols <- gsub(pattern="Countries", replacement="temp_col_location", x=input_cols)
  input_cols <- gsub(pattern="Country", replacement="temp_col_location", x=input_cols)
  
  input_cols <- gsub(pattern="temp_col_focus_location", replacement="Country_Focus_combcol", x=input_cols)
  input_cols <- gsub(pattern="temp_col_location", replacement="Country_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
  input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  
  input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
  input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Other_manager_details_concatenate <- data.frame(pull_trim=NA,pull_trim2=NA,
                                                rbind.fill(Other_manager_details_concatenate0),
                                                stringsAsFactors=FALSE)

Other_manager_details_concatenate[,"pull_trim"] <- Other_manager_details_concatenate[,"pull"]
Other_manager_details_concatenate[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", Other_manager_details_concatenate[,"pull_trim"])

Other_manager_details_concatenate[,"pull_trim2"] <- Other_manager_details_concatenate[,"pull"]
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Fund_Details", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Manager_Details", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Providers", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Provider", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Industry_Focus", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Country_Focus", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])
Other_manager_details_concatenate[,"pull_trim2"] <- gsub(pattern="_Countries", replacement="", x=Other_manager_details_concatenate[,"pull_trim2"])

Other_manager_details_concatenate <- Other_manager_details_concatenate[,c("pull_trim","pull_trim2","pull",
                                                                          colnames(Other_manager_details_concatenate)[!(colnames(Other_manager_details_concatenate) %in% c("pull_trim","pull_trim2","pull"))])]

#Get colnames for all input files

Other_manager_details_common_colnames <- ldply(.data=Other_manager_details_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Other_manager_details_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                                    reshape(Other_manager_details_common_colnames[,!(colnames(Other_manager_details_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                                    Totals=NA,stringsAsFactors=FALSE)
Other_manager_details_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Other_manager_details_common_colnames),1) 
Other_manager_details_common_colnames_num_cols <- colnames(Other_manager_details_common_colnames)[!(colnames(Other_manager_details_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Other_manager_details_common_colnames[,"Totals"] <- rowMeans(Other_manager_details_common_colnames[,Other_manager_details_common_colnames_num_cols],na.rm=TRUE)
Other_manager_details_common_colnames <- Other_manager_details_common_colnames[order(Other_manager_details_common_colnames[,"Totals"],Other_manager_details_common_colnames[,"Order_All_Org"]),]
Other_manager_details_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Other_manager_details_common_colnames),1)
Other_manager_details_common_colnames[,"Totals"] <- rowSums(!is.na(Other_manager_details_common_colnames[,Other_manager_details_common_colnames_num_cols]))

Other_manager_details_good_cols <- Other_manager_details_common_colnames[Other_manager_details_common_colnames[,"Totals"]==length(Other_manager_details_common_colnames_num_cols),"cols"]
Other_manager_details_bad_cols <- Other_manager_details_common_colnames[Other_manager_details_common_colnames[,"Totals"]!=length(Other_manager_details_common_colnames_num_cols),"cols"]

#Other_manager_details_common_colnames <- rbind(Other_manager_details_common_colnames[Other_manager_details_common_colnames[,"Totals"]==length(Other_manager_details_common_colnames_num_cols),],
#                                 Other_manager_details_common_colnames[Other_manager_details_common_colnames[,"Totals"]!=length(Other_manager_details_common_colnames_num_cols),])
Other_manager_details_common_colnames <- rbind(Other_manager_details_common_colnames[Other_manager_details_common_colnames[,"cols"] %in% Other_manager_details_good_cols,],
                                               Other_manager_details_common_colnames[Other_manager_details_common_colnames[,"cols"] %in% Other_manager_details_bad_cols,])
Other_manager_details_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Other_manager_details_common_colnames),1)
colnames(Other_manager_details_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Other_manager_details_common_colnames))
row.names(Other_manager_details_common_colnames) <- seq(nrow(Other_manager_details_common_colnames))

rm2(Other_manager_details_concatenate0)

for(i in which(sapply(Other_manager_details_concatenate,class)=="character"))
{
  Other_manager_details_concatenate[[i]] = trim(Other_manager_details_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Other_manager_details_concatenate))
{
  Other_manager_details_concatenate[,i] <- unknownToNA(Other_manager_details_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Other_manager_details_concatenate[,i] <- ifelse(is.na(Other_manager_details_concatenate[,i]),NA, Other_manager_details_concatenate[,i])
} 
rm2(i)

Other_manager_details_concatenate  <- Other_manager_details_concatenate[order(Other_manager_details_concatenate[,"pull"],Other_manager_details_concatenate[,"Fund_ID"],Other_manager_details_concatenate[,"Fund_Name"]),]
row.names(Other_manager_details_concatenate) <- seq(nrow(Other_manager_details_concatenate))

#rm2(Other_manager_details)
rm2(Other_manager_details_common_colnames_num_cols,Other_manager_details_good_cols,Other_manager_details_bad_cols)


###############################################################################
cat("SECTION: IMPORT SERVICE PROVIDER", "\n")
###############################################################################

Other_service_provider_concatenate0 <- alply(.data=Other_service_provider, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Other_service_provider[1,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input_cols <- gsub(pattern="Countries_Focus", replacement="temp_col_focus_location", x=input_cols)
  input_cols <- gsub(pattern="Country_Focus", replacement="temp_col_focus_location", x=input_cols)
  
  input_cols <- gsub(pattern="Countries", replacement="temp_col_location", x=input_cols)
  input_cols <- gsub(pattern="Country", replacement="temp_col_location", x=input_cols)
  
  input_cols <- gsub(pattern="temp_col_focus_location", replacement="Country_Focus_combcol", x=input_cols)
  input_cols <- gsub(pattern="temp_col_location", replacement="Country_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
  input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  
  input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
  input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Other_service_provider_concatenate <- data.frame(pull_trim=NA,pull_trim2=NA,
                                                 rbind.fill(Other_service_provider_concatenate0),
                                                 stringsAsFactors=FALSE)

Other_service_provider_concatenate[,"pull_trim"] <- Other_service_provider_concatenate[,"pull"]
Other_service_provider_concatenate[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", Other_service_provider_concatenate[,"pull_trim"])

Other_service_provider_concatenate[,"pull_trim2"] <- Other_service_provider_concatenate[,"pull"]
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Fund_Details", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Manager_Details", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Providers", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Provider", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Industry_Focus", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Country_Focus", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])
Other_service_provider_concatenate[,"pull_trim2"] <- gsub(pattern="_Countries", replacement="", x=Other_service_provider_concatenate[,"pull_trim2"])

Other_service_provider_concatenate <- Other_service_provider_concatenate[,c("pull_trim","pull_trim2","pull",
                                                                            colnames(Other_service_provider_concatenate)[!(colnames(Other_service_provider_concatenate) %in% c("pull_trim","pull_trim2","pull"))])]

#Get colnames for all input files

Other_service_provider_common_colnames <- ldply(.data=Other_service_provider_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Other_service_provider_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                                     reshape(Other_service_provider_common_colnames[,!(colnames(Other_service_provider_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                                     Totals=NA,stringsAsFactors=FALSE)
Other_service_provider_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Other_service_provider_common_colnames),1) 
Other_service_provider_common_colnames_num_cols <- colnames(Other_service_provider_common_colnames)[!(colnames(Other_service_provider_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Other_service_provider_common_colnames[,"Totals"] <- rowMeans(Other_service_provider_common_colnames[,Other_service_provider_common_colnames_num_cols],na.rm=TRUE)
Other_service_provider_common_colnames <- Other_service_provider_common_colnames[order(Other_service_provider_common_colnames[,"Totals"],Other_service_provider_common_colnames[,"Order_All_Org"]),]
Other_service_provider_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Other_service_provider_common_colnames),1)
Other_service_provider_common_colnames[,"Totals"] <- rowSums(!is.na(Other_service_provider_common_colnames[,Other_service_provider_common_colnames_num_cols]))

Other_service_provider_good_cols <- Other_service_provider_common_colnames[Other_service_provider_common_colnames[,"Totals"]==length(Other_service_provider_common_colnames_num_cols),"cols"]
Other_service_provider_bad_cols <- Other_service_provider_common_colnames[Other_service_provider_common_colnames[,"Totals"]!=length(Other_service_provider_common_colnames_num_cols),"cols"]

#Other_service_provider_common_colnames <- rbind(Other_service_provider_common_colnames[Other_service_provider_common_colnames[,"Totals"]==length(Other_service_provider_common_colnames_num_cols),],
#                                 Other_service_provider_common_colnames[Other_service_provider_common_colnames[,"Totals"]!=length(Other_service_provider_common_colnames_num_cols),])
Other_service_provider_common_colnames <- rbind(Other_service_provider_common_colnames[Other_service_provider_common_colnames[,"cols"] %in% Other_service_provider_good_cols,],
                                                Other_service_provider_common_colnames[Other_service_provider_common_colnames[,"cols"] %in% Other_service_provider_bad_cols,])
Other_service_provider_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Other_service_provider_common_colnames),1)
colnames(Other_service_provider_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Other_service_provider_common_colnames))
row.names(Other_service_provider_common_colnames) <- seq(nrow(Other_service_provider_common_colnames))

rm2(Other_service_provider_concatenate0)

for(i in which(sapply(Other_service_provider_concatenate,class)=="character"))
{
  Other_service_provider_concatenate[[i]] = trim(Other_service_provider_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Other_service_provider_concatenate))
{
  Other_service_provider_concatenate[,i] <- unknownToNA(Other_service_provider_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Other_service_provider_concatenate[,i] <- ifelse(is.na(Other_service_provider_concatenate[,i]),NA, Other_service_provider_concatenate[,i])
} 
rm2(i)

Other_service_provider_concatenate  <- Other_service_provider_concatenate[order(Other_service_provider_concatenate[,"pull"],Other_service_provider_concatenate[,"Fund_ID"],Other_service_provider_concatenate[,"Fund_Name"]),]
row.names(Other_service_provider_concatenate) <- seq(nrow(Other_service_provider_concatenate))

#rm2(Other_service_provider)
rm2(Other_service_provider_common_colnames_num_cols,Other_service_provider_good_cols,Other_service_provider_bad_cols)


###############################################################################
cat("SECTION: IMPORT INDUSTRY FOCUS", "\n")
###############################################################################

Other_industry_focus_concatenate0 <- alply(.data=Other_industry_focus, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Other_industry_focus[1,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input_cols <- gsub(pattern="Countries_Focus", replacement="temp_col_focus_location", x=input_cols)
  input_cols <- gsub(pattern="Country_Focus", replacement="temp_col_focus_location", x=input_cols)
  
  input_cols <- gsub(pattern="Countries", replacement="temp_col_location", x=input_cols)
  input_cols <- gsub(pattern="Country", replacement="temp_col_location", x=input_cols)
  
  input_cols <- gsub(pattern="temp_col_focus_location", replacement="Country_Focus_combcol", x=input_cols)
  input_cols <- gsub(pattern="temp_col_location", replacement="Country_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
  input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  
  input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
  input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Other_industry_focus_concatenate <- data.frame(pull_trim=NA,pull_trim2=NA,
                                               rbind.fill(Other_industry_focus_concatenate0),
                                               stringsAsFactors=FALSE)

Other_industry_focus_concatenate[,"pull_trim"] <- Other_industry_focus_concatenate[,"pull"]
Other_industry_focus_concatenate[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", Other_industry_focus_concatenate[,"pull_trim"])

Other_industry_focus_concatenate[,"pull_trim2"] <- Other_industry_focus_concatenate[,"pull"]
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Fund_Details", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Manager_Details", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Providers", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Provider", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Industry_Focus", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Country_Focus", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])
Other_industry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Countries", replacement="", x=Other_industry_focus_concatenate[,"pull_trim2"])

Other_industry_focus_concatenate <- Other_industry_focus_concatenate[,c("pull_trim","pull_trim2","pull",
                                                                        colnames(Other_industry_focus_concatenate)[!(colnames(Other_industry_focus_concatenate) %in% c("pull_trim","pull_trim2","pull"))])]

#Get colnames for all input files

Other_industry_focus_common_colnames <- ldply(.data=Other_industry_focus_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Other_industry_focus_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                                   reshape(Other_industry_focus_common_colnames[,!(colnames(Other_industry_focus_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                                   Totals=NA,stringsAsFactors=FALSE)
Other_industry_focus_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Other_industry_focus_common_colnames),1) 
Other_industry_focus_common_colnames_num_cols <- colnames(Other_industry_focus_common_colnames)[!(colnames(Other_industry_focus_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Other_industry_focus_common_colnames[,"Totals"] <- rowMeans(Other_industry_focus_common_colnames[,Other_industry_focus_common_colnames_num_cols],na.rm=TRUE)
Other_industry_focus_common_colnames <- Other_industry_focus_common_colnames[order(Other_industry_focus_common_colnames[,"Totals"],Other_industry_focus_common_colnames[,"Order_All_Org"]),]
Other_industry_focus_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Other_industry_focus_common_colnames),1)
Other_industry_focus_common_colnames[,"Totals"] <- rowSums(!is.na(Other_industry_focus_common_colnames[,Other_industry_focus_common_colnames_num_cols]))

Other_industry_focus_good_cols <- Other_industry_focus_common_colnames[Other_industry_focus_common_colnames[,"Totals"]==length(Other_industry_focus_common_colnames_num_cols),"cols"]
Other_industry_focus_bad_cols <- Other_industry_focus_common_colnames[Other_industry_focus_common_colnames[,"Totals"]!=length(Other_industry_focus_common_colnames_num_cols),"cols"]

#Other_industry_focus_common_colnames <- rbind(Other_industry_focus_common_colnames[Other_industry_focus_common_colnames[,"Totals"]==length(Other_industry_focus_common_colnames_num_cols),],
#                                 Other_industry_focus_common_colnames[Other_industry_focus_common_colnames[,"Totals"]!=length(Other_industry_focus_common_colnames_num_cols),])
Other_industry_focus_common_colnames <- rbind(Other_industry_focus_common_colnames[Other_industry_focus_common_colnames[,"cols"] %in% Other_industry_focus_good_cols,],
                                              Other_industry_focus_common_colnames[Other_industry_focus_common_colnames[,"cols"] %in% Other_industry_focus_bad_cols,])
Other_industry_focus_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Other_industry_focus_common_colnames),1)
colnames(Other_industry_focus_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Other_industry_focus_common_colnames))
row.names(Other_industry_focus_common_colnames) <- seq(nrow(Other_industry_focus_common_colnames))

rm2(Other_industry_focus_concatenate0)

for(i in which(sapply(Other_industry_focus_concatenate,class)=="character"))
{
  Other_industry_focus_concatenate[[i]] = trim(Other_industry_focus_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Other_industry_focus_concatenate))
{
  Other_industry_focus_concatenate[,i] <- unknownToNA(Other_industry_focus_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Other_industry_focus_concatenate[,i] <- ifelse(is.na(Other_industry_focus_concatenate[,i]),NA, Other_industry_focus_concatenate[,i])
} 
rm2(i)

Other_industry_focus_concatenate  <- Other_industry_focus_concatenate[order(Other_industry_focus_concatenate[,"pull"],Other_industry_focus_concatenate[,"Fund_ID"],Other_industry_focus_concatenate[,"Fund_Name"]),]
row.names(Other_industry_focus_concatenate) <- seq(nrow(Other_industry_focus_concatenate))

#rm2(Other_industry_focus)
rm2(Other_industry_focus_common_colnames_num_cols,Other_industry_focus_good_cols,Other_industry_focus_bad_cols)


###############################################################################
cat("SECTION: IMPORT COUNTRY FOCUS", "\n")
###############################################################################

Other_country_focus_concatenate0 <- alply(.data=Other_country_focus, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Other_country_focus[1,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input_cols <- gsub(pattern="Countries_Focus", replacement="temp_col_focus_location", x=input_cols)
  input_cols <- gsub(pattern="Country_Focus", replacement="temp_col_focus_location", x=input_cols)
  
  input_cols <- gsub(pattern="Countries", replacement="temp_col_location", x=input_cols)
  input_cols <- gsub(pattern="Country", replacement="temp_col_location", x=input_cols)
  
  input_cols <- gsub(pattern="temp_col_focus_location", replacement="Country_Focus_combcol", x=input_cols)
  input_cols <- gsub(pattern="temp_col_location", replacement="Country_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
  input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  
  input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
  input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Other_country_focus_concatenate <- data.frame(pull_trim=NA,pull_trim2=NA,
                                              rbind.fill(Other_country_focus_concatenate0),
                                              stringsAsFactors=FALSE)

Other_country_focus_concatenate[,"pull_trim"] <- Other_country_focus_concatenate[,"pull"]
Other_country_focus_concatenate[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", Other_country_focus_concatenate[,"pull_trim"])

Other_country_focus_concatenate[,"pull_trim2"] <- Other_country_focus_concatenate[,"pull"]
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Fund_Details", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Manager_Details", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Providers", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Provider", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Industry_Focus", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Country_Focus", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])
Other_country_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Countries", replacement="", x=Other_country_focus_concatenate[,"pull_trim2"])

Other_country_focus_concatenate <- Other_country_focus_concatenate[,c("pull_trim","pull_trim2","pull",
                                                                      colnames(Other_country_focus_concatenate)[!(colnames(Other_country_focus_concatenate) %in% c("pull_trim","pull_trim2","pull"))])]

#Get colnames for all input files

Other_country_focus_common_colnames <- ldply(.data=Other_country_focus_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Other_country_focus_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                                  reshape(Other_country_focus_common_colnames[,!(colnames(Other_country_focus_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                                  Totals=NA,stringsAsFactors=FALSE)
Other_country_focus_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Other_country_focus_common_colnames),1) 
Other_country_focus_common_colnames_num_cols <- colnames(Other_country_focus_common_colnames)[!(colnames(Other_country_focus_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Other_country_focus_common_colnames[,"Totals"] <- rowMeans(Other_country_focus_common_colnames[,Other_country_focus_common_colnames_num_cols],na.rm=TRUE)
Other_country_focus_common_colnames <- Other_country_focus_common_colnames[order(Other_country_focus_common_colnames[,"Totals"],Other_country_focus_common_colnames[,"Order_All_Org"]),]
Other_country_focus_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Other_country_focus_common_colnames),1)
Other_country_focus_common_colnames[,"Totals"] <- rowSums(!is.na(Other_country_focus_common_colnames[,Other_country_focus_common_colnames_num_cols]))

Other_country_focus_good_cols <- Other_country_focus_common_colnames[Other_country_focus_common_colnames[,"Totals"]==length(Other_country_focus_common_colnames_num_cols),"cols"]
Other_country_focus_bad_cols <- Other_country_focus_common_colnames[Other_country_focus_common_colnames[,"Totals"]!=length(Other_country_focus_common_colnames_num_cols),"cols"]

#Other_country_focus_common_colnames <- rbind(Other_country_focus_common_colnames[Other_country_focus_common_colnames[,"Totals"]==length(Other_country_focus_common_colnames_num_cols),],
#                                 Other_country_focus_common_colnames[Other_country_focus_common_colnames[,"Totals"]!=length(Other_country_focus_common_colnames_num_cols),])
Other_country_focus_common_colnames <- rbind(Other_country_focus_common_colnames[Other_country_focus_common_colnames[,"cols"] %in% Other_country_focus_good_cols,],
                                             Other_country_focus_common_colnames[Other_country_focus_common_colnames[,"cols"] %in% Other_country_focus_bad_cols,])
Other_country_focus_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Other_country_focus_common_colnames),1)
colnames(Other_country_focus_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Other_country_focus_common_colnames))
row.names(Other_country_focus_common_colnames) <- seq(nrow(Other_country_focus_common_colnames))

rm2(Other_country_focus_concatenate0)

for(i in which(sapply(Other_country_focus_concatenate,class)=="character"))
{
  Other_country_focus_concatenate[[i]] = trim(Other_country_focus_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Other_country_focus_concatenate))
{
  Other_country_focus_concatenate[,i] <- unknownToNA(Other_country_focus_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Other_country_focus_concatenate[,i] <- ifelse(is.na(Other_country_focus_concatenate[,i]),NA, Other_country_focus_concatenate[,i])
} 
rm2(i)

Other_country_focus_concatenate  <- Other_country_focus_concatenate[order(Other_country_focus_concatenate[,"pull"],Other_country_focus_concatenate[,"Fund_ID"],Other_country_focus_concatenate[,"Fund_Name"]),]
row.names(Other_country_focus_concatenate) <- seq(nrow(Other_country_focus_concatenate))

#rm2(Other_country_focus)
rm2(Other_country_focus_common_colnames_num_cols,Other_country_focus_good_cols,Other_country_focus_bad_cols)


###############################################################################
cat("SECTION: IMPORT OTHER", "\n")
###############################################################################

Other_noncountry_focus_concatenate0 <- alply(.data=Other_noncountry_focus, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Other_noncountry_focus[1,]
  # directory_in <- output_directory
  # unknowns <- unknowns_strings
  
  input <- data.frame(pull=NA,
                      read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input_cols <- gsub(pattern="Countries_Focus", replacement="temp_col_focus_location", x=input_cols)
  input_cols <- gsub(pattern="Country_Focus", replacement="temp_col_focus_location", x=input_cols)
  
  input_cols <- gsub(pattern="Countries", replacement="temp_col_location", x=input_cols)
  input_cols <- gsub(pattern="Country", replacement="temp_col_location", x=input_cols)
  
  input_cols <- gsub(pattern="temp_col_focus_location", replacement="Country_Focus_combcol", x=input_cols)
  input_cols <- gsub(pattern="temp_col_location", replacement="Country_combcol", x=input_cols)
  
  input_cols <- gsub(pattern="_combcol_combcol", replacement="_combcol", x=input_cols)
  
  colnames(input) <- input_cols
  
  rm(input_cols)
  
  input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b-%y")
  input[,"Date_Added"] <- as.character(input[,"Date_Added"])
  
  input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b-%y")
  input[,"Dead_Date"] <- as.character(input[,"Dead_Date"])
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Other_noncountry_focus_concatenate <- data.frame(pull_trim=NA,pull_trim2=NA,
                                                 rbind.fill(Other_noncountry_focus_concatenate0),
                                                 stringsAsFactors=FALSE)

Other_noncountry_focus_concatenate[,"pull_trim"] <- Other_noncountry_focus_concatenate[,"pull"]
Other_noncountry_focus_concatenate[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", Other_noncountry_focus_concatenate[,"pull_trim"])

Other_noncountry_focus_concatenate[,"pull_trim2"] <- Other_noncountry_focus_concatenate[,"pull"]
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Fund_Details", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Manager_Details", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Providers", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Service_Provider", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Industry_Focus", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Country_Focus", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])
Other_noncountry_focus_concatenate[,"pull_trim2"] <- gsub(pattern="_Countries", replacement="", x=Other_noncountry_focus_concatenate[,"pull_trim2"])

Other_noncountry_focus_concatenate <- Other_noncountry_focus_concatenate[,c("pull_trim","pull_trim2","pull",
                                                                            colnames(Other_noncountry_focus_concatenate)[!(colnames(Other_noncountry_focus_concatenate) %in% c("pull_trim","pull_trim2","pull"))])]

#Get colnames for all input files

Other_noncountry_focus_common_colnames <- ldply(.data=Other_noncountry_focus_concatenate0, .fun = function(x){
  return(data.frame(cols=colnames(x),order=seq(1,length(colnames(x)),1), stringsAsFactors=FALSE))})
Other_noncountry_focus_common_colnames <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                                     reshape(Other_noncountry_focus_common_colnames[,!(colnames(Other_noncountry_focus_common_colnames) %in% c("pull","yr","month","file_org","file_clean"))], direction="wide",idvar=c("cols"),timevar="file_name"),
                                                     Totals=NA,stringsAsFactors=FALSE)
Other_noncountry_focus_common_colnames[,"Order_All_Org"] <- seq(1,nrow(Other_noncountry_focus_common_colnames),1) 
Other_noncountry_focus_common_colnames_num_cols <- colnames(Other_noncountry_focus_common_colnames)[!(colnames(Other_noncountry_focus_common_colnames) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
Other_noncountry_focus_common_colnames[,"Totals"] <- rowMeans(Other_noncountry_focus_common_colnames[,Other_noncountry_focus_common_colnames_num_cols],na.rm=TRUE)
Other_noncountry_focus_common_colnames <- Other_noncountry_focus_common_colnames[order(Other_noncountry_focus_common_colnames[,"Totals"],Other_noncountry_focus_common_colnames[,"Order_All_Org"]),]
Other_noncountry_focus_common_colnames[,"Order_All_Pos"] <- seq(1,nrow(Other_noncountry_focus_common_colnames),1)
Other_noncountry_focus_common_colnames[,"Totals"] <- rowSums(!is.na(Other_noncountry_focus_common_colnames[,Other_noncountry_focus_common_colnames_num_cols]))

Other_noncountry_focus_good_cols <- Other_noncountry_focus_common_colnames[Other_noncountry_focus_common_colnames[,"Totals"]==length(Other_noncountry_focus_common_colnames_num_cols),"cols"]
Other_noncountry_focus_bad_cols <- Other_noncountry_focus_common_colnames[Other_noncountry_focus_common_colnames[,"Totals"]!=length(Other_noncountry_focus_common_colnames_num_cols),"cols"]

#Other_noncountry_focus_common_colnames <- rbind(Other_noncountry_focus_common_colnames[Other_noncountry_focus_common_colnames[,"Totals"]==length(Other_noncountry_focus_common_colnames_num_cols),],
#                                 Other_noncountry_focus_common_colnames[Other_noncountry_focus_common_colnames[,"Totals"]!=length(Other_noncountry_focus_common_colnames_num_cols),])
Other_noncountry_focus_common_colnames <- rbind(Other_noncountry_focus_common_colnames[Other_noncountry_focus_common_colnames[,"cols"] %in% Other_noncountry_focus_good_cols,],
                                                Other_noncountry_focus_common_colnames[Other_noncountry_focus_common_colnames[,"cols"] %in% Other_noncountry_focus_bad_cols,])
Other_noncountry_focus_common_colnames[,"Order_All_Tot"] <- seq(1,nrow(Other_noncountry_focus_common_colnames),1)
colnames(Other_noncountry_focus_common_colnames) <- gsub(pattern="order.", replacement="", x=colnames(Other_noncountry_focus_common_colnames))
row.names(Other_noncountry_focus_common_colnames) <- seq(nrow(Other_noncountry_focus_common_colnames))

rm2(Other_noncountry_focus_concatenate0)

for(i in which(sapply(Other_noncountry_focus_concatenate,class)=="character"))
{
  Other_noncountry_focus_concatenate[[i]] = trim(Other_noncountry_focus_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Other_noncountry_focus_concatenate))
{
  Other_noncountry_focus_concatenate[,i] <- unknownToNA(Other_noncountry_focus_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Other_noncountry_focus_concatenate[,i] <- ifelse(is.na(Other_noncountry_focus_concatenate[,i]),NA, Other_noncountry_focus_concatenate[,i])
} 
rm2(i)

Other_noncountry_focus_concatenate  <- Other_noncountry_focus_concatenate[order(Other_noncountry_focus_concatenate[,"pull"],Other_noncountry_focus_concatenate[,"Fund_ID"],Other_noncountry_focus_concatenate[,"Fund_Name"]),]
row.names(Other_noncountry_focus_concatenate) <- seq(nrow(Other_noncountry_focus_concatenate))

#rm2(Other_noncountry_focus)
rm2(Other_noncountry_focus_common_colnames_num_cols,Other_noncountry_focus_good_cols,Other_noncountry_focus_bad_cols)


###############################################################################
cat("SECTION: MERGE DATA", "\n")
###############################################################################

#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

Merge_IDs_cols_keep <- c("pull_trim","pull","Fund_ID","Dead_Date","yr","month","date","bad_min","bad_max","AUM")                        
Merge_IDs <- data.frame(read.csv(file=paste(final_folder_path,"\\","EurekahedgeHF_NAV_AUM_Ret",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)[,Merge_IDs_cols_keep],
                        stringsAsFactors=FALSE)

rm2(Merge_IDs_cols_keep)

colnames(Merge_IDs)[match("AUM",names(Merge_IDs))] <- "pull_trim2"
                         
#Merge_IDs <- data.frame(read.csv(file=paste(final_folder_path,"\\","EurekahedgeHF_NAV_AUM_Ret",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)
#Merge_IDs_cols_drop <- c("min_date","max_date","Monthly_Ret","Monthly_Ret2","Yearly_Ret2","AUM")
#Merge_IDs <- data.frame(pull_trim2=NA,Merge_IDs[,!(colnames(Merge_IDs) %in% Merge_IDs_cols_drop)],stringsAsFactors=FALSE)
#rm2(Merge_IDs_cols_drop)

Merge_IDs[,"pull_trim2"] <- Merge_IDs[,"pull"]
Merge_IDs[,"pull_trim2"] <- gsub(pattern="_NAV_AUM", replacement="", x=Merge_IDs[,"pull_trim2"])

Merge_IDs <- Merge_IDs[,c("pull_trim","pull_trim2","pull",
                          colnames(Merge_IDs)[!(colnames(Merge_IDs) %in% c("pull_trim","pull_trim2","pull"))])]

Merge_IDs[,"pull_trim"] <- as.character(Merge_IDs[,"pull_trim"])

#View(Merge_IDs[1:1000,])

column_files <- c("Other_manager_details_common_colnames","Other_service_provider_common_colnames",
                  "Other_industry_focus_common_colnames","Other_country_focus_common_colnames",
                  "Other_noncountry_focus_common_colnames")

#Get common Tables

common_files0 <- ldply(.data=column_files, .fun = function(x){
  
  # x <- column_files[[1]]
  
  temp <- colnames(get(x[[1]]))
  temp_trim <- temp[!(temp %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","cols","Totals"))]
  temp_trim2 <- data.frame(pull_trim=NA,pull_trim2=NA,pull=NA,file=temp_trim,type_overall=NA,type_colnames=x[[1]],type_concatenate=NA,flag=1,stringsAsFactors=FALSE)
  
  temp_trim2[,"pull"] <- temp_trim2[,"file"]
  temp_trim2[,"pull"] <- gsub(pattern="(.csv|.CSV)", replacement="", temp_trim2[,"pull"])
  
  temp_trim2[,"pull_trim"] <- temp_trim2[,"pull"]
  temp_trim2[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", temp_trim2[,"pull_trim"])
  
  temp_trim2[,"pull_trim2"] <- temp_trim2[,"pull"]
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Fund_Details", replacement="", x=temp_trim2[,"pull_trim2"])
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Manager_Details", replacement="", x=temp_trim2[,"pull_trim2"])
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Service_Providers", replacement="", x=temp_trim2[,"pull_trim2"])
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Service_Provider", replacement="", x=temp_trim2[,"pull_trim2"])
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Industry_Focus", replacement="", x=temp_trim2[,"pull_trim2"])
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Country_Focus", replacement="", x=temp_trim2[,"pull_trim2"])
  temp_trim2[,"pull_trim2"] <- gsub(pattern="_Countries", replacement="", x=temp_trim2[,"pull_trim2"])
  
  temp_trim2[,"type_overall"] <- temp_trim2[,"type_colnames"]
  temp_trim2[,"type_overall"] <- gsub(pattern="_common_colnames", replacement="", x=temp_trim2[,"type_overall"])
  
  temp_trim2[,"type_concatenate"] <- temp_trim2[,"type_colnames"]
  temp_trim2[,"type_concatenate"] <- gsub(pattern="_common_colnames", replacement="_concatenate", x=temp_trim2[,"type_concatenate"])
  
  temp_trim2 <- temp_trim2[order(temp_trim2[,"pull_trim"],temp_trim2[,"pull_trim2"]),]
  row.names(temp_trim2) <- seq(nrow(temp_trim2))
  
  return(temp_trim2)
  
}, .progress = "text")

common_files <- data.frame(reshape(common_files0[,!(colnames(common_files0) %in% c("pull","type_colnames","file","flag"))],
                                   direction = "wide", idvar=c("pull_trim","pull_trim2"), timevar="type_overall"),stringsAsFactors=FALSE)
common_files <- common_files[order(common_files[,"pull_trim"],common_files[,"pull_trim2"]),]
row.names(common_files) <- seq(nrow(common_files))

rm2(common_files0,column_files)


### Merge each pull by pull_trim2

Other_concatenate0 <- adply(.data=common_files, .margins=1, .fun = function(x,ids,merge_data){
  
  # x <- common_files[1,]
  # x <- common_files[6,]
  # x <- common_files[10,]
  # ids <- c("pull_trim","pull_trim2")
  # merge_data <- "Merge_IDs"
  
  require(data.table)
  
  cat(x[,"pull_trim2"], "\n")
  
  files_temp1 <- as.vector(t(x[,!(colnames(x) %in% c(ids))]))
  files_temp1_trim <- data.frame(file=files_temp1[!is.na(files_temp1)],row_str=NA,col_str=NA,file_str=NA,stringsAsFactors=FALSE)
  files_temp1_trim[,"row_str"] <- paste(files_temp1_trim[,"file"],"[",",","'pull_trim2'","]","==","'",x[,"pull_trim2"],"'",sep="")
  files_temp1_trim[,"col_str"] <- paste("!(colnames(",files_temp1_trim[,"file"],") %in% c('pull','Dead_Date'))",sep="")
  files_temp1_trim[,"file_str"] <- paste(files_temp1_trim[,"file"],"[",files_temp1_trim[,"row_str"],",",files_temp1_trim[,"col_str"],"]",sep="")
  
  rm(files_temp1)
  
  files_temp2 <- merge_data
  files_temp2_trim <- data.frame(file=files_temp2[!is.na(files_temp2)],row_str=NA,col_str=NA,file_str=NA,stringsAsFactors=FALSE)
  files_temp2_trim[,"row_str"] <- paste(files_temp2_trim[,"file"],"[",",","'pull_trim2'","]","==","'",x[,"pull_trim2"],"'",sep="")
  files_temp2_trim[,"col_str"] <- paste("!(colnames(",files_temp2_trim[,"file"],") %in% c('pull'))",sep="")
  files_temp2_trim[,"file_str"] <- paste(files_temp2_trim[,"file"],"[",files_temp2_trim[,"row_str"],",",files_temp2_trim[,"col_str"],"]",sep="")
  
  rm(files_temp2)
  
  files_temp_all <- rbind(files_temp2_trim,files_temp1_trim)
  
  rm(files_temp2_trim,files_temp1_trim)
  
  merge_temp <- eval(parse(text=files_temp_all[1,"file_str"]))
  
  #for (i in 1:nrow(files_temp_all)) {print(paste(files_temp_all[i,"file"],": nrow = ",nrow(eval(parse(text=files_temp_all[i,"file_str"]))),sep=""))}
  
  for (i in 2:nrow(files_temp_all)) {
    
    # i <- 2
    # i <- 3
    # i <- 4
    
    data_temp <- data.frame(file_str=rbind("merge_temp",files_temp_all[i,"file_str"]),stringsAsFactors=FALSE)
    
    common_ids1a <- adply(.data=data_temp, .margins=1, .fun = function(x){
      
      temp_order <- data.frame(cols=colnames(eval(parse(text=x[,"file_str"]))),order=NA,stringsAsFactors=FALSE)
      temp_order[,"order"] <- seq(1,nrow(temp_order))
      return(temp_order)
      
    }, .expand = FALSE, .progress = "none")
    
    common_ids1b <- ddply(.data=common_ids1a[,!(colnames(common_ids1a) %in% c("X1"))], .variables=c("cols"),.fun = function(x){
      
      return(data.frame(freq=nrow(x),avg_order=mean(x[,"order"]),stringsAsFactors=FALSE))
      
    }, .progress = "none")
    
    rm(common_ids1a)
    
    common_ids1b <- common_ids1b[order(common_ids1b[,"avg_order"],common_ids1b[,"freq"],common_ids1b[,"cols"]),]
    row.names(common_ids1b) <- seq(nrow(common_ids1b))

    common_ids1 <- common_ids1b[common_ids1b[,"freq"]==nrow(data_temp),] 
    
    rm(common_ids1b)
    
    #data_str1 <- paste("list(",paste(data_temp[,"file_str"], sep = "", collapse = ","),")",sep="")
    #merge_temp <- join_all(eval(parse(text=data_str1)), by = common_ids1[,"cols"], type = "left", match = "all")
    
    merge_temp <- merge(data.table(merge_temp, key=common_ids1[,"cols"]),
                        eval(parse(text=paste("data.table(",data_temp[2,"file_str"],",key=common_ids1[,'cols'])",sep=""))),
                        by.x=common_ids1[,"cols"], by.y=common_ids1[,"cols"], 
                        all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))

    invisible(gc(verbose = FALSE, reset = TRUE))
    
    rm(data_temp,common_ids1)
  }
  rm(i,files_temp_all)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
  order_ids <- c("pull_trim","pull_trim2",
                 "Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason",
                 "date","yr","month","bad_min","bad_max")
  
  merge_temp <- as.data.frame(merge_temp,stringsAsFactors=FALSE)
  
  merge_temp <- merge_temp[,c(order_ids,colnames(merge_temp)[!(colnames(merge_temp) %in% c(order_ids))])]
  
  merge_temp <- merge_temp[order(merge_temp[,"Fund_ID"],merge_temp[,"date"],
                                   merge_temp[,"pull_trim"],merge_temp[,"pull_trim2"]),]
  row.names(merge_temp) <- seq(nrow(merge_temp))

  rm(order_ids)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
  return(merge_temp)
  
}, ids=c("pull_trim","pull_trim2"), merge_data="Merge_IDs", .expand = TRUE,.progress = "text")

rm2(Merge_IDs)
rm2("Other_manager_details_concatenate","Other_service_provider_concatenate")
rm2("Other_industry_focus_concatenate","Other_country_focus_concatenate")
rm2("Other_noncountry_focus_concatenate")

Other_concatenate_drop_cols <- colnames(common_files)[!(colnames(common_files) %in% c("pull_trim","pull_trim2"))]
Other_concatenate <- Other_concatenate0[,!(colnames(Other_concatenate0) %in% Other_concatenate_drop_cols)]

rm2(Other_concatenate0,Other_concatenate_drop_cols)

Other_concatenate <- Other_concatenate[order(Other_concatenate[,"Fund_ID"],Other_concatenate[,"date"],
                                 Other_concatenate[,"pull_trim"],Other_concatenate[,"pull_trim2"]),]
row.names(Other_concatenate) <- seq(nrow(Other_concatenate))

invisible(gc(verbose = FALSE, reset = TRUE))

#str(Other_concatenate)
#View(Other_concatenate[1:1000,])


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)

write.csv(Other_manager_details_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Manager_Details_cols",".csv",sep=""),row.names=FALSE)
write.csv(Other_service_provider_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Service_Provider_cols",".csv",sep=""),row.names=FALSE)
write.csv(Other_industry_focus_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Industry_Focus_cols",".csv",sep=""),row.names=FALSE)
write.csv(Other_country_focus_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Country_Focus_cols",".csv",sep=""),row.names=FALSE)
write.csv(Other_noncountry_focus_common_colnames, file=paste(common_col_folder_path,"//","EurekahedgeHF_Noncountry_cols",".csv",sep=""),row.names=FALSE)

rm2(common_col_folder_path)
rm2(Other_manager_details_common_colnames,Other_service_provider_common_colnames)
rm2(Other_industry_focus_common_colnames,Other_country_focus_common_colnames)
rm2(Other_noncountry_focus_common_colnames)

write.csv(Other_concatenate, file=paste(final_folder_path,"//","EurekahedgeHF_Other",".csv",sep=""),row.names=FALSE)

rm2(Other_concatenate,final_folder_path)
rm2(Other_manager_details,Other_service_provider)
rm2(Other_industry_focus,Other_country_focus)
rm2(Other_noncountry_focus)

rm2(Other_input,common_files)
