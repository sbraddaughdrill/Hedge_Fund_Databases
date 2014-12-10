# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Clean.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Clean the seperated Eurekahedge Data
#
#
# Run this file after all the sheets are seperated.
# See Macro in 'Sep' folder.
# Place files in 'Sep' folder.
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
cat("SECTION: GET LIST OF EACH FILE TYPE", "\n")
###############################################################################

sep_folder_path <- paste(output_directory, "Sep", sep = "//", collapse = "//")  

sep_folder_files0 <- data.frame(files=list.files(path=sep_folder_path),file_name=NA,yr=NA,month=NA,output_dir=NA,stringsAsFactors=FALSE)

#sep_folder_files <- sep_folder_files0[!grepl(".TXT|.txt", sep_folder_files0[,"files"]),]
sep_folder_files <- sep_folder_files0[grepl(".CSV|.csv", sep_folder_files0[,"files"]),]

sep_folder_files[,"file_name"] <- sep_folder_files[,"files"]
sep_folder_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=sep_folder_files[,"file_name"])

sep_folder_files[,"yr"] <- sep_folder_files[,"files"]
sep_folder_files[,"yr"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", x=sep_folder_files[,"yr"])
sep_folder_files[,"yr"] <- substr(sep_folder_files[,"yr"], 1, 4)
#sep_folder_files[,"yr"] <- as.integer(sep_folder_files[,"yr"])

sep_folder_files[,"month"] <- sep_folder_files[,"files"]
sep_folder_files[,"month"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", x=sep_folder_files[,"month"])
sep_folder_files[,"month"] <- substr(sep_folder_files[,"month"], 5, 6)
#sep_folder_files[,"month"] <- as.integer(sep_folder_files[,"month"])

sep_folder_files[,"output_dir"] <- paste(output_directory,"\\",sep_folder_files[,"yr"],sep_folder_files[,"month"],sep="")

files_fund_details <- sep_folder_files[grepl("Fund_Details", sep_folder_files[,"files"]),]
files_fund_details <- files_fund_details[order(files_fund_details[,"output_dir"],decreasing=c(T)),]
row.names(files_fund_details) <- seq(nrow(files_fund_details))

files_instruments <- sep_folder_files[grepl("Instruments_Traded", sep_folder_files[,"files"]),]
files_instruments <- files_instruments[order(files_instruments[,"output_dir"],decreasing=c(T)),]
row.names(files_instruments) <- seq(nrow(files_instruments))

files_NAV_AUM <- sep_folder_files[grepl("NAV_AUM", sep_folder_files[,"files"]),]
files_NAV_AUM <- files_NAV_AUM[order(files_NAV_AUM[,"output_dir"],decreasing=c(T)),]
row.names(files_NAV_AUM) <- seq(nrow(files_NAV_AUM))

if(nrow(sep_folder_files)!=(nrow(files_fund_details)+nrow(files_instruments)+nrow(files_NAV_AUM))){ print("Not all filenames found") }

rm2(sep_folder_files0,sep_folder_files)

#Check to see if common_cols folder exists.  If not, create it.
common_col_folder_path <- paste(output_directory, "Common_Cols", sep = "//", collapse = "//")  
create_directory(common_col_folder_path,remove=1)


###############################################################################
cat("SECTION: CLEAN FUND DETAILS FILES", "\n")
###############################################################################

fund_details_cat_cols_concatenate0 <- adply(.data=files_fund_details, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- files_fund_details[1,]
  # x <- files_fund_details[3,]
  # x <- files_fund_details[6,]
  # x <- files_fund_details[9,]
  # x <- files_fund_details[10,]
  # x <- files_fund_details[12,]
  # x <- files_fund_details[13,]
  # x <- files_fund_details[14,]
  # x <- files_fund_details[15,]
  # directory_in <- sep_folder_path
  # unknowns <- unknowns_strings
  
  #Create output directory
  output_folder_path <- paste(x[,"output_dir"], "", sep = "//", collapse = "//")  
  create_directory(output_folder_path,remove=1)
  
  outfile_prefix <- x[,"file_name"]
  
  input <- data.frame(read.csv(file=paste(directory_in,"\\",x[,"files"],sep=""),header=FALSE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input <- input[,colSums(is.na(input))<nrow(input)]
  
  # Get Columns and Category Names
  cat_col <- data.frame(order=NA,col=t(input[2,]),bad_col=NA,category_org=t(input[1,]),
                        category_id=NA,stringsAsFactors=FALSE)
  colnames(cat_col) <- c("order","col","bad_col","category_org","category_id")
  
  cat_col[,"order"] <- seq(nrow(cat_col))
  
  #cat_col <- cat_col[order(cat_col[,"order"]),]
  #row.names(cat_col) <- seq(nrow(cat_col))
  
  cat_col[,"bad_col"] <- ifelse(is.na(cat_col[,c("col")]),1,0)
  
  cat_col[,"category_id"] <- ifelse(is.na(cat_col[,c("category_org")]),0,1)
  cat_col[,"category_id"] <- cumsum(cat_col[,"category_id"])
  
  cat_col[,"category_id"] <- ifelse(is.na(cat_col[,c("col")]),NA,cat_col[,"category_id"])
  
  # Clean Column Names
  cat_col[,"col"] <- gsub(pattern="[[:punct:]]", replacement="", x=cat_col[,"col"])
  
  cat_col[,"col"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"col"])
  cat_col[,"col"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"col"])
  cat_col[,"col"] <- gsub("^\\s+|\\s+$", "", cat_col[,"col"])
  
  cat_col[,"col"] <- gsub(pattern=" ", replacement="_", x=cat_col[,"col"])
  
  # Clean Category Names
  cat_col[,"category_org"] <- gsub(pattern="[[:punct:]]", replacement="", x=cat_col[,"category_org"])
  
  cat_col[,"category_org"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"category_org"])
  cat_col[,"category_org"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"category_org"])
  cat_col[,"category_org"] <- gsub("^\\s+|\\s+$", "", cat_col[,"category_org"])
  
  cat_col[,"category_org"] <- gsub(pattern=" ", replacement="_", x=cat_col[,"category_org"])
  
  # Expand Category Names
  
  cat_id <- data.frame(category=unique(cat_col[!is.na(cat_col[,c("category_org")]),c("category_org")]),
                       category_id=NA,stringsAsFactors=FALSE)
  cat_id[,"category_id"] <- seq(nrow(cat_id))
  
  cat_col_id <- merge(cat_col, cat_id, 
                      by.x=c("category_id"),by.y=c("category_id"), 
                      all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))
  
  cat_col_id <- cat_col_id[order(cat_col_id[,"order"]),]
  row.names(cat_col_id) <- seq(nrow(cat_col_id))
  
  cat_col_trim <- cat_col_id[,!(colnames(cat_col_id) %in% c("category_id","category_org"))]
  
  cat_col_trim[,"category"] <- ifelse((cat_col_trim[,c("bad_col")]==0 & is.na(cat_col_trim[,c("category")])),
                                      "Common_IDs",cat_col_trim[,"category"])
  
  rm(cat_col,cat_col_id)
  
  ## REPLACE COLUMM NAMES
  
  input_trim <- input[3:nrow(input),]
  
  rm(input)
  
  colnames(input_trim) <- paste(cat_col_trim[,"col"],"",sep="")
  
  input_trim[,"Fund_ID"] <- as.integer(input_trim[,"Fund_ID"])
  input_trim <- input_trim[order(input_trim[,"Fund_ID"]),]
  row.names(input_trim) <- seq(nrow(input_trim))
  
  ## CHECK FOR BAD ROWS
  #colnames(input_trim)
  #aa0 <- unique(input_trim[,c(which(cat_col_trim[,"bad_col"]==1))])
  #aaa1 <- input_trim[row.names(aa0),c(1,which(cat_col_trim[,"bad_col"]==1))]
  #which(cat_col_trim[,"bad_col"]==1)
  #aa1 <- input_trim[(!is.na(input_trim[,75])|!is.na(input_trim[,76])),c(1,which(cat_col_trim[,"bad_col"]==1))]
  
  
  ## LOOP OVER CATEGORIES AND OUTPUT COMMONID AND CATEGORY
  
  out_files <- adply(.data=cat_id, .margins=1, .fun = function(y,data,col_cat,id_cat,directory_out,file_out){
    
    # y <- cat_id[1,]
    # data <- input_trim
    # col_cat <- cat_col_trim
    # id_cat <- "Common_IDs"
    # directory_out <- output_folder_path
    # file_out <- outfile_prefix
    
    cols_keep <- col_cat[col_cat[,"category"] %in% c(id_cat,y[,"category"]),"col"]
    
    x_out <- data[,colnames(data) %in% cols_keep]
    
    outname <- paste(directory_out,file_out,"_",y[,"category"],".csv",sep="")
    
    write.csv(x_out, file=outname,row.names=FALSE)
    
    return(outname)
    
  },data=input_trim,col_cat=cat_col_trim,id_cat="Common_IDs",directory_out=output_folder_path,file_out=outfile_prefix)
  
  colnames(out_files) <- c(colnames(cat_id),"file_clean")
  
  cat_col_trim2 <- merge(cat_col_trim,out_files[,(colnames(out_files) %in% c("category","file_clean"))], 
                         by.x=c("category"),by.y=c("category"), 
                         all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))
  
  cat_col_trim2 <- cat_col_trim2[order(cat_col_trim2[,"order"]),]
  
  row.names(cat_col_trim2) <- seq(nrow(cat_col_trim2))
  
  rm(cat_id,input_trim)
  
  rm(output_folder_path,outfile_prefix)
  
  gc()
  
  ##RETURN THE COLNAMES
  return(cat_col_trim2)
  
}, directory_in=sep_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

fund_details_cat_cols_concatenate0 <- fund_details_cat_cols_concatenate0[,c("files","file_name","file_clean",
                                                                            colnames(fund_details_cat_cols_concatenate0)[!(colnames(fund_details_cat_cols_concatenate0) %in% c("files","file_name","file_clean"))])]

### Find empty cols
#fund_details_empty_cols <- fund_details_cat_cols_concatenate0[is.na(fund_details_cat_cols_concatenate[,"col"]),]
fund_details_empty_cols <- fund_details_cat_cols_concatenate0[fund_details_cat_cols_concatenate0[,"bad_col"]==1,]

fund_details_cat_cols_concatenate <- fund_details_cat_cols_concatenate0[,!(colnames(fund_details_cat_cols_concatenate0) %in% c("file_clean"))]

write.csv(fund_details_cat_cols_concatenate, file=paste(common_col_folder_path,"//","EurekahedgeHF_Fund_Details_common_cols_org",".csv",sep=""),row.names=FALSE)

files_fund_details_out0 <- unique(fund_details_cat_cols_concatenate0[,c(colnames(files_fund_details),"file_clean")])

colnames(files_fund_details_out0)[match("files",names(files_fund_details_out0))] <- "file_org"
files_fund_details_out0[,"file_org"] <- paste(sep_folder_path,"//",files_fund_details_out0[,"file_org"],sep="")
files_fund_details_out1 <- files_fund_details_out0[!is.na(files_fund_details_out0[,"file_clean"]),]

files_fund_details_out2 <- unique(files_fund_details_out1[,c("yr","month","file_org","file_clean")])
row.names(files_fund_details_out2) <- seq(nrow(files_fund_details_out2))

#write.csv(files_fund_details_out, file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_files",".csv",sep=""),row.names=FALSE)

files_fund_details_out3 <- data.frame(pull=NA,file_name=NA,files_fund_details_out2,stringsAsFactors=FALSE)
files_fund_details_out3[,"file_name"] <- files_fund_details_out3[,"file_clean"]
files_fund_details_out3[,"file_name"] <- gsub("\\\\","/",files_fund_details_out3[,"file_name"])
files_fund_details_out3[,"file_name"] <- gsub("//","/",files_fund_details_out3[,"file_name"])
files_fund_details_out3[,"file_name"] <- gsub("//","/",files_fund_details_out3[,"file_name"])

files_fund_details_out3[,"file_name"] <- encodeString(files_fund_details_out3[,"file_name"])

#files_fund_details_out3[,"pull"] <- regexpr("/[^/]*$", files_fund_details_out3[,"file_name"])
files_fund_details_out3[,"pull"] <- sapply(gregexpr("\\/", files_fund_details_out3[,"file_name"]), tail, 1)
files_fund_details_out3[,"file_name"] <- substr(files_fund_details_out3[,"file_name"],files_fund_details_out3[,"pull"]+1,nchar(files_fund_details_out3[,"file_name"]))

files_fund_details_out3[,"pull"] <- files_fund_details_out3[,"file_name"]
files_fund_details_out3[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=files_fund_details_out3[,"pull"])

rm2(files_fund_details,files_fund_details_out0,files_fund_details_out1,files_fund_details_out2)

files_fund_details_all <- files_fund_details_out3

files_fund_details_stats <- files_fund_details_all[grep("Statistics",files_fund_details_all[,"pull"],invert=FALSE),]
files_fund_details_nonstats <- files_fund_details_all[grep("Statistics",files_fund_details_all[,"pull"],invert=TRUE),]

files_fund_details_fund_details <- files_fund_details_nonstats[grep("Fund_Details_Fund_Details",files_fund_details_nonstats[,"pull"],invert=FALSE),]
files_fund_details_nonfund_details <- files_fund_details_nonstats[grep("Fund_Details_Fund_Details",files_fund_details_nonstats[,"pull"],invert=TRUE),]

files_fund_details_fee_redemption <- files_fund_details_nonfund_details[grep("Fee_and_Redemption_Structure",files_fund_details_nonfund_details[,"pull"],invert=FALSE),]
files_fund_details_nonfee_redemption <- files_fund_details_nonfund_details[grep("Fee_and_Redemption_Structure",files_fund_details_nonfund_details[,"pull"],invert=TRUE),]

files_fund_details_profile_strategy <- files_fund_details_nonfee_redemption[grep("Profile_Strategy",files_fund_details_nonfee_redemption[,"pull"],invert=FALSE),]
files_fund_details_nonprofile_strategy <- files_fund_details_nonfee_redemption[grep("Profile_Strategy",files_fund_details_nonfee_redemption[,"pull"],invert=TRUE),]

files_fund_details_identifier <- files_fund_details_nonprofile_strategy[grep("Identifier",files_fund_details_nonprofile_strategy[,"pull"],invert=FALSE),]
files_fund_details_nonidentifier <- files_fund_details_nonprofile_strategy[grep("Identifier",files_fund_details_nonprofile_strategy[,"pull"],invert=TRUE),]

# Service Provider, Manager_Details, Industry Focus, & Country Focus
files_fund_details_other <- files_fund_details_nonidentifier

rm2(files_fund_details_all,files_fund_details_out3)
rm2(files_fund_details_nonstats,files_fund_details_nonfund_details,files_fund_details_nonfee_redemption)
rm2(files_fund_details_nonprofile_strategy,files_fund_details_nonidentifier)

#write.csv(files_fund_details_stats[,!(colnames(files_fund_details_stats) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Stats_files",".csv",sep=""),row.names=FALSE)
write.csv(files_fund_details_stats,file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Stats_files",".csv",sep=""),row.names=FALSE)

#write.csv(files_fund_details_fund_details[,!(colnames(files_fund_details_fund_details) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Fund_Details_files",".csv",sep=""),row.names=FALSE)
write.csv(files_fund_details_fund_details,file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Fund_Details_files",".csv",sep=""),row.names=FALSE)

#write.csv(files_fund_details_fee_redemption[,!(colnames(files_fund_details_fee_redemption) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Fee_Redemption_files",".csv",sep=""),row.names=FALSE)
write.csv(files_fund_details_fee_redemption,file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Fee_Redemption_files",".csv",sep=""),row.names=FALSE)

#write.csv(files_fund_details_profile_strategy[,!(colnames(files_fund_details_profile_strategy) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Profile_Strategy_files",".csv",sep=""),row.names=FALSE)
write.csv(files_fund_details_profile_strategy,file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Profile_Strategy_files",".csv",sep=""),row.names=FALSE)

#write.csv(files_fund_details_identifier[,!(colnames(files_fund_details_identifier) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Identifier_files",".csv",sep=""),row.names=FALSE)
write.csv(files_fund_details_identifier,file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Identifier_files",".csv",sep=""),row.names=FALSE)

#write.csv(files_fund_details_other[,!(colnames(files_fund_details_other) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Other_files",".csv",sep=""),row.names=FALSE)
write.csv(files_fund_details_other,file=paste(output_directory,"//","EurekahedgeHF_Fund_Details_Other_files",".csv",sep=""),row.names=FALSE)

rm2(files_fund_details_stats,files_fund_details_fund_details,files_fund_details_fee_redemption)
rm2(files_fund_details_profile_strategy,files_fund_details_identifier,files_fund_details_other)


### Find common cols
fund_details_common_cols <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                       reshape(fund_details_cat_cols_concatenate[,!(colnames(fund_details_cat_cols_concatenate) %in% c("files","yr","month","output_dir","bad_col"))], 
                                               direction="wide",idvar=c("category","col"),timevar=c("file_name")),
                                       Totals=NA,stringsAsFactors=FALSE)

fund_details_common_cols[,"Order_All_Org"] <- seq(1,nrow(fund_details_common_cols),1) 
fund_details_common_cols_num_cols <- colnames(fund_details_common_cols)[!(colnames(fund_details_common_cols) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","category","col","Totals"))]
fund_details_common_cols[,"Totals"] <- rowMeans(fund_details_common_cols[,fund_details_common_cols_num_cols],na.rm=TRUE)
fund_details_common_cols <- fund_details_common_cols[order(fund_details_common_cols[,"Totals"],fund_details_common_cols[,"Order_All_Org"]),]
fund_details_common_cols[,"Order_All_Pos"] <- seq(1,nrow(fund_details_common_cols),1)
fund_details_common_cols[,"Totals"] <- rowSums(!is.na(fund_details_common_cols[,fund_details_common_cols_num_cols]))

fund_details_good_cols <- fund_details_common_cols[fund_details_common_cols[,"Totals"]==length(fund_details_common_cols_num_cols),"col"]
fund_details_bad_cols <- fund_details_common_cols[fund_details_common_cols[,"Totals"]!=length(fund_details_common_cols_num_cols),"col"]

fund_details_common_cols <- rbind(fund_details_common_cols[fund_details_common_cols[,"col"] %in% fund_details_good_cols,],
                                  fund_details_common_cols[fund_details_common_cols[,"col"] %in% fund_details_bad_cols,])
fund_details_common_cols[,"Order_All_Tot"] <- seq(1,nrow(fund_details_common_cols),1)
colnames(fund_details_common_cols) <- gsub(pattern="order.", replacement="", x=colnames(fund_details_common_cols))
row.names(fund_details_common_cols) <- seq(nrow(fund_details_common_cols))

write.csv(fund_details_common_cols, file=paste(common_col_folder_path,"//","EurekahedgeHF_Fund_Details_common_cols_org_melt",".csv",sep=""),row.names=FALSE)

rm2(fund_details_cat_cols_concatenate0,fund_details_cat_cols_concatenate,fund_details_empty_cols)
rm2(fund_details_common_cols,fund_details_good_cols,fund_details_bad_cols)


###############################################################################
cat("SECTION: CLEAN NAV AND AUM FILES", "\n")
###############################################################################

nav_aum_cat_cols_concatenate0 <- adply(.data=files_NAV_AUM, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- files_NAV_AUM[1,]
  # x <- files_NAV_AUM[4,]
  # directory_in <- sep_folder_path
  # unknowns <- unknowns_strings
  
  #Create output directory
  output_folder_path <- paste(x[,"output_dir"], "", sep = "//", collapse = "//")  
  create_directory(output_folder_path,remove=1)
  
  outfile_prefix <- x[,"file_name"]
  
  input <- data.frame(read.csv(file=paste(directory_in,"\\",x[,"files"],sep=""),header=FALSE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input <- input[,colSums(is.na(input))<nrow(input)]
  
  # Get Columns and Category Names
  cat_col <- data.frame(order=NA,col=t(input[1,]),bad_col=NA,stringsAsFactors=FALSE)
  colnames(cat_col) <- c("order","col","bad_col")
  
  cat_col[,"order"] <- seq(nrow(cat_col))
  
  #cat_col <- cat_col[order(cat_col[,"order"]),]
  #row.names(cat_col) <- seq(nrow(cat_col))
  
  cat_col[1,"col"] <- "Ret_AUM"
  
  cat_col[,"bad_col"] <- ifelse(is.na(cat_col[,c("col")]),1,0)
  
  # Clean Column Names
  cat_col[,"col"] <- gsub(pattern="[[:punct:]]", replacement="", x=cat_col[,"col"])
  
  cat_col[,"col"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"col"])
  cat_col[,"col"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"col"])
  cat_col[,"col"] <- gsub("^\\s+|\\s+$", "", cat_col[,"col"])
  
  cat_col[,"col"] <- gsub(pattern=" ", replacement="_", x=cat_col[,"col"])
  
  ## REPLACE COLUMM NAMES
  
  input_trim <- input[2:nrow(input),]
  
  colnames(input_trim) <- paste(cat_col[,"col"],"",sep="")
  
  rm(input)
  
  input_trim[,"Fund_ID"] <- as.integer(input_trim[,"Fund_ID"])
  input_trim <- input_trim[order(input_trim[,"Fund_ID"]),]
  row.names(input_trim) <- seq(nrow(input_trim))
  
  write.csv(input_trim, file=paste(output_folder_path,outfile_prefix,".csv",sep=""),row.names=FALSE)
  
  cat_col_trim <- merge(data.frame(x,cat_id=1,stringsAsFactors=FALSE),
                        data.frame(cat_col,cat_id=1,stringsAsFactors=FALSE), 
                        by.x=c("cat_id"),by.y=c("cat_id"), 
                        all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))
  
  cat_col_trim <- cat_col_trim[order(cat_col_trim[,"order"]),]
  
  #cat_col_trim2 <- cat_col_trim[,!(colnames(cat_col_trim) %in% c("cat_id"))]
  colnames(cat_col_trim)[match("cat_id",names(cat_col_trim))] <- "file_clean"
  cat_col_trim[,"file_clean"] <- paste(output_folder_path,outfile_prefix,".csv",sep="")
  
  row.names(cat_col_trim) <- seq(nrow(cat_col_trim))
  
  rm(input_trim)
  
  rm(output_folder_path,outfile_prefix)
  
  gc()
  
  ##RETURN THE COLNAMES
  return(cat_col_trim)
  
}, directory_in=sep_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

nav_aum_cat_cols_concatenate0 <- nav_aum_cat_cols_concatenate0[,c("files","file_name","file_clean",
                                                                  colnames(nav_aum_cat_cols_concatenate0)[!(colnames(nav_aum_cat_cols_concatenate0) %in% c("files","file_name","file_clean"))])]

### Find empty cols
nav_aum_empty_cols <- nav_aum_cat_cols_concatenate0[nav_aum_cat_cols_concatenate0[,"bad_col"]==1,]

nav_aum_cat_cols_concatenate <- nav_aum_cat_cols_concatenate0[,!(colnames(nav_aum_cat_cols_concatenate0) %in% c("file_clean"))]

write.csv(nav_aum_cat_cols_concatenate, file=paste(common_col_folder_path,"//","EurekahedgeHF_NAV_AUM_common_cols_org",".csv",sep=""),row.names=FALSE)

files_NAV_AUM_out0 <- unique(nav_aum_cat_cols_concatenate0[,c(colnames(files_NAV_AUM),"file_clean")])

colnames(files_NAV_AUM_out0)[match("files",names(files_NAV_AUM_out0))] <- "file_org"
files_NAV_AUM_out0[,"file_org"] <- paste(sep_folder_path,"//",files_NAV_AUM_out0[,"file_org"],sep="")
files_NAV_AUM_out1 <- files_NAV_AUM_out0[!is.na(files_NAV_AUM_out0[,"file_clean"]),]

files_NAV_AUM_out2 <- unique(files_NAV_AUM_out1[,c("yr","month","file_org","file_clean")])
row.names(files_NAV_AUM_out2) <- seq(nrow(files_NAV_AUM_out2))

#write.csv(files_NAV_AUM_out, file=paste(output_directory,"//","EurekahedgeHF_NAV_AUM_files",".csv",sep=""),row.names=FALSE)

files_NAV_AUM_out3 <- data.frame(pull=NA,file_name=NA,files_NAV_AUM_out2,stringsAsFactors=FALSE)
files_NAV_AUM_out3[,"file_name"] <- files_NAV_AUM_out3[,"file_clean"]
files_NAV_AUM_out3[,"file_name"] <- gsub("\\\\","/",files_NAV_AUM_out3[,"file_name"])
files_NAV_AUM_out3[,"file_name"] <- gsub("//","/",files_NAV_AUM_out3[,"file_name"])
files_NAV_AUM_out3[,"file_name"] <- gsub("//","/",files_NAV_AUM_out3[,"file_name"])

files_NAV_AUM_out3[,"file_name"] <- encodeString(files_NAV_AUM_out3[,"file_name"])

#files_NAV_AUM_out3[,"pull"] <- regexpr("/[^/]*$", files_NAV_AUM_out3[,"file_name"])
files_NAV_AUM_out3[,"pull"] <- sapply(gregexpr("\\/", files_NAV_AUM_out3[,"file_name"]), tail, 1)
files_NAV_AUM_out3[,"file_name"] <- substr(files_NAV_AUM_out3[,"file_name"],files_NAV_AUM_out3[,"pull"]+1,nchar(files_NAV_AUM_out3[,"file_name"]))

files_NAV_AUM_out3[,"pull"] <- files_NAV_AUM_out3[,"file_name"]
files_NAV_AUM_out3[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=files_NAV_AUM_out3[,"pull"])

rm2(files_NAV_AUM,files_NAV_AUM_out0,files_NAV_AUM_out1,files_NAV_AUM_out2)

files_NAV_AUM_all <- files_NAV_AUM_out3

files_NAV_AUM_nav_aum <- files_NAV_AUM_all[grep("NAV_AUM",files_NAV_AUM_all[,"pull"],invert=FALSE),]
files_NAV_AUM_nonnav_aum <- files_NAV_AUM_all[grep("NAV_AUM",files_NAV_AUM_all[,"pull"],invert=TRUE),]

rm2(files_NAV_AUM_all,files_NAV_AUM_out3)
rm2(files_NAV_AUM_nonnav_aum)

#write.csv(files_NAV_AUM_nav_aum[,!(colnames(files_NAV_AUM_nav_aum) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_NAV_AUM_files",".csv",sep=""),row.names=FALSE)
write.csv(files_NAV_AUM_nav_aum,file=paste(output_directory,"//","EurekahedgeHF_NAV_AUM_files",".csv",sep=""),row.names=FALSE)

rm2(files_NAV_AUM_nav_aum)


### Find common cols
nav_aum_common_cols <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                  reshape(nav_aum_cat_cols_concatenate[,!(colnames(nav_aum_cat_cols_concatenate) %in% c("files","yr","month","output_dir","bad_col"))], 
                                          direction="wide",idvar=c("col"),timevar=c("file_name")),
                                  Totals=NA,stringsAsFactors=FALSE)

nav_aum_common_cols[,"Order_All_Org"] <- seq(1,nrow(nav_aum_common_cols),1) 
nav_aum_common_cols_num_cols <- colnames(nav_aum_common_cols)[!(colnames(nav_aum_common_cols) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","col","Totals"))]
nav_aum_common_cols[,"Totals"] <- rowMeans(nav_aum_common_cols[,nav_aum_common_cols_num_cols],na.rm=TRUE)
nav_aum_common_cols <- nav_aum_common_cols[order(nav_aum_common_cols[,"Totals"],nav_aum_common_cols[,"Order_All_Org"]),]
nav_aum_common_cols[,"Order_All_Pos"] <- seq(1,nrow(nav_aum_common_cols),1)
nav_aum_common_cols[,"Totals"] <- rowSums(!is.na(nav_aum_common_cols[,nav_aum_common_cols_num_cols]))

nav_aum_good_cols <- nav_aum_common_cols[nav_aum_common_cols[,"Totals"]==length(nav_aum_common_cols_num_cols),"col"]
nav_aum_bad_cols <- nav_aum_common_cols[nav_aum_common_cols[,"Totals"]!=length(nav_aum_common_cols_num_cols),"col"]

nav_aum_common_cols <- rbind(nav_aum_common_cols[nav_aum_common_cols[,"col"] %in% nav_aum_good_cols,],
                             nav_aum_common_cols[nav_aum_common_cols[,"col"] %in% nav_aum_bad_cols,])
nav_aum_common_cols[,"Order_All_Tot"] <- seq(1,nrow(nav_aum_common_cols),1)
colnames(nav_aum_common_cols) <- gsub(pattern="order.", replacement="", x=colnames(nav_aum_common_cols))
row.names(nav_aum_common_cols) <- seq(nrow(nav_aum_common_cols))

write.csv(nav_aum_common_cols, file=paste(common_col_folder_path,"//","EurekahedgeHF_NAV_AUM_common_cols_org_melt",".csv",sep=""),row.names=FALSE)

rm2(nav_aum_cat_cols_concatenate0,nav_aum_cat_cols_concatenate,nav_aum_empty_cols)
rm2(nav_aum_common_cols,nav_aum_good_cols,nav_aum_bad_cols)


###############################################################################
cat("SECTION: CLEAN INSTRUMENTS TRADED FILES", "\n")
###############################################################################

instruments_cat_cols_concatenate0 <- adply(.data=files_instruments, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- files_instruments[1,]
  # x <- files_instruments[9,]
  # directory_in <- sep_folder_path
  # unknowns <- unknowns_strings
  
  #Create output directory
  output_folder_path <- paste(x[,"output_dir"], "", sep = "//", collapse = "//")  
  create_directory(output_folder_path,remove=1)
  
  outfile_prefix <- x[,"file_name"]
  
  input <- data.frame(read.csv(file=paste(directory_in,"\\",x[,"files"],sep=""),header=FALSE,na.strings="NA",stringsAsFactors=FALSE),
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
  
  input <- input[,colSums(is.na(input))<nrow(input)]
  
  #aa <- input[(!is.na(input[,5])),]
  
  # Get Columns and Category Names
  cat_col <- data.frame(order=NA,col=t(input[1,]),bad_col=NA,stringsAsFactors=FALSE)
  colnames(cat_col) <- c("order","col","bad_col")
  
  cat_col[,"order"] <- seq(nrow(cat_col))
  
  #cat_col <- cat_col[order(cat_col[,"order"]),]
  #row.names(cat_col) <- seq(nrow(cat_col))
  
  cat_col[,"bad_col"] <- ifelse(is.na(cat_col[,c("col")]),1,0)
  
  # Clean Column Names
  cat_col[,"col"] <- gsub(pattern="[[:punct:]]", replacement="", x=cat_col[,"col"])
  
  cat_col[,"col"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"col"])
  cat_col[,"col"] <- gsub(pattern=" {2,}", replacement=" ", x=cat_col[,"col"])
  cat_col[,"col"] <- gsub("^\\s+|\\s+$", "", cat_col[,"col"])
  
  cat_col[,"col"] <- gsub(pattern=" ", replacement="_", x=cat_col[,"col"])
  
  ## REPLACE COLUMM NAMES
  
  input_trim <- input[2:nrow(input),]
  
  colnames(input_trim) <- paste(cat_col[,"col"],"",sep="")
  
  rm(input)
  
  input_trim[,"Fund_ID"] <- as.integer(input_trim[,"Fund_ID"])
  input_trim <- input_trim[order(input_trim[,"Fund_ID"]),]
  row.names(input_trim) <- seq(nrow(input_trim))
  
  write.csv(input_trim, file=paste(output_folder_path,outfile_prefix,".csv",sep=""),row.names=FALSE)
  
  cat_col_trim <- merge(data.frame(x,cat_id=1,stringsAsFactors=FALSE),
                        data.frame(cat_col,cat_id=1,stringsAsFactors=FALSE), 
                        by.x=c("cat_id"),by.y=c("cat_id"), 
                        all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))
  
  cat_col_trim <- cat_col_trim[order(cat_col_trim[,"order"]),]
  
  #cat_col_trim2 <- cat_col_trim[,!(colnames(cat_col_trim) %in% c("cat_id"))]
  colnames(cat_col_trim)[match("cat_id",names(cat_col_trim))] <- "file_clean"
  cat_col_trim[,"file_clean"] <- paste(output_folder_path,outfile_prefix,".csv",sep="")
  
  row.names(cat_col_trim) <- seq(nrow(cat_col_trim))
  
  rm(input_trim)
  
  rm(output_folder_path,outfile_prefix)
  
  gc()
  
  ##RETURN THE COLNAMES
  return(cat_col_trim)
  
  
}, directory_in=sep_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

instruments_cat_cols_concatenate0 <- instruments_cat_cols_concatenate0[,c("files","file_name","file_clean",
                                                                          colnames(instruments_cat_cols_concatenate0)[!(colnames(instruments_cat_cols_concatenate0) %in% c("files","file_name","file_clean"))])]

### Find empty cols
instruments_empty_cols <- instruments_cat_cols_concatenate0[instruments_cat_cols_concatenate0[,"bad_col"]==1,]

instruments_cat_cols_concatenate <- instruments_cat_cols_concatenate0[,!(colnames(instruments_cat_cols_concatenate0) %in% c("file_clean"))]

write.csv(instruments_cat_cols_concatenate, file=paste(common_col_folder_path,"//","EurekahedgeHF_Instruments_Traded_common_cols_org",".csv",sep=""),row.names=FALSE)

files_instruments_out0 <- unique(instruments_cat_cols_concatenate0[,c(colnames(files_instruments),"file_clean")])

colnames(files_instruments_out0)[match("files",names(files_instruments_out0))] <- "file_org"
files_instruments_out0[,"file_org"] <- paste(sep_folder_path,"//",files_instruments_out0[,"file_org"],sep="")
files_instruments_out1 <- files_instruments_out0[!is.na(files_instruments_out0[,"file_clean"]),]

files_instruments_out2 <- unique(files_instruments_out1[,c("yr","month","file_org","file_clean")])
row.names(files_instruments_out2) <- seq(nrow(files_instruments_out2))

#write.csv(files_instruments_out2, file=paste(output_directory,"//","EurekahedgeHF_Instruments_files",".csv",sep=""),row.names=FALSE)

files_instruments_out3 <- data.frame(pull=NA,file_name=NA,files_instruments_out2,stringsAsFactors=FALSE)
files_instruments_out3[,"file_name"] <- files_instruments_out3[,"file_clean"]
files_instruments_out3[,"file_name"] <- gsub("\\\\","/",files_instruments_out3[,"file_name"])
files_instruments_out3[,"file_name"] <- gsub("//","/",files_instruments_out3[,"file_name"])
files_instruments_out3[,"file_name"] <- gsub("//","/",files_instruments_out3[,"file_name"])

files_instruments_out3[,"file_name"] <- encodeString(files_instruments_out3[,"file_name"])

#files_instruments_out3[,"pull"] <- regexpr("/[^/]*$", files_instruments_out3[,"file_name"])
files_instruments_out3[,"pull"] <- sapply(gregexpr("\\/", files_instruments_out3[,"file_name"]), tail, 1)
files_instruments_out3[,"file_name"] <- substr(files_instruments_out3[,"file_name"],files_instruments_out3[,"pull"]+1,nchar(files_instruments_out3[,"file_name"]))

files_instruments_out3[,"pull"] <- files_instruments_out3[,"file_name"]
files_instruments_out3[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=files_instruments_out3[,"pull"])

rm2(files_instruments,files_instruments_out0,files_instruments_out1,files_instruments_out2)

files_instruments_all <- files_instruments_out3

files_instruments_instruments <- files_instruments_all[grep("Instrument",files_instruments_all[,"pull"],invert=FALSE),]
files_instruments_noninstruments <- files_instruments_all[grep("Instrument",files_instruments_all[,"pull"],invert=TRUE),]

rm2(files_instruments_all,files_instruments_out3)
rm2(files_instruments_noninstruments)

#write.csv(files_instruments_instruments[,!(colnames(files_instruments_instruments) %in% c("pull","file_name"))],
#          file=paste(output_directory,"//","EurekahedgeHF_Instruments_files",".csv",sep=""),row.names=FALSE)
write.csv(files_instruments_instruments,file=paste(output_directory,"//","EurekahedgeHF_Instruments_files",".csv",sep=""),row.names=FALSE)

rm2(files_instruments_instruments)


### Find common cols
instruments_common_cols <- data.frame(Order_All_Org=NA,Order_All_Pos=NA,Order_All_Tot=NA,
                                      reshape(instruments_cat_cols_concatenate[,!(colnames(instruments_cat_cols_concatenate) %in% c("files","yr","month","output_dir","bad_col"))], 
                                              direction="wide",idvar=c("col"),timevar=c("file_name")),
                                      Totals=NA,stringsAsFactors=FALSE)

instruments_common_cols[,"Order_All_Org"] <- seq(1,nrow(instruments_common_cols),1) 
instruments_common_cols_num_cols <- colnames(instruments_common_cols)[!(colnames(instruments_common_cols) %in% c("Order_All_Org","Order_All_Pos","Order_All_Tot","col","Totals"))]
instruments_common_cols[,"Totals"] <- rowMeans(instruments_common_cols[,instruments_common_cols_num_cols],na.rm=TRUE)
instruments_common_cols <- instruments_common_cols[order(instruments_common_cols[,"Totals"],instruments_common_cols[,"Order_All_Org"]),]
instruments_common_cols[,"Order_All_Pos"] <- seq(1,nrow(instruments_common_cols),1)
instruments_common_cols[,"Totals"] <- rowSums(!is.na(instruments_common_cols[,instruments_common_cols_num_cols]))

instruments_good_cols <- instruments_common_cols[instruments_common_cols[,"Totals"]==length(instruments_common_cols_num_cols),"col"]
instruments_bad_cols <- instruments_common_cols[instruments_common_cols[,"Totals"]!=length(instruments_common_cols_num_cols),"col"]

instruments_common_cols <- rbind(instruments_common_cols[instruments_common_cols[,"col"] %in% instruments_good_cols,],
                                 instruments_common_cols[instruments_common_cols[,"col"] %in% instruments_bad_cols,])
instruments_common_cols[,"Order_All_Tot"] <- seq(1,nrow(instruments_common_cols),1)
colnames(instruments_common_cols) <- gsub(pattern="order.", replacement="", x=colnames(instruments_common_cols))
row.names(instruments_common_cols) <- seq(nrow(instruments_common_cols))

write.csv(instruments_common_cols, file=paste(common_col_folder_path,"//","EurekahedgeHF_Instruments_Traded_common_cols_org_melt",".csv",sep=""),row.names=FALSE)

rm2(instruments_cat_cols_concatenate0,instruments_cat_cols_concatenate,instruments_empty_cols)
rm2(instruments_common_cols,instruments_good_cols,instruments_bad_cols)

rm2(sep_folder_path,common_col_folder_path)
