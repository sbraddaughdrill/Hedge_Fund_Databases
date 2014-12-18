# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Expand_Fields_part1.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Expand Eurekahedge Fields
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

external_packages <- c("cwhmisc","data.table","DataCombine","formatR","gdata","gtools","Hmisc","lubridate","memisc",
                       "mitools","plm","plyr","reshape2","splitstackshape","stringi","stringr","taRifx","tm","zoo")
invisible(unlist(sapply(external_packages,load_external_packages, repo_str=repo, simplify=FALSE, USE.NAMES=FALSE)))
installed_packages <- list_installed_packages(external_packages)

rm(external_packages,installed_packages,repo)


###############################################################################
cat("SECTION: DEFINE DIRECTORIES", "\n")
###############################################################################

#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

#Check to see if final folder exists.  If not, create it.
final_folder_expand_path <- paste(output_directory, "Final_Expand1", sep = "//", collapse = "//")  
create_directory(final_folder_expand_path,remove=1)

final_folder_files0 <- data.frame(files=list.files(path=final_folder_path),file_name=NA,import=NA,stringsAsFactors=FALSE)

#final_folder_files <- final_folder_files0[!grepl(".TXT|.txt", final_folder_files0[,"files"]),]
final_folder_files <- final_folder_files0[grepl(".CSV|.csv", final_folder_files0[,"files"]),]

rm2(final_folder_files0)

final_folder_files[,"file_name"] <- final_folder_files[,"files"]
final_folder_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=final_folder_files[,"file_name"])

final_folder_files[,"import"] <- ifelse(grepl("(Stats|Fund_Detail|Fee_and_Redemption|Profile_Strategy|Identifier)",final_folder_files[,"file_name"]),1,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(NAV_AUM_Ret)",final_folder_files[,"file_name"]),2,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(Other)",final_folder_files[,"file_name"]),3,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(Instruments_Traded)",final_folder_files[,"file_name"]),4,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(is.na(final_folder_files[,"import"]),0,final_folder_files[,"import"])

invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: IMPORT FILES", "\n")
###############################################################################

a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(1),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[2,]
  # x <- final_folder_files[3,]
  # x <- final_folder_files[4,] 
  
  # directory_in <- final_folder_path
  # unknowns <- unknowns_strings
  
  #temp_cols <- c("date","yr","month","bad_min","bad_max")
  
  input <- data.table(pull_trim=NA,
                      #matrix(NA, ncol=length(temp_cols), nrow=1, dimnames=list(c(), temp_cols)),
                      read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE))
  
  #set(input, j=which(colnames(input) %in% temp_cols), value=NULL )
  #rm(temp_cols)
  
  setkeyv(input,NULL)
  setnames(input,"pull","pull_trim2")
  
  order_ids <- c("pull_trim","pull_trim2","pull")
  order_ids_trim <- order_ids[order_ids %in% colnames(input)]
  setcolorder(input, c(order_ids_trim,colnames(input)[!(colnames(input) %in% c(order_ids_trim))]))
  
  rm(order_ids,order_ids_trim)
  
  input <- input[, pull_trim:=pull_trim2,by=NULL]
  
  cols <- c("pull_trim")
  for (k in cols) {
    set(input, i=NULL, j=k, value=gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", input[[k]], perl=TRUE))
    set(input, i=NULL, j=k, value=as.character(input[[k]]))
    
    #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
  }
  rm(k,cols)
  
  cols <- c("pull_trim2")
  for (k in cols) {
    set(input, i=NULL, j=k, value=gsub(pattern="_Fund_Details_Fund_Details", replacement="", input[[k]], perl=TRUE))
    set(input, i=NULL, j=k, value=gsub(pattern="_Fund_Details_Fee_and_Redemption_Structure", replacement="", input[[k]], perl=TRUE))
    set(input, i=NULL, j=k, value=gsub(pattern="_Fund_Details_Unique_Identifiers", replacement="", input[[k]], perl=TRUE))
    set(input, i=NULL, j=k, value=gsub(pattern="_Fund_Details_Profile_Strategy_Description", replacement="", input[[k]], perl=TRUE))
    set(input, i=NULL, j=k, value=gsub(pattern="_Fund_Details_Statistics", replacement="", input[[k]], perl=TRUE))
    set(input, i=NULL, j=k, value=gsub(pattern="_Instruments_Traded", replacement="", input[[k]], perl=TRUE))
    
    #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
  }
  rm(k,cols)
  
  for (k in which(sapply(input,class)=="character")) 
  {
    set(input, i=NULL, j=k, value=gsub("^\\s+|\\s+$", "", input[[k]], perl=TRUE))
  }
  rm(k)
  
  for (k in colnames(input)) 
  {
    #k <- 1
    set(input, i=NULL, j=k, value=unknownToNA(input[[k]], unknown=unknowns,force=TRUE))
    set(input, i=NULL, j=k, value=ifelse(is.na(input[[k]]),NA,input[[k]]))
  }
  rm(k)
  
  input[,droprow := rowSums(is.na(.SD))<ncol(input),.SDcols=colnames(input)]
  #set(input,i=NULL,j=which(colnames(input)==c("droprow")),value=NULL)
  input <- input[(droprow)][,droprow:=NULL][]
  
  input <- input[,which(unlist(lapply(input, function(x)!all(is.na(x))))),with=FALSE]
  
  setorderv(input, c("Fund_ID","pull_trim","pull_trim2"),c(1,1,1))
  
  if("Date_Added" %in% colnames(input)){input <- input[, Date_Added:=as.yearmon(Date_Added,format="%b %Y")]} 
  if("Dead_Date" %in% colnames(input)){input <- input[, Dead_Date:=as.yearmon(Dead_Date,format="%b %Y")]} 
  if("Inception_Date" %in% colnames(input)){input <- input[, Inception_Date:=as.yearmon(Inception_Date,format="%b %Y")]} 
  
  assign(x[,"file_name"], input, envir = .GlobalEnv)
  
  rm(input)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
}, directory_in=final_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))

a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(3),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[6,]
  
  # directory_in <- final_folder_path
  # unknowns <- unknowns_strings
  
  input <- data.table(read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE))
  setkeyv(input,NULL)
  
  cols <- c("pull_trim")
  for (k in cols) {
    set(input, i=NULL, j=k, value=as.character(input[[k]]))
  }
  rm(k,cols)
  
  for (k in which(sapply(input,class)=="character")) 
  {
    set(input, i=NULL, j=k, value=gsub("^\\s+|\\s+$", "", input[[k]], perl=TRUE))
  }
  rm(k)
  
  for (k in colnames(input)) 
  {
    #k <- 1
    set(input, i=NULL, j=k, value=unknownToNA(input[[k]], unknown=unknowns,force=TRUE))
    set(input, i=NULL, j=k, value=ifelse(is.na(input[[k]]),NA,input[[k]]))
  }
  rm(k)
  
  input[,droprow := rowSums(is.na(.SD))<ncol(input),.SDcols=colnames(input)]
  #set(input,i=NULL,j=which(colnames(input)==c("droprow")),value=NULL)
  input <- input[(droprow)][,droprow:=NULL][]
  
  input <- input[,which(unlist(lapply(input, function(x)!all(is.na(x))))),with=FALSE]
  
  setorderv(input, c("Fund_ID","pull_trim","pull_trim2"),c(1,1,1))
  
  if("Date_Added" %in% colnames(input)){input <- input[, Date_Added:=as.yearmon(Date_Added,format="%b %Y")]} 
  if("Dead_Date" %in% colnames(input)){input <- input[, Dead_Date:=as.yearmon(Dead_Date,format="%b %Y")]} 
  if("Inception_Date" %in% colnames(input)){input <- input[, Inception_Date:=as.yearmon(Inception_Date,format="%b %Y")]} 
  
  assign(x[,"file_name"], input, envir = .GlobalEnv)
  
  rm(input)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
}, directory_in=final_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))


# ###############################################################################
# cat("SECTION: IMPORT DEAD DATES", "\n")
# ###############################################################################
# 
# #Import Dead Dates
# 
# Dead_Dates <- read.csv(file=paste(final_folder_path,"//","EurekahedgeHF_Dead_Dates_fixed",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)
# 
# for(i in which(sapply(Dead_Dates,class)=="character"))
# {
#   Dead_Dates[[i]] = trim(Dead_Dates[[i]])
# }
# rm(i)
# for (i in 1:ncol(Dead_Dates))
# {
#   Dead_Dates[,i] <- unknownToNA(Dead_Dates[,i], unknown=unknowns_strings,force=TRUE)
#   Dead_Dates[,i] <- ifelse(is.na(Dead_Dates[,i]),NA,Dead_Dates[,i])
# } 
# rm(i)
# 
# Dead_Dates <- Dead_Dates[rowSums(is.na(Dead_Dates[,1:ncol(Dead_Dates)]))<ncol(Dead_Dates),]
# row.names(Dead_Dates) <- seq(nrow(Dead_Dates))
# 
# Dead_Dates <- Dead_Dates[,colSums(is.na(Dead_Dates))<nrow(Dead_Dates)]
# 
# Dead_Dates  <- Dead_Dates[order(Dead_Dates[,"Fund_ID"],Dead_Dates[,"Fund_Name"],Dead_Dates[,"pull_trim"],Dead_Dates[,"pull"]),]
# row.names(Dead_Dates) <- seq(nrow(Dead_Dates))
# 
# Dead_Dates[,"Dead_Date"] <- as.Date(Dead_Dates[,"Dead_Date"],format="%Y-%m-%d")
# Dead_Dates[,"Dead_Date"] <- as.yearmon(Dead_Dates[,"Dead_Date"],format="%b %Y")
# 
# 
# ###############################################################################
# cat("SECTION: FIX DEAD DATES", "\n")
# ###############################################################################
# 
# #Fix Dates
# 
# EurekahedgeHF_Stats_noreturns_deaddates <- EurekahedgeHF_Stats_noreturns
# EurekahedgeHF_Fund_Details_deaddates <- EurekahedgeHF_Fund_Details
# EurekahedgeHF_Fee_and_Redemption_deaddates <- EurekahedgeHF_Fee_and_Redemption
# EurekahedgeHF_Profile_Strategy_deaddates <- EurekahedgeHF_Profile_Strategy
# EurekahedgeHF_Identifiers_deaddates <- EurekahedgeHF_Identifiers
# EurekahedgeHF_Other_deaddates <- EurekahedgeHF_Other
# 
# rm2(EurekahedgeHF_Stats_noreturns,EurekahedgeHF_Fund_Details,EurekahedgeHF_Fee_and_Redemption)
# rm2(EurekahedgeHF_Profile_Strategy,EurekahedgeHF_Identifiers,EurekahedgeHF_Other)
# 
# fix_dead_dates_temp1 <- list(data=c("EurekahedgeHF_Stats_noreturns_deaddates"),col=c("Dead_Date"))
# fix_dead_dates_temp2 <- list(data=c("EurekahedgeHF_Fund_Details_deaddates"),col=c("Dead_Date"))
# fix_dead_dates_temp3 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_deaddates"),col=c("Dead_Date"))
# fix_dead_dates_temp4 <- list(data=c("EurekahedgeHF_Profile_Strategy_deaddates"),col=c("Dead_Date"))
# fix_dead_dates_temp5 <- list(data=c("EurekahedgeHF_Identifiers_deaddates"),col=c("Dead_Date"))
# fix_dead_dates_temp6 <- list(data=c("EurekahedgeHF_Other_deaddates"),col=c("Dead_Date"))
# 
# fix_dead_dates_all <- list(fix_dead_dates_temp1,fix_dead_dates_temp2,fix_dead_dates_temp3,
#                            fix_dead_dates_temp4,fix_dead_dates_temp5,fix_dead_dates_temp6)
# 
# rm2(fix_dead_dates_temp1,fix_dead_dates_temp2,fix_dead_dates_temp3)
# rm2(fix_dead_dates_temp4,fix_dead_dates_temp5,fix_dead_dates_temp6)
# 
# l_ply(.data=fix_dead_dates_all, .fun = function(x,dead_dates_fixed){
#   
#   # x <- fix_dead_dates_all[[1]]
#   # dead_dates_fixed <- Dead_Dates
#   
#   dead_dates_fixed_trim <- unique(dead_dates_fixed[,!(colnames(dead_dates_fixed) %in% c("pull"))])
#   
#   data_temp <- get(x[[1]])
#   col_temp <- x[[2]]
#   
#   col_order <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
#   col_order[,"order"] <- seq(1,nrow(col_order))
#   
#   col_order_org <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
#   col_order_org[,"order"] <- seq(1,nrow(col_order_org))
#   col_order_org[,"col"] <- ifelse(col_order_org[,"col"] %in% col_temp,paste(col_order_org[,"col"],"_org",sep=""),col_order_org[,"col"])
#   
#   col_order_all <- unique(rbind(col_order_org,col_order))
#   
#   col_order_all <- col_order_all[order(col_order_all[,"order"]),]
#   row.names(col_order_all) <- seq(nrow(col_order_all))
#   
#   #Rename original columns
#   data_temp <- rename.vars(data_temp, col_temp, paste(col_temp,"_org",sep=""),info=FALSE)
#   
#   data_temp2 <- merge(data_temp,dead_dates_fixed_trim,
#                       by.x=c("Fund_ID","Fund_Name"), 
#                       by.y=c("Fund_ID","Fund_Name"), 
#                       all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))
#   
#   data_temp2  <- data_temp2[order(data_temp2[,"Fund_ID"],data_temp2[,"Fund_Name"],data_temp2[,"pull"]),]
#   row.names(data_temp2) <- seq(nrow(data_temp2))
#   
#   #Get column order
#   
#   #data_temp2 <- data_temp2[,sort(colnames(data_temp2), decreasing = FALSE)]
#   data_temp2 <- data_temp2[,col_order_all[,"col"]]
#   
#   assign(x[[1]], data_temp2, envir = .GlobalEnv)
#   
#   gc()
#   
# },dead_dates_fixed=Dead_Dates, .progress = "text")
# 
# rm2(fix_dead_dates_all)


###############################################################################
cat("SECTION: STRIP COMMENTS FROM VARIABLES", "\n")
###############################################################################

strip_comments_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns"),data_out=c("EurekahedgeHF_Stats_noreturns_comments"),
                             col=NULL)
strip_comments_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details"),data_out=c("EurekahedgeHF_Fund_Details_comments"),
                             col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange"))
strip_comments_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption"),data_out=c("EurekahedgeHF_Fee_and_Redemption_comments"),
                             col=c("Management_Fee","Performance_Fee","Other_Fee"))
strip_comments_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy"),data_out=c("EurekahedgeHF_Profile_Strategy_comments"),
                             col=NULL)
strip_comments_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers"),data_out=c("EurekahedgeHF_Identifiers_comments"),
                             col=NULL)
strip_comments_temp6 <- list(data_in=c("EurekahedgeHF_Other"),data_out=c("EurekahedgeHF_Other_comments"),
                             col=NULL)

strip_comments_all0 <- list(strip_comments_temp1,strip_comments_temp2,strip_comments_temp3,
                            strip_comments_temp4,strip_comments_temp5,strip_comments_temp6)

rm2(strip_comments_temp1,strip_comments_temp2,strip_comments_temp3,strip_comments_temp4,strip_comments_temp5,strip_comments_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

strip_comments_all1 <- llply(.data=strip_comments_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(strip_comments_all0)

strip_comments_all1 <- strip_comments_all1[order(-sapply(strip_comments_all1,"[[","size"),
                                                 -sapply(strip_comments_all1,"[[","ncol"),
                                                 -sapply(strip_comments_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=strip_comments_all1, .fun = function(x){
  
  # x <- strip_comments_all1[[1]]
  # x <- strip_comments_all1[[2]]
  # x <- strip_comments_all1[[3]]
  # x <- strip_comments_all1[[4]]
  # x <- strip_comments_all1[[5]]
  # x <- strip_comments_all1[[6]]
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    col_temp <- x[[which(names(x)==c("col"))]]
    
    col_order <- data.frame(col=colnames(get(x[[which(names(x)==c("data_in"))]])),order=NA,stringsAsFactors=FALSE)
    col_order[,"order"] <- seq(1,nrow(col_order))
    
    col_order_org <- data.frame(col=colnames(get(x[[which(names(x)==c("data_in"))]])),order=NA,stringsAsFactors=FALSE)
    col_order_org[,"order"] <- seq(1,nrow(col_order_org))
    col_order_org[,"col"] <- ifelse(col_order_org[,"col"] %in% col_temp, paste(col_order_org[,"col"],"_org",sep=""),col_order_org[,"col"])
    
    col_order_comments <- data.frame(col=colnames(get(x[[which(names(x)==c("data_in"))]])),order=NA,stringsAsFactors=FALSE)
    col_order_comments[,"order"] <- seq(1,nrow(col_order_comments))
    col_order_comments[,"col"] <- ifelse(col_order_comments[,"col"] %in% col_temp,paste(col_order_comments[,"col"],"_comments",sep=""),col_order_comments[,"col"])
    
    col_order_all <- unique(rbind(rbind(col_order_org,col_order),col_order_comments))
    
    col_order_all <- col_order_all[order(col_order_all[,"order"]),]
    row.names(col_order_all) <- seq(nrow(col_order_all))
    
    strip_cols <- c(paste(col_temp,"_temp",sep=""), paste(col_temp,"_comments",sep=""))
    
    data_temp2 <-  data.table(get(x[[which(names(x)==c("data_in"))]]), matrix(NA, ncol=length(strip_cols), nrow=1, dimnames=list(c(), strip_cols)))
    setkeyv(data_temp2,NULL)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    setnames(data_temp2, old=col_temp, new=paste(col_temp,"_org",sep=""))
    setnames(data_temp2, old=paste(col_temp,"_temp",sep=""), new=col_temp)
    
    setcolorder(data_temp2, col_order_all[,"col"])
    
    data_temp2 <- strip_comments(data_temp2,col_temp)
    
    data_temp2 <- create_noncomments(data_temp2,col_temp)
    
    assign(x[[which(names(x)==c("data_out"))]], data_temp2, envir = .GlobalEnv)
    
    rm(data_temp2)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
}, .progress = "text")

rm2(strip_comments_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: CHECK FOR UNKNOWNS", "\n")
###############################################################################

#EurekahedgeHF_Stats_noreturns_unknowns <- EurekahedgeHF_Stats_noreturns_comments
#EurekahedgeHF_Fund_Details_unknowns <- EurekahedgeHF_Fund_Details_comments
#EurekahedgeHF_Fee_and_Redemption_unknowns <- EurekahedgeHF_Fee_and_Redemption_comments
#EurekahedgeHF_Profile_Strategy_unknowns <- EurekahedgeHF_Profile_Strategy_comments
#EurekahedgeHF_Identifiers_unknowns <- EurekahedgeHF_Identifiers_comments
#EurekahedgeHF_Other_unknowns <- EurekahedgeHF_Other_comments
#EurekahedgeHF_Instruments_Traded_unknowns <- EurekahedgeHF_Instruments_Traded_comments

#rm2(EurekahedgeHF_Stats_noreturns_comments,EurekahedgeHF_Fund_Details_comments,EurekahedgeHF_Fee_and_Redemption_comments)
#rm2(EurekahedgeHF_Profile_Strategy_comments,EurekahedgeHF_Identifiers_comments,EurekahedgeHF_Other_comments)
#rm2(EurekahedgeHF_Instruments_Traded_comments)

#unknowns_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_unknowns"),
#                       col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange","Exchange_Name"))
#unknowns_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_unknowns"),
#                       col=c("Management_Fee","Performance_Fee","Other_Fee"))

#unknowns_all <- list(unknowns_temp1,unknowns_temp2)

#rm2(unknowns_temp1,unknowns_temp2)
#invisible(gc(verbose = FALSE, reset = TRUE))

unknowns_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_comments"),data_out=c("EurekahedgeHF_Stats_noreturns_unknowns"),
                       col=NULL)
unknowns_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_comments"),data_out=c("EurekahedgeHF_Fund_Details_unknowns"),
                       col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange","Exchange_Name"))
unknowns_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_comments"),data_out=c("EurekahedgeHF_Fee_and_Redemption_unknowns"),
                       col=c("Management_Fee","Performance_Fee","Other_Fee"))
unknowns_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_comments"),data_out=c("EurekahedgeHF_Profile_Strategy_unknowns"),
                       col=NULL)
unknowns_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_comments"),data_out=c("EurekahedgeHF_Identifiers_unknowns"),
                       col=NULL)
unknowns_temp6 <- list(data_in=c("EurekahedgeHF_Other_comments"),data_out=c("EurekahedgeHF_Other_unknowns"),
                       col=NULL)

unknowns_all0 <- list(unknowns_temp1,unknowns_temp2,unknowns_temp3,
                      unknowns_temp4,unknowns_temp5,unknowns_temp6)

rm2(unknowns_temp1,unknowns_temp2,unknowns_temp3,unknowns_temp4,unknowns_temp5,unknowns_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

unknowns_all1 <- llply(.data=unknowns_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(unknowns_all0)

unknowns_all1 <- unknowns_all1[order(-sapply(unknowns_all1,"[[","size"),
                                     -sapply(unknowns_all1,"[[","ncol"),
                                     -sapply(unknowns_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=unknowns_all1, .fun = function(x,unknowns){
  
  # x <- unknowns_all1[[1]]
  # unknowns <- unknowns_strings
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    data_temp <- data.table(get(x[[which(names(x)==c("data_in"))]]))
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    col_temp <- x[[which(names(x)==c("col"))]]
    
    for (k in col_temp) 
    {
      #k <- 1
      set(data_temp, i=NULL, j=k, value=unknownToNA(data_temp[[k]], unknown=unknowns,force=TRUE))
      set(data_temp, i=NULL, j=k, value=ifelse(is.na(data_temp[[k]]),NA,data_temp[[k]]))
    }
    rm(k)
    
    assign(x[[which(names(x)==c("data_out"))]], data_temp, envir = .GlobalEnv)
    
    rm(data_temp,col_temp)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
  #   
  #   data_temp <- data.table(get(x[[1]]))
  #   col_temp <- x[[2]]
  #   
  #   for (k in col_temp) 
  #   {
  #     #k <- 1
  #     set(data_temp, i=NULL, j=k, value=unknownToNA(data_temp[[k]], unknown=unknowns,force=TRUE))
  #     set(data_temp, i=NULL, j=k, value=ifelse(is.na(data_temp[[k]]),NA,data_temp[[k]]))
  #   }
  #   rm(k)
  #   
  #   assign(x[[1]], data_temp, envir = .GlobalEnv)
  #   
  #   rm(data_temp,col_temp)
  #   
  #   invisible(gc(verbose = FALSE, reset = TRUE))
  
},unknowns=unknowns_strings, .progress = "text")

rm2(unknowns_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: CHANGE NOT SPECIFIED PHRASES TO NA", "\n")
###############################################################################

NA_Phrases <- c("NA","N/A","N\\A","NOT APPLICABLE","NOT APPILCABLE","NOT DEFINED","NOT DISCLOSED","NOT DISLCOSED","NOT DISLOSED","UNDISCLOSED",
                "TO BE ADVISED","TO BE ADVISE","TBA","SEE PROSPECTUS FOR FULL DETAILS","UPON REQUEST",
                "SUBJECT TO MANAGER'S DISCRETION")

# EurekahedgeHF_Stats_noreturns_na_phrases <- EurekahedgeHF_Stats_noreturns_unknowns
# EurekahedgeHF_Fund_Details_na_phrases <- EurekahedgeHF_Fund_Details_unknowns
# EurekahedgeHF_Fee_and_Redemption_na_phrases <- EurekahedgeHF_Fee_and_Redemption_unknowns
# EurekahedgeHF_Profile_Strategy_na_phrases <- EurekahedgeHF_Profile_Strategy_unknowns
# EurekahedgeHF_Identifiers_na_phrases <- EurekahedgeHF_Identifiers_unknowns
# EurekahedgeHF_Other_na_phrases <- EurekahedgeHF_Other_unknowns
# #EurekahedgeHF_Instruments_Traded_na_phrases <- EurekahedgeHF_Instruments_Traded_unknowns
# 
# rm2(EurekahedgeHF_Stats_noreturns_unknowns,EurekahedgeHF_Fund_Details_unknowns,EurekahedgeHF_Fee_and_Redemption_unknowns)
# rm2(EurekahedgeHF_Profile_Strategy_unknowns,EurekahedgeHF_Identifiers_unknowns,EurekahedgeHF_Other_unknowns)
# #rm2(EurekahedgeHF_Instruments_Traded_unknowns)
# 
# na_phrases_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_na_phrases"),
#                          col=c("Dividend_Policy","Domicile","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange","Exchange_Name",
#                                "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
# na_phrases_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_na_phrases"),
#                          col=c("Management_Fee","Performance_Fee","Other_Fee"))
# 
# na_phrases_all <- list(na_phrases_temp1,na_phrases_temp2)
# 
# rm2(na_phrases_temp1,na_phrases_temp2)
# invisible(gc(verbose = FALSE, reset = TRUE))

na_phrases_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_unknowns"),data_out=c("EurekahedgeHF_Stats_noreturns_na_phrases"),
                         col=NULL)
na_phrases_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_unknowns"),data_out=c("EurekahedgeHF_Fund_Details_na_phrases"),
                         col=c("Dividend_Policy","Domicile","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange","Exchange_Name",
                               "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
na_phrases_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_unknowns"),data_out=c("EurekahedgeHF_Fee_and_Redemption_na_phrases"),
                         col=c("Management_Fee","Performance_Fee","Other_Fee"))
na_phrases_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_unknowns"),data_out=c("EurekahedgeHF_Profile_Strategy_na_phrases"),
                         col=NULL)
na_phrases_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_unknowns"),data_out=c("EurekahedgeHF_Identifiers_na_phrases"),
                         col=NULL)
na_phrases_temp6 <- list(data_in=c("EurekahedgeHF_Other_unknowns"),data_out=c("EurekahedgeHF_Other_na_phrases"),
                         col=NULL)

na_phrases_all0 <- list(na_phrases_temp1,na_phrases_temp2,na_phrases_temp3,
                        na_phrases_temp4,na_phrases_temp5,na_phrases_temp6)

rm2(na_phrases_temp1,na_phrases_temp2,na_phrases_temp3,na_phrases_temp4,na_phrases_temp5,na_phrases_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

na_phrases_all1 <- llply(.data=na_phrases_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(na_phrases_all0)

na_phrases_all1 <- na_phrases_all1[order(-sapply(na_phrases_all1,"[[","size"),
                                         -sapply(na_phrases_all1,"[[","ncol"),
                                         -sapply(na_phrases_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=na_phrases_all1, .fun = function(x,phrases){
  
  # x <- na_phrases_all1[[1]]
  # phrases <- NA_Phrases
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], 
           not_specified_to_na(data.table(get(x[[which(names(x)==c("data_in"))]])),x[[which(names(x)==c("col"))]],phrases), 
           envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
  #assign(x[[1]], not_specified_to_na(data.table(get(x[[1]])),x[[2]],phrases), envir = .GlobalEnv)
  
  #invisible(gc(verbose = FALSE, reset = TRUE))
  
},phrases=NA_Phrases, .progress = "text")

rm2(na_phrases_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: CHANGE NO PHRASES TO 'NO'", "\n")
###############################################################################

NO_Phrases <- c("NIL","NONE","NONE AFTER 12 MONTHS","NONE AFTER 1ST YEAR","NO DIVIDEND","NON DIVIDEND","LITTLE OR NO")

# EurekahedgeHF_Stats_noreturns_no_phrases <- EurekahedgeHF_Stats_noreturns_na_phrases
# EurekahedgeHF_Fund_Details_no_phrases <- EurekahedgeHF_Fund_Details_na_phrases
# EurekahedgeHF_Fee_and_Redemption_no_phrases <- EurekahedgeHF_Fee_and_Redemption_na_phrases
# EurekahedgeHF_Profile_Strategy_no_phrases <- EurekahedgeHF_Profile_Strategy_na_phrases
# EurekahedgeHF_Identifiers_no_phrases <- EurekahedgeHF_Identifiers_na_phrases
# EurekahedgeHF_Other_no_phrases <- EurekahedgeHF_Other_na_phrases
# #EurekahedgeHF_Instruments_Traded_no_phrases <- EurekahedgeHF_Instruments_Traded_na_phrases
# 
# rm2(EurekahedgeHF_Stats_noreturns_na_phrases,EurekahedgeHF_Fund_Details_na_phrases,EurekahedgeHF_Fee_and_Redemption_na_phrases)
# rm2(EurekahedgeHF_Profile_Strategy_na_phrases,EurekahedgeHF_Identifiers_na_phrases,EurekahedgeHF_Other_na_phrases)
# #rm2(EurekahedgeHF_Instruments_Traded_na_phrases)
# 
# no_phrases_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_no_phrases"),
#                          col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange",
#                                "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
# no_phrases_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_no_phrases"),
#                          col=c("Management_Fee","Performance_Fee","Other_Fee"))
# 
# no_phrases_all <- list(no_phrases_temp1,no_phrases_temp2)
# 
# rm2(no_phrases_temp1,no_phrases_temp2)
# invisible(gc(verbose = FALSE, reset = TRUE))


no_phrases_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_na_phrases"),data_out=c("EurekahedgeHF_Stats_noreturns_no_phrases"),
                         col=NULL)
no_phrases_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_na_phrases"),data_out=c("EurekahedgeHF_Fund_Details_no_phrases"),
                         col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange",
                               "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
no_phrases_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_na_phrases"),data_out=c("EurekahedgeHF_Fee_and_Redemption_no_phrases"),
                         col=c("Management_Fee","Performance_Fee","Other_Fee"))
no_phrases_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_na_phrases"),data_out=c("EurekahedgeHF_Profile_Strategy_no_phrases"),
                         col=NULL)
no_phrases_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_na_phrases"),data_out=c("EurekahedgeHF_Identifiers_no_phrases"),
                         col=NULL)
no_phrases_temp6 <- list(data_in=c("EurekahedgeHF_Other_na_phrases"),data_out=c("EurekahedgeHF_Other_no_phrases"),
                         col=NULL)

no_phrases_all0 <- list(no_phrases_temp1,no_phrases_temp2,no_phrases_temp3,
                        no_phrases_temp4,no_phrases_temp5,no_phrases_temp6)

rm2(no_phrases_temp1,no_phrases_temp2,no_phrases_temp3,no_phrases_temp4,no_phrases_temp5,no_phrases_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

no_phrases_all1 <- llply(.data=no_phrases_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(no_phrases_all0)

no_phrases_all1 <- no_phrases_all1[order(-sapply(no_phrases_all1,"[[","size"),
                                         -sapply(no_phrases_all1,"[[","ncol"),
                                         -sapply(no_phrases_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=no_phrases_all1, .fun = function(x,phrases){
  
  # x <- no_phrases_all1[[1]]
  # phrases <- NO_Phrases
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], 
           no_to_no(data.table(get(x[[which(names(x)==c("data_in"))]])),x[[which(names(x)==c("col"))]],phrases), 
           envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
  #assign(x[[1]], no_to_no(data.table(get(x[[1]])),x[[2]],phrases), envir = .GlobalEnv)
  
  #invisible(gc(verbose = FALSE, reset = TRUE))
  
},phrases=NO_Phrases, .progress = "text")

rm2(no_phrases_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: CHANGE YES PHRASES TO 'YES'", "\n")
###############################################################################

YES_Phrases <- c("RARELY","OCCASIONALLY")

# EurekahedgeHF_Stats_noreturns_yes_phrases <- EurekahedgeHF_Stats_noreturns_no_phrases
# EurekahedgeHF_Fund_Details_yes_phrases <- EurekahedgeHF_Fund_Details_no_phrases
# EurekahedgeHF_Fee_and_Redemption_yes_phrases <- EurekahedgeHF_Fee_and_Redemption_no_phrases
# EurekahedgeHF_Profile_Strategy_yes_phrases <- EurekahedgeHF_Profile_Strategy_no_phrases
# EurekahedgeHF_Identifiers_yes_phrases <- EurekahedgeHF_Identifiers_no_phrases
# EurekahedgeHF_Other_yes_phrases <- EurekahedgeHF_Other_no_phrases
# #EurekahedgeHF_Instruments_Traded_yes_phrases <- EurekahedgeHF_Instruments_Traded_no_phrases
# 
# rm2(EurekahedgeHF_Stats_noreturns_no_phrases,EurekahedgeHF_Fund_Details_no_phrases,EurekahedgeHF_Fee_and_Redemption_no_phrases)
# rm2(EurekahedgeHF_Profile_Strategy_no_phrases,EurekahedgeHF_Identifiers_no_phrases,EurekahedgeHF_Other_no_phrases)
# #rm2(EurekahedgeHF_Instruments_Traded_no_phrases)
# 
# yes_phrases_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_yes_phrases"),
#                           col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange",
#                                 "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
# yes_phrases_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_yes_phrases"),
#                           col=c("Management_Fee","Performance_Fee","Other_Fee"))
# 
# yes_phrases_all <- list(yes_phrases_temp1,yes_phrases_temp2)
# 
# rm2(yes_phrases_temp1,yes_phrases_temp2)
# invisible(gc(verbose = FALSE, reset = TRUE))


yes_phrases_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_no_phrases"),data_out=c("EurekahedgeHF_Stats_noreturns_yes_phrases"),
                          col=NULL)
yes_phrases_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_no_phrases"),data_out=c("EurekahedgeHF_Fund_Details_yes_phrases"),
                          col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange",
                                "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
yes_phrases_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_no_phrases"),data_out=c("EurekahedgeHF_Fee_and_Redemption_yes_phrases"),
                          col=c("Management_Fee","Performance_Fee","Other_Fee"))
yes_phrases_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_no_phrases"),data_out=c("EurekahedgeHF_Profile_Strategy_yes_phrases"),
                          col=NULL)
yes_phrases_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_no_phrases"),data_out=c("EurekahedgeHF_Identifiers_yes_phrases"),
                          col=NULL)
yes_phrases_temp6 <- list(data_in=c("EurekahedgeHF_Other_no_phrases"),data_out=c("EurekahedgeHF_Other_yes_phrases"),
                          col=NULL)

yes_phrases_all0 <- list(yes_phrases_temp1,yes_phrases_temp2,yes_phrases_temp3,
                         yes_phrases_temp4,yes_phrases_temp5,yes_phrases_temp6)

rm2(yes_phrases_temp1,yes_phrases_temp2,yes_phrases_temp3,yes_phrases_temp4,yes_phrases_temp5,yes_phrases_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

yes_phrases_all1 <- llply(.data=yes_phrases_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(yes_phrases_all0)

yes_phrases_all1 <- yes_phrases_all1[order(-sapply(yes_phrases_all1,"[[","size"),
                                           -sapply(yes_phrases_all1,"[[","ncol"),
                                           -sapply(yes_phrases_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=yes_phrases_all1, .fun = function(x,phrases){
  
  # x <- yes_phrases_all1[[1]]
  # phrases <- NO_Phrases
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], 
           yes_to_yes(data.table(get(x[[which(names(x)==c("data_in"))]])),x[[which(names(x)==c("col"))]],phrases), 
           envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
  #assign(x[[1]], yes_to_yes(data.table(get(x[[1]])),x[[2]],phrases), envir = .GlobalEnv)
  
  #invisible(gc(verbose = FALSE, reset = TRUE))
  
},phrases=YES_Phrases, .progress = "text")

rm2(yes_phrases_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: CONVERT FEES TO NUMERIC", "\n")
###############################################################################

# EurekahedgeHF_Stats_noreturns_fees_numeric <- EurekahedgeHF_Stats_noreturns_yes_phrases
# EurekahedgeHF_Fund_Details_fees_numeric <- EurekahedgeHF_Fund_Details_yes_phrases
# EurekahedgeHF_Fee_and_Redemption_fees_numeric <- EurekahedgeHF_Fee_and_Redemption_yes_phrases
# EurekahedgeHF_Profile_Strategy_fees_numeric <- EurekahedgeHF_Profile_Strategy_yes_phrases
# EurekahedgeHF_Identifiers_fees_numeric <- EurekahedgeHF_Identifiers_yes_phrases
# EurekahedgeHF_Other_fees_numeric <- EurekahedgeHF_Other_yes_phrases
# #EurekahedgeHF_Instruments_Traded_fees_numeric <- EurekahedgeHF_Instruments_Traded_yes_phrases
# 
# rm2(EurekahedgeHF_Stats_noreturns_yes_phrases,EurekahedgeHF_Fund_Details_yes_phrases,EurekahedgeHF_Fee_and_Redemption_yes_phrases)
# rm2(EurekahedgeHF_Profile_Strategy_yes_phrases,EurekahedgeHF_Identifiers_yes_phrases,EurekahedgeHF_Other_yes_phrases)
# #rm2(EurekahedgeHF_Instruments_Traded_yes_phrases)
# 
# fee_numeric_temp1 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_fees_numeric"),
#                           col=c("Management_Fee","Performance_Fee","Other_Fee"))
# 
# fee_numeric_all <- list(fee_numeric_temp1)
# 
# rm2(fee_numeric_temp1)
# invisible(gc(verbose = FALSE, reset = TRUE))


fees_numeric_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_yes_phrases"),data_out=c("EurekahedgeHF_Stats_noreturns_fees_numeric"),
                           col=NULL)
fees_numeric_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_yes_phrases"),data_out=c("EurekahedgeHF_Fund_Details_fees_numeric"),
                           col=NULL)
fees_numeric_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_yes_phrases"),data_out=c("EurekahedgeHF_Fee_and_Redemption_fees_numeric"),
                           col=c("Management_Fee","Performance_Fee","Other_Fee"))
fees_numeric_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_yes_phrases"),data_out=c("EurekahedgeHF_Profile_Strategy_fees_numeric"),
                           col=NULL)
fees_numeric_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_yes_phrases"),data_out=c("EurekahedgeHF_Identifiers_fees_numeric"),
                           col=NULL)
fees_numeric_temp6 <- list(data_in=c("EurekahedgeHF_Other_yes_phrases"),data_out=c("EurekahedgeHF_Other_fees_numeric"),
                           col=NULL)

fees_numeric_all0 <- list(fees_numeric_temp1,fees_numeric_temp2,fees_numeric_temp3,
                          fees_numeric_temp4,fees_numeric_temp5,fees_numeric_temp6)

rm2(fees_numeric_temp1,fees_numeric_temp2,fees_numeric_temp3,fees_numeric_temp4,fees_numeric_temp5,fees_numeric_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

fees_numeric_all1 <- llply(.data=fees_numeric_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(fees_numeric_all0)

fees_numeric_all1 <- fees_numeric_all1[order(-sapply(fees_numeric_all1,"[[","size"),
                                             -sapply(fees_numeric_all1,"[[","ncol"),
                                             -sapply(fees_numeric_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=fees_numeric_all1, .fun = function(x,unknowns){
  
  # x <- fees_numeric_all1[[1]]
  # unknowns <- unknowns_strings
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    data_temp <- data.table(get(x[[which(names(x)==c("data_in"))]]))
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    col_temp <- x[[which(names(x)==c("col"))]]
    
    for (i in 1:length(col_temp)) {
      
      # i <- 1
      # i <- 2
      # i <- 3
      
      data_temp_convert <-  data.table(matrix(NA, ncol=7, nrow=nrow(data_temp), 
                                              dimnames=list(c(), c("org","no_char","no_punct","first_num","remove_blanks","flag","convert"))))
      
      cols <- c("org")
      for (k in cols) {
        
        set(data_temp_convert, i=NULL, j=k, value=as.character(data_temp_convert[[k]]))
        set(data_temp_convert, i=NULL, j=k, value=data_temp[[which(colnames(data_temp) %in% col_temp[i])]])
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" {2,}", replacement=" ", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="^\\s+|\\s+$", replacement="", data_temp_convert[[k]], perl=TRUE))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      #data_temp_convert <- unique(data_temp_convert)
      
      #Remove Characters
      data_temp_convert[, no_char := org, by = NULL]
      
      cols <- c("no_char")
      for (k in cols) {
        
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="([[:alpha:]])", replacement="", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" {2,}", replacement=" ", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="^\\s+|\\s+$", replacement="", data_temp_convert[[k]], perl=TRUE))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      
      #Remove Punctuation (except % and $)
      data_temp_convert[, no_punct := no_char, by = NULL]
      
      cols <- c("no_punct")
      for (k in cols) {
        
        #set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="([[:punct:]])", replacement="", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="(:|,|;|=|-|&|\\+|_|\\?|!|/)", replacement="", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="\\. ", replacement="\\.", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" \\. ", replacement="\\.", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="\\.{2,}", replacement="\\.", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" {2,}", replacement=" ", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="^\\s+|\\s+$", replacement="", data_temp_convert[[k]], perl=TRUE))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      
      #Get first number   
      data_temp_convert[, first_num := no_punct, by = NULL]
      
      cols <- c("first_num")
      for (k in cols) {
        
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" .*$", replacement="", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" {2,}", replacement=" ", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="^\\s+|\\s+$", replacement="", data_temp_convert[[k]], perl=TRUE))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      
      #Remove Blanks
      data_temp_convert[, remove_blanks := first_num, by = NULL]
      
      cols <- c("remove_blanks")
      for (k in cols) {
        
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="(\\$|%)", replacement="", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern=" {2,}", replacement=" ", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=gsub(pattern="^\\s+|\\s+$", replacement="", data_temp_convert[[k]], perl=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=unknownToNA(data_temp_convert[[k]], unknown=unknowns,force=TRUE))
        set(data_temp_convert, i=NULL, j=k, value=ifelse(is.na(data_temp_convert[[k]]),NA,data_temp_convert[[k]]))
        set(data_temp_convert, i=NULL, j=k, value=as.numeric(data_temp_convert[[k]]))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      #Flag number type
      # decimal: flag=1
      # percent: flag=2
      # bp:      flag=3
      # fixed:   flag=4
      
      cols <- c("flag")
      for (k in cols) {
        
        set(data_temp_convert, i=NULL, j=k, value=ifelse((data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]>=0 & data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]<1.0),1,
                                                         ifelse((data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]>=1.0 & data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]<100.0),2,
                                                                ifelse((data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]>=100.0 & data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]<10000),3,
                                                                       ifelse((data_temp_convert[[which(colnames(data_temp_convert) %in% "remove_blanks")]]>=10000),4,NA)))))
        set(data_temp_convert, i=NULL, j=k, value=ifelse(grepl("%",data_temp_convert[[which(colnames(data_temp_convert) %in% "first_num")]]),2,data_temp_convert[[k]]))
        set(data_temp_convert, i=NULL, j=k, value=ifelse(grepl("\\$",data_temp_convert[[which(colnames(data_temp_convert) %in% "first_num")]]),4,data_temp_convert[[k]]))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      
      #Convert to decimals  
      data_temp_convert[, convert := remove_blanks, by = NULL]
      
      cols <- c("convert")
      for (k in cols) {
        
        set(data_temp_convert, i=NULL, j=k, value=ifelse(data_temp_convert[[k]]==2,data_temp_convert[[which(colnames(data_temp_convert) %in% "convert")]]/100,data_temp_convert[[which(colnames(data_temp_convert) %in% "convert")]]))
        set(data_temp_convert, i=NULL, j=k, value=ifelse(data_temp_convert[[k]]==3,data_temp_convert[[which(colnames(data_temp_convert) %in% "convert")]]/10000,data_temp_convert[[which(colnames(data_temp_convert) %in% "convert")]]))
        set(data_temp_convert, i=NULL, j=k, value=ifelse(data_temp_convert[[k]]==4,NA,data_temp_convert[[which(colnames(data_temp_convert) %in% "convert")]]))
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      cols <- col_temp[i]
      for (k in cols) {
        
        set(data_temp, i=NULL, j=k, value=data_temp_convert[[which(colnames(data_temp_convert) %in% "convert")]])
        
        #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
      }
      rm(k,cols)
      
      rm(data_temp_convert)
      
      invisible(gc(verbose = FALSE, reset = TRUE))
    }
    
    assign(x[[which(names(x)==c("data_out"))]], data_temp, envir = .GlobalEnv)
    
    rm(data_temp)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
},unknowns=unknowns_strings, .progress = "text")

rm2(fees_numeric_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: CHANGE Y/N TO BINARY", "\n")
###############################################################################

# EurekahedgeHF_Stats_noreturns_yn_binary <- EurekahedgeHF_Stats_noreturns_fees_numeric
# EurekahedgeHF_Fund_Details_yn_binary <- EurekahedgeHF_Fund_Details_fees_numeric
# EurekahedgeHF_Fee_and_Redemption_yn_binary <- EurekahedgeHF_Fee_and_Redemption_fees_numeric
# EurekahedgeHF_Profile_Strategy_yn_binary <- EurekahedgeHF_Profile_Strategy_fees_numeric
# EurekahedgeHF_Identifiers_yn_binary <- EurekahedgeHF_Identifiers_fees_numeric
# EurekahedgeHF_Other_yn_binary <- EurekahedgeHF_Other_fees_numeric
# #EurekahedgeHF_Instruments_Traded_yn_binary <- EurekahedgeHF_Instruments_Traded_fees_numeric
# 
# rm2(EurekahedgeHF_Stats_noreturns_fees_numeric,EurekahedgeHF_Fund_Details_fees_numeric,EurekahedgeHF_Fee_and_Redemption_fees_numeric)
# rm2(EurekahedgeHF_Profile_Strategy_fees_numeric,EurekahedgeHF_Identifiers_fees_numeric,EurekahedgeHF_Other_fees_numeric)
# #rm2(EurekahedgeHF_Instruments_Traded_fees_numeric)
# 
# yn_binary_temp1 <- list(data=c("EurekahedgeHF_Stats_noreturns_yn_binary"),
#                         col=c("Flagship","Closed","Limited","Dead"))
# yn_binary_temp2 <- list(data=c("EurekahedgeHF_Fund_Details_yn_binary"),
#                         col=c("Flagship","Closed","Limited","Dead",
#                               "Invest_In_Private_Placements","Managed_Accounts_Offered","UCITS_combcol",
#                               "Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange"))
# yn_binary_temp3 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_yn_binary"),
#                         col=c("Flagship","Closed","Limited","Dead"))
# yn_binary_temp4 <- list(data=c("EurekahedgeHF_Profile_Strategy_yn_binary"),
#                         col=c("Flagship","Closed","Limited","Dead"))
# yn_binary_temp5 <- list(data=c("EurekahedgeHF_Identifiers_yn_binary"),
#                         col=c("Flagship","Closed","Limited","Dead"))
# yn_binary_temp6 <- list(data=c("EurekahedgeHF_Other_yn_binary"),
#                         col=c("Flagship","Closed","Limited","Dead"))
# 
# yn_binary_all <- list(yn_binary_temp1,yn_binary_temp2,yn_binary_temp3,
#                       yn_binary_temp4,yn_binary_temp5,yn_binary_temp6)
# 
# rm2(yn_binary_temp1,yn_binary_temp2,yn_binary_temp3)
# rm2(yn_binary_temp4,yn_binary_temp5,yn_binary_temp6)

yn_binary_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_fees_numeric"),data_out=c("EurekahedgeHF_Stats_noreturns_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_fees_numeric"),data_out=c("EurekahedgeHF_Fund_Details_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead",
                              "Invest_In_Private_Placements","Managed_Accounts_Offered","UCITS_combcol",
                              "Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange"))
yn_binary_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_fees_numeric"),data_out=c("EurekahedgeHF_Fee_and_Redemption_yn_binary"),
                        col=c("Management_Fee","Performance_Fee","Other_Fee"))
yn_binary_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_fees_numeric"),data_out=c("EurekahedgeHF_Profile_Strategy_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_fees_numeric"),data_out=c("EurekahedgeHF_Identifiers_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp6 <- list(data_in=c("EurekahedgeHF_Other_fees_numeric"),data_out=c("EurekahedgeHF_Other_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))

yn_binary_all0 <- list(yn_binary_temp1,yn_binary_temp2,yn_binary_temp3,
                       yn_binary_temp4,yn_binary_temp5,yn_binary_temp6)

rm2(yn_binary_temp1,yn_binary_temp2,yn_binary_temp3,yn_binary_temp4,yn_binary_temp5,yn_binary_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

yn_binary_all1 <- llply(.data=yn_binary_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(yn_binary_all0)

yn_binary_all1 <- yn_binary_all1[order(-sapply(yn_binary_all1,"[[","size"),
                                       -sapply(yn_binary_all1,"[[","ncol"),
                                       -sapply(yn_binary_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))

### Fix Dividend Policy

yn_binary_dividend_policy <-  data.frame(matrix(NA, ncol=3, nrow=nrow(EurekahedgeHF_Fund_Details_fees_numeric), 
                                                dimnames=list(c(), c("org","replace","final"))),stringsAsFactors=FALSE)
yn_binary_dividend_policy[,"org"] <- EurekahedgeHF_Fund_Details_fees_numeric[,"Dividend_Policy"]

#yn_binary_dividend_policy <- unique(yn_binary_dividend_policy)

yn_binary_dividend_policy[,"replace"] <- yn_binary_dividend_policy[,"org"]
yn_binary_dividend_policy[,"replace"] <- ifelse(is.na(yn_binary_dividend_policy[,"replace"]),
                                                NA,ifelse(grepl("(Yes|No)",yn_binary_dividend_policy[,"replace"]),
                                                          yn_binary_dividend_policy[,"replace"],"REPLACED"))

yn_binary_dividend_policy[,"final"] <- yn_binary_dividend_policy[,"replace"]
# yn_binary_dividend_policy[,"final"] <- ifelse(is.na(yn_binary_dividend_policy[,"final"]),
#                                               NA,ifelse(grepl("(REPLACED)",yn_binary_dividend_policy[,"final"]),
#                                                         "Yes",yn_binary_dividend_policy[,"final"]))
yn_binary_dividend_policy[,"final"] <- ifelse(is.na(yn_binary_dividend_policy[,"final"]),
                                              NA,ifelse(grepl("(REPLACED)",yn_binary_dividend_policy[,"final"]),
                                                        NA,yn_binary_dividend_policy[,"final"]))

#EurekahedgeHF_Fund_Details_fees_numeric[,"Dividend_Policy"] <- yn_binary_dividend_policy[,"final"]

cols <- c("Dividend_Policy")
for (k in cols) {
  
  set(EurekahedgeHF_Fund_Details_fees_numeric, i=NULL, j=k, value=yn_binary_dividend_policy[[which(colnames(yn_binary_dividend_policy) %in% "final")]])
  
  #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
}
rm(k,cols)

rm2(yn_binary_dividend_policy)
invisible(gc(verbose = FALSE, reset = TRUE))


### Fix Fund Closed

yn_binary_fund_closed <-  data.frame(matrix(NA, ncol=3, nrow=nrow(EurekahedgeHF_Fund_Details_fees_numeric), 
                                            dimnames=list(c(), c("org","replace","final"))),stringsAsFactors=FALSE)
yn_binary_fund_closed[,"org"] <- EurekahedgeHF_Fund_Details_fees_numeric[,"Fund_Closed"]

#yn_binary_fund_closed <- unique(yn_binary_fund_closed)

yn_binary_fund_closed[,"replace"] <- yn_binary_fund_closed[,"org"]
yn_binary_fund_closed[,"replace"] <- ifelse(is.na(yn_binary_fund_closed[,"replace"]),
                                            NA,ifelse(grepl("(Closed|closed|Hard Closed|Hard closed|Hard-Closed|Hard-closed|Soft Closed|Soft closed|Soft-Closed|Soft-closed)",yn_binary_fund_closed[,"replace"]),
                                                      "REPLACED",yn_binary_fund_closed[,"replace"]))

yn_binary_fund_closed[,"final"] <- yn_binary_fund_closed[,"replace"]
yn_binary_fund_closed[,"final"] <- ifelse(is.na(yn_binary_fund_closed[,"final"]),
                                          NA,ifelse(grepl("(REPLACED)",yn_binary_fund_closed[,"final"]),
                                                    "Yes",yn_binary_fund_closed[,"final"]))
# yn_binary_fund_closed[,"final"] <- ifelse(is.na(yn_binary_fund_closed[,"final"]),
#                                           NA,ifelse(grepl("(REPLACED)",yn_binary_fund_closed[,"final"]),
#                                                     NA,yn_binary_fund_closed[,"final"]))

#EurekahedgeHF_Fund_Details_fees_numeric[,"Fund_Closed"] <- yn_binary_fund_closed[,"final"]

cols <- c("Fund_Closed")
for (k in cols) {
  
  set(EurekahedgeHF_Fund_Details_fees_numeric, i=NULL, j=k, value=yn_binary_fund_closed[[which(colnames(yn_binary_fund_closed) %in% "final")]])
  
  #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
}
rm(k,cols)

rm2(yn_binary_fund_closed)
invisible(gc(verbose = FALSE, reset = TRUE))


### Fix High Water Mark

yn_binary_high_water_mark <-  data.frame(matrix(NA, ncol=3, nrow=nrow(EurekahedgeHF_Fund_Details_fees_numeric), 
                                                dimnames=list(c(), c("org","replace","final"))),stringsAsFactors=FALSE)
yn_binary_high_water_mark[,"org"] <- EurekahedgeHF_Fund_Details_fees_numeric[,"High_Water_Mark"]

#yn_binary_high_water_mark <- unique(yn_binary_high_water_mark)

yn_binary_high_water_mark[,"replace"] <- yn_binary_high_water_mark[,"org"]
yn_binary_high_water_mark[,"replace"] <- ifelse(is.na(yn_binary_high_water_mark[,"replace"]),
                                                NA,ifelse(grepl("(Yes|No)",yn_binary_high_water_mark[,"replace"]),
                                                          yn_binary_high_water_mark[,"replace"],"REPLACED"))

yn_binary_high_water_mark[,"final"] <- yn_binary_high_water_mark[,"replace"]
# yn_binary_high_water_mark[,"final"] <- ifelse(is.na(yn_binary_high_water_mark[,"final"]),
#                                               NA,ifelse(grepl("(REPLACED)",yn_binary_high_water_mark[,"final"]),
#                                                         "Yes",yn_binary_high_water_mark[,"final"]))
yn_binary_high_water_mark[,"final"] <- ifelse(is.na(yn_binary_high_water_mark[,"final"]),
                                              NA,ifelse(grepl("(REPLACED)",yn_binary_high_water_mark[,"final"]),
                                                        NA,yn_binary_high_water_mark[,"final"]))

#EurekahedgeHF_Fund_Details_fees_numeric[,"High_Water_Mark"] <- yn_binary_high_water_mark[,"final"]

cols <- c("High_Water_Mark")
for (k in cols) {
  
  set(EurekahedgeHF_Fund_Details_fees_numeric, i=NULL, j=k, value=yn_binary_high_water_mark[[which(colnames(yn_binary_high_water_mark) %in% "final")]])
  
  #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
}
rm(k,cols)

rm2(yn_binary_high_water_mark)
invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=yn_binary_all1, .fun = function(x,unknowns){
  
  # x <- yn_binary_all1[[1]]
  # x <- yn_binary_all1[[2]]
  # unknowns <- unknowns_strings
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    col_temp <- x[[which(names(x)==c("col"))]]
    
    col_order <- data.frame(col=colnames(get(x[[which(names(x)==c("data_in"))]])),order=NA,stringsAsFactors=FALSE)
    col_order[,"order"] <- seq(1,nrow(col_order))
    
    col_order_bin <- data.frame(col=colnames(get(x[[which(names(x)==c("data_in"))]])),order=NA,stringsAsFactors=FALSE)
    col_order_bin[,"order"] <- seq(1,nrow(col_order_bin))
    col_order_bin[,"col"] <- ifelse(col_order_bin[,"col"] %in% col_temp,paste(col_order_bin[,"col"],"_bin",sep=""),col_order_bin[,"col"])
    
    col_order_all <- unique(rbind(col_order_bin,col_order))
    
    col_order_all <- col_order_all[order(col_order_all[,"order"]),]
    row.names(col_order_all) <- seq(nrow(col_order_all))
    
    bin_cols <- paste(col_temp,"_bin",sep="")
    
    data_temp2 <-  data.table(get(x[[which(names(x)==c("data_in"))]]), matrix(NA, ncol=length(bin_cols), nrow=1, dimnames=list(c(), bin_cols)))
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    cols <- col_temp
    for (k in cols) {
      
      set(data_temp2, i=NULL, j=paste(k,"_bin",sep = ""), value=data_temp2[[k]])
      
      #cat("Loop: ",which(cols==k)," of ",length(cols), "\n")
    }
    rm(k,cols)
    
    setcolorder(data_temp2, col_order_all[,"col"])
    
    data_temp2 <- yn_to_binary(data_temp2,bin_cols)
    
    cols <- bin_cols
    for (k in cols) 
    {
      #k <- 1
      set(data_temp2, i=NULL, j=k, value=unknownToNA(data_temp2[[k]], unknown=unknowns,force=TRUE))
      set(data_temp2, i=NULL, j=k, value=ifelse(is.na(data_temp2[[k]]),NA,data_temp2[[k]]))
    }
    rm(k,cols)
    
    data_temp2[,droprow := rowSums(is.na(.SD))<ncol(data_temp2),.SDcols=colnames(data_temp2)]
    #set(data_temp2,i=NULL,j=which(colnames(data_temp2)==c("droprow")),value=NULL)
    data_temp2 <- data_temp2[(droprow)][,droprow:=NULL][]
    
    assign(x[[which(names(x)==c("data_out"))]], data_temp2, envir = .GlobalEnv)
    
    rm(data_temp2)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
},unknowns=unknowns_strings, .progress = "text")

rm2(yn_binary_all1)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: OUTPUT", "\n")
###############################################################################

#Check to see if common_cols folder exists.  If not, create it.
final_expand_folder_path <- paste(output_directory, "Final_Expand", sep = "//", collapse = "//")  
create_directory(final_expand_folder_path,remove=1)

write.csv(EurekahedgeHF_Stats_noreturns_yn_binary, file=paste(final_expand_folder_path,"//","EurekahedgeHF_Stats_noreturns_part1",".csv",sep=""),row.names=FALSE)
write.csv(EurekahedgeHF_Fund_Details_yn_binary, file=paste(final_expand_folder_path,"//","EurekahedgeHF_Fund_Details_part1",".csv",sep=""),row.names=FALSE)
write.csv(EurekahedgeHF_Fee_and_Redemption_yn_binary, file=paste(final_expand_folder_path,"//","EurekahedgeHF_Fee_and_Redemption_part1",".csv",sep=""),row.names=FALSE)
write.csv(EurekahedgeHF_Profile_Strategy_yn_binary, file=paste(final_expand_folder_path,"//","EurekahedgeHF_Profile_Strategy_part1",".csv",sep=""),row.names=FALSE)
write.csv(EurekahedgeHF_Identifiers_yn_binary, file=paste(final_expand_folder_path,"//","EurekahedgeHF_Identifiers_part1",".csv",sep=""),row.names=FALSE)
write.csv(EurekahedgeHF_Other_yn_binary, file=paste(final_expand_folder_path,"//","EurekahedgeHF_Other_part1",".csv",sep=""),row.names=FALSE)

rm2(final_expand_folder_path)
rm2(EurekahedgeHF_Stats_noreturns_yn_binary,EurekahedgeHF_Fund_Details_yn_binary)
rm2(EurekahedgeHF_Fee_and_Redemption_yn_binary,EurekahedgeHF_Profile_Strategy_yn_binary)
rm2(EurekahedgeHF_Identifiers_yn_binary,EurekahedgeHF_Other_yn_binary)

