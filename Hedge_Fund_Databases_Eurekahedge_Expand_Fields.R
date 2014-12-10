# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Expand_Fields.R
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

external_packages <- c("compare","cwhmisc","data.table","DataCombine","fastmatch","foreign","formatR","gdata",
                       "gtools","Hmisc","installr","knitr","koRpus","lmtest","lubridate","markdown","memisc","mitools",
                       "pander","pbapply","plm","plyr","psych","quantreg","R.oo","R2wd","reporttools","reshape2","rms","RSQLite",
                       "sandwich","sqldf","stargazer","stringr","texreg","tm","UsingR","xtable","zoo")
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
final_folder_expand_path <- paste(output_directory, "Final_Expand", sep = "//", collapse = "//")  
create_directory(final_folder_expand_path,remove=1)

final_folder_files0 <- data.frame(files=list.files(path=final_folder_path),file_name=NA,import=NA,stringsAsFactors=FALSE)

#final_folder_files <- final_folder_files0[!grepl(".TXT|.txt", final_folder_files0[,"files"]),]
final_folder_files <- final_folder_files0[grepl(".CSV|.csv", final_folder_files0[,"files"]),]

rm2(final_folder_files0)

final_folder_files[,"file_name"] <- final_folder_files[,"files"]
final_folder_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=final_folder_files[,"file_name"])

final_folder_files[,"import"] <- ifelse(grepl("(Stats|Fund_Detail|Fee_and_Redemption|Profile_Strategy|Identifier|Instruments_Traded)",final_folder_files[,"file_name"]),1,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(Other)",final_folder_files[,"file_name"]),2,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(NAV_AUM_Ret)",final_folder_files[,"file_name"]),3,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(is.na(final_folder_files[,"import"]),0,final_folder_files[,"import"])


###############################################################################
cat("SECTION: IMPORT FILES", "\n")
###############################################################################

a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(1),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[2,]
  # x <- final_folder_files[3,]
  # x <- final_folder_files[4,] 
  
  # directory_in <- final_folder_path
  # unknowns <- unknowns_strings
  
  #input <- read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)
  
  input <- data.frame(pull_trim=NA,
                      read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      stringsAsFactors=FALSE)
  
  colnames(input)[match("pull",names(input))] <- "pull_trim2"
  
  input[,"pull_trim"] <- input[,"pull_trim2"] 
  input[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", input[,"pull_trim"])
  
  input[,"pull_trim"] <- as.character(input[,"pull_trim"])
  
  input[,"pull_trim2"] <- gsub(pattern="_Fund_Details_Fund_Details", replacement="", x=input[,"pull_trim2"])
  input[,"pull_trim2"] <- gsub(pattern="_Fund_Details_Fee_and_Redemption_Structure", replacement="", x=input[,"pull_trim2"])
  input[,"pull_trim2"] <- gsub(pattern="_Fund_Details_Unique_Identifiers", replacement="", x=input[,"pull_trim2"])
  input[,"pull_trim2"] <- gsub(pattern="_Fund_Details_Profile_Strategy_Description", replacement="", x=input[,"pull_trim2"])
  input[,"pull_trim2"] <- gsub(pattern="_Fund_Details_Statistics", replacement="", x=input[,"pull_trim2"])
  input[,"pull_trim2"] <- gsub(pattern="_Instruments_Traded", replacement="", x=input[,"pull_trim2"])
  
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
  
  #colnames(input) <- tolower(colnames(input))
  
  #input  <- input[order(input[,"Fund_ID"],input[,"Fund_Name"],input[,"pull"]),]
  input  <- input[order(input[,"Fund_ID"],input[,"pull_trim"],input[,"pull_trim2"]),]
  row.names(input) <- seq(nrow(input))
  
  if("Date_Added" %in% colnames(input)){input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b %Y")} 
  if("Dead_Date" %in% colnames(input)){input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b %Y")} 
  
  assign(x[,"file_name"], input, envir = .GlobalEnv)
  
  rm(input)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
}, directory_in=final_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))

a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(2),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[8,]
  
  # directory_in <- final_folder_path
  # unknowns <- unknowns_strings
  
  #input <- read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)
  
  input <- data.frame(read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      stringsAsFactors=FALSE)
  
  input[,"pull_trim"] <- as.character(input[,"pull_trim"])
  
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
  
  #colnames(input) <- tolower(colnames(input))
  
  #input  <- input[order(input[,"Fund_ID"],input[,"Fund_Name"],input[,"pull"]),]
  input  <- input[order(input[,"Fund_ID"],input[,"pull_trim"],input[,"pull_trim2"]),]
  row.names(input) <- seq(nrow(input))
  
  if("Date_Added" %in% colnames(input)){input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b %Y")} 
  if("Dead_Date" %in% colnames(input)){input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b %Y")} 
  
  assign(x[,"file_name"], input, envir = .GlobalEnv)
  
  rm(input)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
}, directory_in=final_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))


a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(3),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[8,]
  
  # directory_in <- final_folder_path
  # unknowns <- unknowns_strings
  
  #input <- read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)
  
  input_cols_keep <- c("pull_trim","pull","Fund_ID","Dead_Date","yr","month","date","bad_min","bad_max","AUM")                        
  input <- data.frame(read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)[,input_cols_keep],
                      stringsAsFactors=FALSE)
  
  #input <- data.frame(read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
  #                    stringsAsFactors=FALSE)
  
  colnames(input)[match("AUM",names(input))] <- "pull_trim2"
  
  input[,"pull_trim2"] <- input[,"pull"]
  input[,"pull_trim2"] <- gsub(pattern="_NAV_AUM", replacement="", x=input[,"pull_trim2"])
  
  input <- input[,c("pull_trim","pull_trim2","pull",
                    colnames(input)[!(colnames(input) %in% c("pull_trim","pull_trim2","pull"))])]
  
  input[,"pull_trim"] <- as.character(input[,"pull_trim"])
  
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
  
  #colnames(input) <- tolower(colnames(input))
  
  #input  <- input[order(input[,"Fund_ID"],input[,"Fund_Name"],input[,"pull"]),]
  input  <- input[order(input[,"Fund_ID"],input[,"pull_trim"],input[,"pull_trim2"]),]
  row.names(input) <- seq(nrow(input))
  
  if("Date_Added" %in% colnames(input)){input[,"Date_Added"] <- as.yearmon(input[,"Date_Added"],format="%b %Y")} 
  if("Dead_Date" %in% colnames(input)){input[,"Dead_Date"] <- as.yearmon(input[,"Dead_Date"],format="%b %Y")} 
  
  #assign(x[,"file_name"], input, envir = .GlobalEnv)
  assign("Merge_IDs", input, envir = .GlobalEnv)
  
  rm(input)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
}, directory_in=final_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: ADD PULL IDS AND DEAD DATES", "\n")
###############################################################################

EurekahedgeHF_Stats_noreturns_deaddates <- EurekahedgeHF_Stats_noreturns
EurekahedgeHF_Fund_Details_deaddates <- EurekahedgeHF_Fund_Details
EurekahedgeHF_Fee_and_Redemption_deaddates <- EurekahedgeHF_Fee_and_Redemption
EurekahedgeHF_Profile_Strategy_deaddates <- EurekahedgeHF_Profile_Strategy
EurekahedgeHF_Identifiers_deaddates <- EurekahedgeHF_Identifiers
EurekahedgeHF_Instruments_Traded_deaddates <- EurekahedgeHF_Instruments_Traded
EurekahedgeHF_Other_deaddates <- EurekahedgeHF_Other

rm2(EurekahedgeHF_Stats_noreturns,EurekahedgeHF_Fund_Details,EurekahedgeHF_Fee_and_Redemption)
rm2(EurekahedgeHF_Profile_Strategy,EurekahedgeHF_Identifiers,EurekahedgeHF_Instruments_Traded,EurekahedgeHF_Other)

fix_dead_dates_temp1 <- list(data=c("EurekahedgeHF_Stats_noreturns_deaddates"))
fix_dead_dates_temp2 <- list(data=c("EurekahedgeHF_Fund_Details_deaddates"))
fix_dead_dates_temp3 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_deaddates"))
fix_dead_dates_temp4 <- list(data=c("EurekahedgeHF_Profile_Strategy_deaddates"))
fix_dead_dates_temp5 <- list(data=c("EurekahedgeHF_Identifiers_deaddates"))
fix_dead_dates_temp6 <- list(data=c("EurekahedgeHF_Instruments_Traded_deaddates"))

fix_dead_dates_all <- list(fix_dead_dates_temp1,fix_dead_dates_temp2,fix_dead_dates_temp3,
                           fix_dead_dates_temp4,fix_dead_dates_temp5,fix_dead_dates_temp6)

rm2(fix_dead_dates_temp1,fix_dead_dates_temp2,fix_dead_dates_temp3)
rm2(fix_dead_dates_temp4,fix_dead_dates_temp5,fix_dead_dates_temp6)

l_ply(.data=fix_dead_dates_all, .fun = function(x,ids,merge_data){
  
  # x <- fix_dead_dates_all[[1]]
  # ids <- c("pull_trim","pull_trim2")
  # merge_data <- "Merge_IDs"
  
  require(data.table)
  
  #data_name_temp <- x[[1]]
  #data_temp <- get(x[[1]])
  
  cat(x[[1]], "\n")
  
  files_temp1 <- x[[1]]
  files_temp1_trim <- data.frame(file=x[[1]],row_str=NA,col_str=NA,file_str=NA,stringsAsFactors=FALSE)
  #files_temp1_trim[,"row_str"] <- paste(files_temp1_trim[,"file"],"[",",","'pull_trim2'","]","==","'",x[,"pull_trim2"],"'",sep="")
  files_temp1_trim[,"row_str"] <- ""
  files_temp1_trim[,"col_str"] <- paste("!(colnames(",files_temp1_trim[,"file"],") %in% c('Dead_Date'))",sep="")
  files_temp1_trim[,"file_str"] <- paste(files_temp1_trim[,"file"],"[",files_temp1_trim[,"row_str"],",",files_temp1_trim[,"col_str"],"]",sep="")
  
  rm(files_temp1)
  
  files_temp2 <- merge_data
  files_temp2_trim <- data.frame(file=files_temp2[!is.na(files_temp2)],row_str=NA,col_str=NA,file_str=NA,stringsAsFactors=FALSE)
  #files_temp2_trim[,"row_str"] <- paste(files_temp2_trim[,"file"],"[",",","'pull_trim2'","]","==","'",x[,"pull_trim2"],"'",sep="")
  files_temp2_trim[,"row_str"] <- ""
  files_temp2_trim[,"col_str"] <- paste("!(colnames(",files_temp2_trim[,"file"],") %in% c('pull'))",sep="")
  files_temp2_trim[,"file_str"] <- paste(files_temp2_trim[,"file"],"[",files_temp2_trim[,"row_str"],",",files_temp2_trim[,"col_str"],"]",sep="")
  
  rm(files_temp2)
  
  files_temp_all <- rbind(files_temp2_trim,files_temp1_trim)
  
  rm(files_temp2_trim,files_temp1_trim)
  
  merge_temp <- eval(parse(text=files_temp_all[1,"file_str"]))

  #for (i in 1:nrow(files_temp_all)) {print(paste(files_temp_all[i,"file"],": nrow = ",nrow(eval(parse(text=files_temp_all[i,"file_str"]))),sep=""))}
  
  for (i in 2:nrow(files_temp_all)) {
    
    # i <- 2
    
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
  
  assign(x[[1]], merge_temp, envir = .GlobalEnv)
  
  rm(merge_temp,order_ids)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
},ids=c("pull_trim","pull_trim2"), merge_data="Merge_IDs", .progress = "text")

rm2(Merge_IDs,fix_dead_dates_all)


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

EurekahedgeHF_Stats_noreturns_comments <- EurekahedgeHF_Stats_noreturns_deaddates
EurekahedgeHF_Fund_Details_comments <- EurekahedgeHF_Fund_Details_deaddates
EurekahedgeHF_Fee_and_Redemption_comments <- EurekahedgeHF_Fee_and_Redemption_deaddates
EurekahedgeHF_Profile_Strategy_comments <- EurekahedgeHF_Profile_Strategy_deaddates
EurekahedgeHF_Identifiers_comments <- EurekahedgeHF_Identifiers_deaddates
EurekahedgeHF_Other_comments <- EurekahedgeHF_Other_deaddates

rm2(EurekahedgeHF_Stats_noreturns_deaddates,EurekahedgeHF_Fund_Details_deaddates,EurekahedgeHF_Fee_and_Redemption_deaddates)
rm2(EurekahedgeHF_Profile_Strategy_deaddates,EurekahedgeHF_Identifiers_deaddates,EurekahedgeHF_Other_deaddates)

strip_comments_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_comments"),
                             col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange"))
strip_comments_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_comments"),
                             col=c("Management_Fee","Performance_Fee","Other_Fee"))

strip_comments_all <- list(strip_comments_temp1,strip_comments_temp2)

rm2(strip_comments_temp1,strip_comments_temp2)

l_ply(.data=strip_comments_all, .fun = function(x){
  
  # x <- strip_comments_all[[1]]
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  col_order <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order[,"order"] <- seq(1,nrow(col_order))
  
  col_order_org <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order_org[,"order"] <- seq(1,nrow(col_order_org))
  col_order_org[,"col"] <- ifelse(col_order_org[,"col"] %in% col_temp,paste(col_order_org[,"col"],"_org",sep=""),col_order_org[,"col"])
  
  col_order_comments <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order_comments[,"order"] <- seq(1,nrow(col_order_comments))
  col_order_comments[,"col"] <- ifelse(col_order_comments[,"col"] %in% col_temp,paste(col_order_comments[,"col"],"_comments",sep=""),col_order_comments[,"col"])
  
  col_order_all <- unique(rbind(rbind(col_order_org,col_order),col_order_comments))
  
  col_order_all <- col_order_all[order(col_order_all[,"order"]),]
  row.names(col_order_all) <- seq(nrow(col_order_all))
  
  #Rename original columns
  data_temp <- rename.vars(data_temp, col_temp, paste(col_temp,"_org",sep=""),info=FALSE)
  
  strip_cols <- c(col_temp, paste(col_temp,"_comments",sep=""))
  
  data_temp2 <-  data.frame(data_temp, matrix(NA, ncol=length(strip_cols), nrow=nrow(data_temp), dimnames=list(c(), strip_cols)), stringsAsFactors=FALSE)
  
  #Get column order
  
  #data_temp2 <- data_temp2[,sort(colnames(data_temp2), decreasing = FALSE)]
  data_temp2 <- data_temp2[,col_order_all[,"col"]]
  
  data_temp2 <- strip_comments(data_temp2,col_temp)
  data_temp2 <- as.data.frame(data_temp2,stringsAsFactors=FALSE)
  
  #Get text before comments
  data_temp2 <- create_noncomments(data_temp2,col_temp)
  data_temp2 <- as.data.frame(data_temp2,stringsAsFactors=FALSE)
  
  assign(x[[1]], data_temp2, envir = .GlobalEnv)
  
  gc()
  
}, .progress = "text")

rm2(strip_comments_all)


###############################################################################
cat("SECTION: CHECK FOR UNKNOWNS", "\n")
###############################################################################

EurekahedgeHF_Stats_noreturns_unknowns <- EurekahedgeHF_Stats_noreturns_comments
EurekahedgeHF_Fund_Details_unknowns <- EurekahedgeHF_Fund_Details_comments
EurekahedgeHF_Fee_and_Redemption_unknowns <- EurekahedgeHF_Fee_and_Redemption_comments
EurekahedgeHF_Profile_Strategy_unknowns <- EurekahedgeHF_Profile_Strategy_comments
EurekahedgeHF_Identifiers_unknowns <- EurekahedgeHF_Identifiers_comments
EurekahedgeHF_Other_unknowns <- EurekahedgeHF_Other_comments

rm2(EurekahedgeHF_Stats_noreturns_comments,EurekahedgeHF_Fund_Details_comments,EurekahedgeHF_Fee_and_Redemption_comments)
rm2(EurekahedgeHF_Profile_Strategy_comments,EurekahedgeHF_Identifiers_comments,EurekahedgeHF_Other_comments)

unknowns_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_unknowns"),
                       col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange","Exchange_Name"))
unknowns_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_unknowns"),
                       col=c("Management_Fee","Performance_Fee","Other_Fee"))

unknowns_all <- list(unknowns_temp1,unknowns_temp2)

rm2(unknowns_temp1,unknowns_temp2)

l_ply(.data=unknowns_all, .fun = function(x,unknowns){
  
  # x <- unknowns_all[[1]]
  # unknowns <- unknowns_strings
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  data_temp <- data.table(data_temp)[, (col_temp) := llply(.SD, vector_clean_na,unknowns=unknowns,.progress = "text"), .SDcols = col_temp]
  data_temp <- as.data.frame(data_temp, stringsAsFactors=FALSE)
  
  assign(x[[1]], data_temp, envir = .GlobalEnv)
  
  gc()
  
},unknowns=unknowns_strings, .progress = "text")

rm2(unknowns_all)


###############################################################################
cat("SECTION: CHANGE NOT SPECIFIED PHRASES TO NA", "\n")
###############################################################################

NA_Phrases <- c("NA","N/A","N\\A","NOT APPLICABLE","NOT APPILCABLE","NOT DEFINED","NOT DISCLOSED","NOT DISLCOSED","NOT DISLOSED","UNDISCLOSED",
                "TO BE ADVISED","TO BE ADVISE","TBA","SEE PROSPECTUS FOR FULL DETAILS","UPON REQUEST",
                "SUBJECT TO MANAGER'S DISCRETION")

EurekahedgeHF_Stats_noreturns_na_phrases <- EurekahedgeHF_Stats_noreturns_unknowns
EurekahedgeHF_Fund_Details_na_phrases <- EurekahedgeHF_Fund_Details_unknowns
EurekahedgeHF_Fee_and_Redemption_na_phrases <- EurekahedgeHF_Fee_and_Redemption_unknowns
EurekahedgeHF_Profile_Strategy_na_phrases <- EurekahedgeHF_Profile_Strategy_unknowns
EurekahedgeHF_Identifiers_na_phrases <- EurekahedgeHF_Identifiers_unknowns
EurekahedgeHF_Other_na_phrases <- EurekahedgeHF_Other_unknowns

rm2(EurekahedgeHF_Stats_noreturns_unknowns,EurekahedgeHF_Fund_Details_unknowns,EurekahedgeHF_Fee_and_Redemption_unknowns)
rm2(EurekahedgeHF_Profile_Strategy_unknowns,EurekahedgeHF_Identifiers_unknowns,EurekahedgeHF_Other_unknowns)

na_phrases_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_na_phrases"),
                         col=c("Dividend_Policy","Domicile","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange","Exchange_Name",
                               "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
na_phrases_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_na_phrases"),
                         col=c("Management_Fee","Performance_Fee","Other_Fee"))

na_phrases_all <- list(na_phrases_temp1,na_phrases_temp2)

rm2(na_phrases_temp1,na_phrases_temp2)

l_ply(.data=na_phrases_all, .fun = function(x,phrases){
  
  # x <- na_phrases_all[[1]]
  # phrases <- NA_Phrases
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  data_temp <- not_specified_to_na(data_temp,col_temp,phrases)
  data_temp <- as.data.frame(data_temp,stringsAsFactors=FALSE)
  
  assign(x[[1]], data_temp, envir = .GlobalEnv)
  
  gc()
  
},phrases=NA_Phrases, .progress = "text")

rm2(na_phrases_all)


###############################################################################
cat("SECTION: CHANGE NO PHRASES TO 'NO'", "\n")
###############################################################################

NO_Phrases <- c("NIL","NONE","NONE AFTER 12 MONTHS","NONE AFTER 1ST YEAR","NO DIVIDEND","NON DIVIDEND","LITTLE OR NO")

EurekahedgeHF_Stats_noreturns_no_phrases <- EurekahedgeHF_Stats_noreturns_na_phrases
EurekahedgeHF_Fund_Details_no_phrases <- EurekahedgeHF_Fund_Details_na_phrases
EurekahedgeHF_Fee_and_Redemption_no_phrases <- EurekahedgeHF_Fee_and_Redemption_na_phrases
EurekahedgeHF_Profile_Strategy_no_phrases <- EurekahedgeHF_Profile_Strategy_na_phrases
EurekahedgeHF_Identifiers_no_phrases <- EurekahedgeHF_Identifiers_na_phrases
EurekahedgeHF_Other_no_phrases <- EurekahedgeHF_Other_na_phrases

rm2(EurekahedgeHF_Stats_noreturns_na_phrases,EurekahedgeHF_Fund_Details_na_phrases,EurekahedgeHF_Fee_and_Redemption_na_phrases)
rm2(EurekahedgeHF_Profile_Strategy_na_phrases,EurekahedgeHF_Identifiers_na_phrases,EurekahedgeHF_Other_na_phrases)

no_phrases_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_no_phrases"),
                         col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange",
                               "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
no_phrases_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_no_phrases"),
                         col=c("Management_Fee","Performance_Fee","Other_Fee"))

no_phrases_all <- list(no_phrases_temp1,no_phrases_temp2)

rm2(no_phrases_temp1,no_phrases_temp2)

l_ply(.data=no_phrases_all, .fun = function(x,phrases){
  
  # x <- no_phrases_all[[1]]
  # phrases <- NO_Phrases
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  data_temp <- not_specified_to_na(data_temp,col_temp,phrases)
  data_temp <- as.data.frame(data_temp,stringsAsFactors=FALSE)
  
  assign(x[[1]], data_temp, envir = .GlobalEnv)
  
  gc()
  
},phrases=NO_Phrases, .progress = "text")

rm2(no_phrases_all)


###############################################################################
cat("SECTION: CHANGE YES PHRASES TO 'YES'", "\n")
###############################################################################

YES_Phrases <- c("RARELY","OCCASIONALLY")

EurekahedgeHF_Stats_noreturns_yes_phrases <- EurekahedgeHF_Stats_noreturns_no_phrases
EurekahedgeHF_Fund_Details_yes_phrases <- EurekahedgeHF_Fund_Details_no_phrases
EurekahedgeHF_Fee_and_Redemption_yes_phrases <- EurekahedgeHF_Fee_and_Redemption_no_phrases
EurekahedgeHF_Profile_Strategy_yes_phrases <- EurekahedgeHF_Profile_Strategy_no_phrases
EurekahedgeHF_Identifiers_yes_phrases <- EurekahedgeHF_Identifiers_no_phrases
EurekahedgeHF_Other_yes_phrases <- EurekahedgeHF_Other_no_phrases

rm2(EurekahedgeHF_Stats_noreturns_no_phrases,EurekahedgeHF_Fund_Details_no_phrases,EurekahedgeHF_Fee_and_Redemption_no_phrases)
rm2(EurekahedgeHF_Profile_Strategy_no_phrases,EurekahedgeHF_Identifiers_no_phrases,EurekahedgeHF_Other_no_phrases)

yes_phrases_temp1 <- list(data=c("EurekahedgeHF_Fund_Details_yes_phrases"),
                          col=c("Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange",
                                "Fund_Size_USm","Fund_Capacity_USm","Firms_Total_Asset_USm","Total_Asset_in_Hedge_Funds_USm"))
yes_phrases_temp2 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_yes_phrases"),
                          col=c("Management_Fee","Performance_Fee","Other_Fee"))

yes_phrases_all <- list(yes_phrases_temp1,yes_phrases_temp2)

rm2(yes_phrases_temp1,yes_phrases_temp2)

l_ply(.data=yes_phrases_all, .fun = function(x,phrases){
  
  # x <- yes_phrases_all[[1]]
  # phrases <- NO_Phrases
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  data_temp <- not_specified_to_na(data_temp,col_temp,phrases)
  data_temp <- as.data.frame(data_temp,stringsAsFactors=FALSE)
  
  assign(x[[1]], data_temp, envir = .GlobalEnv)
  
  gc()
  
},phrases=YES_Phrases, .progress = "text")

rm2(yes_phrases_all)


###############################################################################
cat("SECTION: CONVERT FEES TO NUMERIC", "\n")
###############################################################################

EurekahedgeHF_Stats_noreturns_fees_numeric <- EurekahedgeHF_Stats_noreturns_yes_phrases
EurekahedgeHF_Fund_Details_fees_numeric <- EurekahedgeHF_Fund_Details_yes_phrases
EurekahedgeHF_Fee_and_Redemption_fees_numeric <- EurekahedgeHF_Fee_and_Redemption_yes_phrases
EurekahedgeHF_Profile_Strategy_fees_numeric <- EurekahedgeHF_Profile_Strategy_yes_phrases
EurekahedgeHF_Identifiers_fees_numeric <- EurekahedgeHF_Identifiers_yes_phrases
EurekahedgeHF_Other_fees_numeric <- EurekahedgeHF_Other_yes_phrases

rm2(EurekahedgeHF_Stats_noreturns_yes_phrases,EurekahedgeHF_Fund_Details_yes_phrases,EurekahedgeHF_Fee_and_Redemption_yes_phrases)
rm2(EurekahedgeHF_Profile_Strategy_yes_phrases,EurekahedgeHF_Identifiers_yes_phrases,EurekahedgeHF_Other_yes_phrases)

fee_numeric_temp1 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_fees_numeric"),
                          col=c("Management_Fee","Performance_Fee","Other_Fee"))

fee_numeric_all <- list(fee_numeric_temp1)

rm2(fee_numeric_temp1)

l_ply(.data=fee_numeric_all, .fun = function(x,unknowns){
  
  # x <- fee_numeric_all[[1]]
  # unknowns <- unknowns_strings
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  for (i in 1:length(col_temp)) {
    
    # i <- 1
    # i <- 2
    # i <- 3
    
    data_temp_convert <-  data.frame(matrix(NA, ncol=7, nrow=nrow(data_temp), 
                                            dimnames=list(c(), c("org","no_char","no_punct","first_num","remove_blanks",
                                                                 "flag","convert"))),stringsAsFactors=FALSE)
    
    data_temp_convert[,"org"] <- data_temp[,col_temp[i]]
    data_temp_convert[,"org"] <- gsub(pattern=" {2,}", replacement=" ", x=data_temp_convert[,"org"])
    data_temp_convert[,"org"] <- gsub("^\\s+|\\s+$", "", data_temp_convert[,"org"])
    
    #data_temp_convert <- unique(data_temp_convert)
    
    #Remove Characters
    data_temp_convert[,"no_char"] <- data_temp_convert[,"org"]
    data_temp_convert[,"no_char"] <- gsub(pattern="([[:alpha:]])", replacement="", x=data_temp_convert[,"no_char"])
    data_temp_convert[,"no_char"] <- gsub(pattern=" {2,}", replacement=" ", x=data_temp_convert[,"no_char"])
    data_temp_convert[,"no_char"] <- gsub("^\\s+|\\s+$", "", data_temp_convert[,"no_char"])
    
    #Remove Punctuation (except % and $)
    data_temp_convert[,"no_punct"] <- data_temp_convert[,"no_char"]
    #data_temp_convert[,"no_punct"] <- gsub(pattern="([[:punct:]])", replacement="", x=data_temp_convert[,"no_punct"])
    data_temp_convert[,"no_punct"] <- gsub(pattern="(:|,|;|=|-|&|\\+|_|?|!|/)", replacement="", x=data_temp_convert[,"no_punct"])
    data_temp_convert[,"no_punct"] <- gsub(pattern="\\. ", replacement="\\.", x=data_temp_convert[,"no_punct"])
    data_temp_convert[,"no_punct"] <- gsub(pattern=" \\. ", replacement="\\.", x=data_temp_convert[,"no_punct"])
    data_temp_convert[,"no_punct"] <- gsub(pattern="\\.{2,}", replacement="\\.", x=data_temp_convert[,"no_punct"])
    data_temp_convert[,"no_punct"] <- gsub(pattern=" {2,}", replacement=" ", x=data_temp_convert[,"no_punct"])
    data_temp_convert[,"no_punct"] <- gsub("^\\s+|\\s+$", "", data_temp_convert[,"no_punct"])
    
    #Get first number
    data_temp_convert[,"first_num"] <- data_temp_convert[,"no_punct"]
    data_temp_convert[,"first_num"] <- gsub(pattern=" .*$", replacement="", x=data_temp_convert[,"first_num"])
    data_temp_convert[,"first_num"] <- gsub(pattern=" {2,}", replacement=" ", x=data_temp_convert[,"first_num"])
    data_temp_convert[,"first_num"] <- gsub("^\\s+|\\s+$", "", data_temp_convert[,"first_num"])
    
    #Remove Blanks
    data_temp_convert[,"remove_blanks"] <- data_temp_convert[,"first_num"]
    data_temp_convert[,"remove_blanks"] <- gsub(pattern="(\\$|%)", replacement="", x=data_temp_convert[,"remove_blanks"])
    data_temp_convert[,"remove_blanks"] <- gsub(pattern=" {2,}", replacement=" ", x=data_temp_convert[,"remove_blanks"])
    data_temp_convert[,"remove_blanks"] <- gsub("^\\s+|\\s+$", "", data_temp_convert[,"remove_blanks"])
    data_temp_convert[,"remove_blanks"] <- unknownToNA(data_temp_convert[,"remove_blanks"], unknown=unknowns,force=TRUE)
    data_temp_convert[,"remove_blanks"] <- ifelse(is.na(data_temp_convert[,"remove_blanks"]),NA,data_temp_convert[,"remove_blanks"])
    data_temp_convert[,"remove_blanks"] <- as.numeric(data_temp_convert[,"remove_blanks"])
    
    #Flag number type
    # decimal: flag=1
    # percent: flag=2
    # bp:      flag=3
    # fixed:   flag=4
    
    data_temp_convert[,"flag"] <- ifelse((data_temp_convert[,"remove_blanks"]>=0 & data_temp_convert[,"remove_blanks"]<1.0),1,
                                         ifelse((data_temp_convert[,"remove_blanks"]>=1.0 & data_temp_convert[,"remove_blanks"]<100.0),2,
                                                ifelse((data_temp_convert[,"remove_blanks"]>=100.0 & data_temp_convert[,"remove_blanks"]<10000),3,
                                                       ifelse((data_temp_convert[,"remove_blanks"]>=10000),4,NA))))
    
    data_temp_convert[,"flag"] <- ifelse(grepl("%",data_temp_convert[,"first_num"]),2,data_temp_convert[,"flag"])
    data_temp_convert[,"flag"] <- ifelse(grepl("\\$",data_temp_convert[,"first_num"]),4,data_temp_convert[,"flag"])
    
    
    #Convert to decimals
    data_temp_convert[,"convert"] <- data_temp_convert[,"remove_blanks"]
    data_temp_convert[,"convert"] <- ifelse(data_temp_convert[,"flag"]==2,data_temp_convert[,"convert"]/100,data_temp_convert[,"convert"])
    data_temp_convert[,"convert"] <- ifelse(data_temp_convert[,"flag"]==3,data_temp_convert[,"convert"]/10000,data_temp_convert[,"convert"])
    data_temp_convert[,"convert"] <- ifelse(data_temp_convert[,"flag"]==4,NA,data_temp_convert[,"convert"])
    
    data_temp[,col_temp[i]] <- data_temp_convert[,"convert"]
    
    rm(data_temp_convert)
  }
  
  assign(x[[1]], data_temp, envir = .GlobalEnv)
  
  gc()
  
},unknowns=unknowns_strings, .progress = "text")

rm2(fee_numeric_all)


###############################################################################
cat("SECTION: CHANGE Y/N TO BINARY", "\n")
###############################################################################

EurekahedgeHF_Stats_noreturns_yn_binary <- EurekahedgeHF_Stats_noreturns_fees_numeric
EurekahedgeHF_Fund_Details_yn_binary <- EurekahedgeHF_Fund_Details_fees_numeric
EurekahedgeHF_Fee_and_Redemption_yn_binary <- EurekahedgeHF_Fee_and_Redemption_fees_numeric
EurekahedgeHF_Profile_Strategy_yn_binary <- EurekahedgeHF_Profile_Strategy_fees_numeric
EurekahedgeHF_Identifiers_yn_binary <- EurekahedgeHF_Identifiers_fees_numeric
EurekahedgeHF_Other_yn_binary <- EurekahedgeHF_Other_fees_numeric

rm2(EurekahedgeHF_Stats_noreturns_fees_numeric,EurekahedgeHF_Fund_Details_fees_numeric,EurekahedgeHF_Fee_and_Redemption_fees_numeric)
rm2(EurekahedgeHF_Profile_Strategy_fees_numeric,EurekahedgeHF_Identifiers_fees_numeric,EurekahedgeHF_Other_fees_numeric)

yn_binary_temp1 <- list(data=c("EurekahedgeHF_Stats_noreturns_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp2 <- list(data=c("EurekahedgeHF_Fund_Details_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead",
                              "Invest_In_Private_Placements","Managed_Accounts_Offered","UCITS_combcol",
                              "Dividend_Policy","Fund_Closed","High_Water_Mark","Hurdle_Rate","Listed_on_Exchange"))
yn_binary_temp3 <- list(data=c("EurekahedgeHF_Fee_and_Redemption_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp4 <- list(data=c("EurekahedgeHF_Profile_Strategy_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp5 <- list(data=c("EurekahedgeHF_Identifiers_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))
yn_binary_temp6 <- list(data=c("EurekahedgeHF_Other_yn_binary"),
                        col=c("Flagship","Closed","Limited","Dead"))

yn_binary_all <- list(yn_binary_temp1,yn_binary_temp2,yn_binary_temp3,
                      yn_binary_temp4,yn_binary_temp5,yn_binary_temp6)

rm2(yn_binary_temp1,yn_binary_temp2,yn_binary_temp3)
rm2(yn_binary_temp4,yn_binary_temp5,yn_binary_temp6)


### Fix Dividend Policy

yn_binary_dividend_policy <-  data.frame(matrix(NA, ncol=3, nrow=nrow(EurekahedgeHF_Fund_Details_yn_binary), 
                                                dimnames=list(c(), c("org","replace","final"))),stringsAsFactors=FALSE)
yn_binary_dividend_policy[,"org"] <- EurekahedgeHF_Fund_Details_yn_binary[,"Dividend_Policy"]

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

EurekahedgeHF_Fund_Details_yn_binary[,"Dividend_Policy"] <- yn_binary_dividend_policy[,"final"]

rm2(yn_binary_dividend_policy)


### Fix Fund Closed

yn_binary_fund_closed <-  data.frame(matrix(NA, ncol=3, nrow=nrow(EurekahedgeHF_Fund_Details_yn_binary), 
                                            dimnames=list(c(), c("org","replace","final"))),stringsAsFactors=FALSE)
yn_binary_fund_closed[,"org"] <- EurekahedgeHF_Fund_Details_yn_binary[,"Fund_Closed"]

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

EurekahedgeHF_Fund_Details_yn_binary[,"Fund_Closed"] <- yn_binary_fund_closed[,"final"]

rm2(yn_binary_fund_closed)


### Fix High Water Mark

yn_binary_high_water_mark <-  data.frame(matrix(NA, ncol=3, nrow=nrow(EurekahedgeHF_Fund_Details_yn_binary), 
                                                dimnames=list(c(), c("org","replace","final"))),stringsAsFactors=FALSE)
yn_binary_high_water_mark[,"org"] <- EurekahedgeHF_Fund_Details_yn_binary[,"High_Water_Mark"]

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

EurekahedgeHF_Fund_Details_yn_binary[,"High_Water_Mark"] <- yn_binary_high_water_mark[,"final"]

rm2(yn_binary_high_water_mark)


l_ply(.data=yn_binary_all, .fun = function(x,unknowns){
  
  # x <- yn_binary_all[[1]]
  # x <- yn_binary_all[[2]]
  # unknowns <- unknowns_strings
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  # aa <- unique(data_temp[,c(col_temp)])
  # bb <- data_temp[data_temp[,"Fund_Closed"] %in% "Hard Closed",]
  # unique(data_temp[,c(col_temp[1])])
  # unique(data_temp[,c(col_temp[2])])
  # unique(data_temp[,c(col_temp[3])])
  # unique(data_temp[,c(col_temp[4])])
  # unique(data_temp[,c(col_temp[5])])
  # unique(data_temp[,c(col_temp[6])])
  # unique(data_temp[,c(col_temp[7])])
  # unique(data_temp[,c(col_temp[8])])
  # unique(data_temp[,c(col_temp[9])])
  # unique(data_temp[,c(col_temp[10])])
  # unique(data_temp[,c(col_temp[11])])
  # unique(data_temp[,c(col_temp[12])])
  
  col_order <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order[,"order"] <- seq(1,nrow(col_order))
  
  col_order_bin <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order_bin[,"order"] <- seq(1,nrow(col_order_bin))
  col_order_bin[,"col"] <- ifelse(col_order_bin[,"col"] %in% col_temp,paste(col_order_bin[,"col"],"_bin",sep=""),col_order_bin[,"col"])
  
  col_order_all <- unique(rbind(col_order_bin,col_order))
  
  col_order_all <- col_order_all[order(col_order_all[,"order"]),]
  row.names(col_order_all) <- seq(nrow(col_order_all))
  
  bin_cols <- paste(col_temp,"_bin",sep="")
  
  data_temp2 <-  data.frame(data_temp, matrix(NA, ncol=length(bin_cols), nrow=nrow(data_temp), dimnames=list(c(), bin_cols)), stringsAsFactors=FALSE)
  
  data_temp2[,bin_cols] <-  data_temp2[,col_temp]
  
  #Get column order
  
  #data_temp2 <- data_temp2[,sort(colnames(data_temp2), decreasing = FALSE)]
  data_temp2 <- data_temp2[,col_order_all[,"col"]]
  
  data_temp2 <- yn_to_binary(data_temp2,bin_cols)
  data_temp2 <- as.data.frame(data_temp2,stringsAsFactors=FALSE)
  
  data_temp2 <- data.table(data_temp2)[, (bin_cols) := llply(.SD, vector_clean_na,unknowns=unknowns,.progress = "text"), .SDcols = bin_cols]
  data_temp2 <- as.data.frame(data_temp2,stringsAsFactors=FALSE)
  
  data_temp2 <- data_temp2[rowSums(is.na(data_temp2[,1:ncol(data_temp2)]))<ncol(data_temp2),]
  
  row.names(data_temp2) <- seq(nrow(data_temp2))
  
  assign(x[[1]], data_temp2, envir = .GlobalEnv)
  
  gc()
  
},unknowns=unknowns_strings, .progress = "text")

rm2(yn_binary_all)


###############################################################################
cat("SECTION: SPLIT (EXPAND) COLUMNS", "\n")
###############################################################################

EurekahedgeHF_Stats_noreturns_split_expand <- EurekahedgeHF_Stats_noreturns_yn_binary
EurekahedgeHF_Fund_Details_split_expand <- EurekahedgeHF_Fund_Details_yn_binary
EurekahedgeHF_Fee_and_Redemption_split_expand <- EurekahedgeHF_Fee_and_Redemption_yn_binary
EurekahedgeHF_Profile_Strategy_split_expand <- EurekahedgeHF_Profile_Strategy_yn_binary
EurekahedgeHF_Identifiers_split_expand <- EurekahedgeHF_Identifiers_yn_binary
EurekahedgeHF_Other_split_expand <- EurekahedgeHF_Other_yn_binary

rm2(EurekahedgeHF_Stats_noreturns_yn_binary,EurekahedgeHF_Fund_Details_yn_binary,EurekahedgeHF_Fee_and_Redemption_yn_binary)
rm2(EurekahedgeHF_Profile_Strategy_yn_binary,EurekahedgeHF_Identifiers_yn_binary,EurekahedgeHF_Other_yn_binary)

split_expand_temp1 <- list(data=c("EurekahedgeHF_Other_split_expand"),
                           col=c("Custodian","Principal_Prime_Broker_combcol","Secondary_Prime_Broker","Synthetic_Prime_Broker",
                                 "Legal_Advisor_Offshore","Legal_Advisor_Onshore","Legal_Advisor"))

split_expand_all <- list(split_expand_temp1)

rm2(split_expand_temp1)

l_ply(.data=split_expand_all, .fun = function(x){
  
  # x <- split_expand_all[[1]]
  
  data_temp <- get(x[[1]])
  col_temp <- x[[2]]
  
  col_order <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order[,"order"] <- seq(1,nrow(col_order))
  
  col_order_org <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order_org[,"order"] <- seq(1,nrow(col_order_org))
  col_order_org[,"col"] <- ifelse(col_order_org[,"col"] %in% col_temp,paste(col_order_org[,"col"],"_org",sep=""),col_order_org[,"col"])
  
  col_order_comments <- data.frame(col=colnames(data_temp),order=NA,stringsAsFactors=FALSE)
  col_order_comments[,"order"] <- seq(1,nrow(col_order_comments))
  col_order_comments[,"col"] <- ifelse(col_order_comments[,"col"] %in% col_temp,paste(col_order_comments[,"col"],"_comments",sep=""),col_order_comments[,"col"])
  
  col_order_all <- unique(rbind(rbind(col_order_org,col_order),col_order_comments))
  
  col_order_all <- col_order_all[order(col_order_all[,"order"]),]
  row.names(col_order_all) <- seq(nrow(col_order_all))
  
  #Rename original columns
  data_temp <- rename.vars(data_temp, col_temp, paste(col_temp,"_org",sep=""),info=FALSE)
  
  strip_cols <- c(col_temp, paste(col_temp,"_comments",sep=""))
  
  data_temp2 <-  data.frame(data_temp, matrix(NA, ncol=length(strip_cols), nrow=nrow(data_temp), dimnames=list(c(), strip_cols)), stringsAsFactors=FALSE)
  
  #Get column order
  
  #data_temp2 <- data_temp2[,sort(colnames(data_temp2), decreasing = FALSE)]
  data_temp2 <- data_temp2[,col_order_all[,"col"]]
  
  data_temp2 <- strip_comments(data_temp2,col_temp)
  data_temp2 <- as.data.frame(data_temp2,stringsAsFactors=FALSE)
  
  #Get text before comments
  data_temp2 <- create_noncomments(data_temp2,col_temp)
  data_temp2 <- as.data.frame(data_temp2,stringsAsFactors=FALSE)
  
  assign(x[[1]], data_temp2, envir = .GlobalEnv)
  
  gc()
  
}, .progress = "text")

rm2(strip_comments_all)



