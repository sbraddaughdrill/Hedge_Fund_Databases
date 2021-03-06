# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Expand_Fields_part3.R
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
final_folder_expand1_path <- paste(output_directory, "Final_Expand1", sep = "//", collapse = "//")  
create_directory(final_folder_expand1_path,remove=1)

final_folder_expand2_path <- paste(output_directory, "Final_Expand2", sep = "//", collapse = "//")  
create_directory(final_folder_expand2_path,remove=1)

final_folder_expand3_path <- paste(output_directory, "Final_Expand3", sep = "//", collapse = "//")  
create_directory(final_folder_expand3_path,remove=1)

### Final Files
final_folder_files0 <- data.frame(files=list.files(path=final_folder_path),file_name=NA,import=NA,stringsAsFactors=FALSE)

final_folder_files <- final_folder_files0[grepl(".CSV|.csv", final_folder_files0[,"files"]),]

rm2(final_folder_files0)

final_folder_files[,"file_name"] <- final_folder_files[,"files"]
final_folder_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=final_folder_files[,"file_name"])

final_folder_files[,"import"] <- ifelse(grepl("(Stats|Fund_Detail|Fee_and_Redemption|Profile_Strategy|Identifier)",final_folder_files[,"file_name"]),1,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(NAV_AUM_Ret)",final_folder_files[,"file_name"]),2,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(Other)",final_folder_files[,"file_name"]),3,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(Instruments_Traded)",final_folder_files[,"file_name"]),4,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(is.na(final_folder_files[,"import"]),0,final_folder_files[,"import"])
final_folder_files[,"import"] <- ifelse(grepl("(_part2)",final_folder_files[,"file_name"]),0,final_folder_files[,"import"])

invisible(gc(verbose = FALSE, reset = TRUE))


# ### Expand Final 1 Files
# final_folder_expand1_files0 <- data.frame(files=list.files(path=final_folder_expand1_path),file_name=NA,import=NA,stringsAsFactors=FALSE)
# 
# final_folder_expand1_files <- final_folder_expand1_files0[grepl(".CSV|.csv", final_folder_expand1_files0[,"files"]),]
# 
# rm2(final_folder_expand1_files0)
# 
# final_folder_expand1_files[,"file_name"] <- final_folder_expand1_files[,"files"]
# final_folder_expand1_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=final_folder_expand1_files[,"file_name"])
# 
# final_folder_expand1_files[,"import"] <- ifelse(grepl("(Stats|Fund_Detail|Fee_and_Redemption|Profile_Strategy|Identifier)",final_folder_expand1_files[,"file_name"]),1,final_folder_expand1_files[,"import"])
# final_folder_expand1_files[,"import"] <- ifelse(grepl("(NAV_AUM_Ret|Instruments_Traded)",final_folder_expand1_files[,"file_name"]),2,final_folder_expand1_files[,"import"])
# final_folder_expand1_files[,"import"] <- ifelse(grepl("(Other)",final_folder_expand1_files[,"file_name"]),3,final_folder_expand1_files[,"import"])
# final_folder_expand1_files[,"import"] <- ifelse(is.na(final_folder_expand1_files[,"import"]),0,final_folder_expand1_files[,"import"])
# final_folder_expand1_files[,"import"] <- ifelse(grepl("(_part2)",final_folder_expand1_files[,"file_name"]),0,final_folder_expand1_files[,"import"])
# 
# invisible(gc(verbose = FALSE, reset = TRUE))


### Expand Final 2 Files
final_folder_expand2_files0 <- data.frame(files=list.files(path=final_folder_expand2_path),file_name=NA,import=NA,stringsAsFactors=FALSE)

final_folder_expand2_files <- final_folder_expand2_files0[grepl(".CSV|.csv", final_folder_expand2_files0[,"files"]),]

rm2(final_folder_expand2_files0)

final_folder_expand2_files[,"file_name"] <- final_folder_expand2_files[,"files"]
final_folder_expand2_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=final_folder_expand2_files[,"file_name"])

final_folder_expand2_files[,"import"] <- ifelse(grepl("(Stats|Fund_Detail|Fee_and_Redemption|Profile_Strategy|Identifier)",final_folder_expand2_files[,"file_name"]),1,final_folder_expand2_files[,"import"])
final_folder_expand2_files[,"import"] <- ifelse(grepl("(NAV_AUM_Ret)",final_folder_expand2_files[,"file_name"]),2,final_folder_expand2_files[,"import"])
final_folder_expand2_files[,"import"] <- ifelse(grepl("(Other)",final_folder_expand2_files[,"file_name"]),3,final_folder_expand2_files[,"import"])
final_folder_expand2_files[,"import"] <- ifelse(grepl("(Instruments_Traded)",final_folder_expand2_files[,"file_name"]),4,final_folder_expand2_files[,"import"])
final_folder_expand2_files[,"import"] <- ifelse(is.na(final_folder_expand2_files[,"import"]),0,final_folder_expand2_files[,"import"])
#final_folder_expand2_files[,"import"] <- ifelse(grepl("(_part2)",final_folder_expand2_files[,"file_name"]),0,final_folder_expand2_files[,"import"])

invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: IMPORT FILES", "\n")
###############################################################################

a_ply(.data=final_folder_expand2_files[final_folder_expand2_files[,"import"] %in% c(1),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_expand2_files[1,]
  # x <- final_folder_expand2_files[3,]
  # x <- final_folder_expand2_files[4,] 
  
  # directory_in <- final_folder_expand2_path
  # unknowns <- unknowns_strings
  
  #temp_cols <- c("date","yr","month","bad_min","bad_max")
  
  input <- data.table(read.csv(file=paste(final_folder_expand2_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE))
  
  setkeyv(input,NULL)
  
  cols <- c("pull_trim")
  for (k in cols) {
    set(input, i=NULL, j=k, value=as.character(input[[k]]))
  }
  rm(k,cols)
  
  order_ids <- c("pull_trim","pull_trim2","pull")
  order_ids_trim <- order_ids[order_ids %in% colnames(input)]
  setcolorder(input, c(order_ids_trim,colnames(input)[!(colnames(input) %in% c(order_ids_trim))]))
  
  rm(order_ids,order_ids_trim)
  
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
  
}, directory_in=final_folder_expand2_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))

a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(2),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[6,]
  
  # directory_in <- final_folder_path
  # unknowns <- unknowns_strings
  
  input_cols_keep <- c("pull_trim","pull","Fund_ID","Dead_Date","yr","month","date","bad_min","bad_max","AUM")  
  
  input <- data.table(read.csv(file=paste(final_folder_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)[,input_cols_keep])
  setkeyv(input,NULL)
  
  cols <- c("pull_trim")
  for (k in cols) {
    set(input, i=NULL, j=k, value=as.character(input[[k]]))
  }
  rm(k,cols)
  
  setnames(input,"AUM","pull_trim2")
  
  input <- input[, pull_trim2:=pull,by=NULL]
  
  cols <- c("pull_trim2")
  for (k in cols) {
    set(input, i=NULL, j=k, value=gsub(pattern="_NAV_AUM", replacement="", input[[k]], perl=TRUE))
  }
  rm(k,cols)
  
  order_ids <- c("pull_trim","pull_trim2","pull")
  order_ids_trim <- order_ids[order_ids %in% colnames(input)]
  setcolorder(input, c(order_ids_trim,colnames(input)[!(colnames(input) %in% c(order_ids_trim))]))
  
  rm(order_ids,order_ids_trim)
  
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
  input <- input[(droprow)][,droprow:=NULL][]
  
  input <- input[,which(unlist(lapply(input, function(x)!all(is.na(x))))),with=FALSE]
  
  setorderv(input, c("Fund_ID","pull_trim","pull_trim2"),c(1,1,1))
  
  if("Date_Added" %in% colnames(input)){input <- input[, Date_Added:=as.yearmon(Date_Added,format="%b %Y")]} 
  if("Dead_Date" %in% colnames(input)){input <- input[, Dead_Date:=as.yearmon(Dead_Date,format="%b %Y")]} 
  if("Inception_Date" %in% colnames(input)){input <- input[, Inception_Date:=as.yearmon(Inception_Date,format="%b %Y")]} 
  
  assign("Merge_IDs", input, envir = .GlobalEnv)
  
  rm(input)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
}, directory_in=final_folder_path, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

invisible(gc(verbose = FALSE, reset = TRUE))

a_ply(.data=final_folder_files[final_folder_files[,"import"] %in% c(4),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_files[5,]
  
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


###############################################################################
cat("SECTION: ADD PULL IDS AND DEAD DATES", "\n")
###############################################################################

fix_dead_dates_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_part2"),data_out=c("EurekahedgeHF_Stats_noreturns_part3"))
fix_dead_dates_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_part2"),data_out=c("EurekahedgeHF_Fund_Details_part3"))
fix_dead_dates_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_part2"),data_out=c("EurekahedgeHF_Fee_and_Redemption_part3"))
fix_dead_dates_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_part2"),data_out=c("EurekahedgeHF_Profile_Strategy_part3"))
fix_dead_dates_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_part2"),data_out=c("EurekahedgeHF_Identifiers_part3"))
fix_dead_dates_temp6 <- list(data_in=c("EurekahedgeHF_Instruments_Traded"),data_out=c("EurekahedgeHF_Instruments_Traded_part3"))

fix_dead_dates_all0 <- list(fix_dead_dates_temp1,fix_dead_dates_temp2,fix_dead_dates_temp3,
                            fix_dead_dates_temp4,fix_dead_dates_temp5,fix_dead_dates_temp6)

rm2(fix_dead_dates_temp1,fix_dead_dates_temp2,fix_dead_dates_temp3)
rm2(fix_dead_dates_temp4,fix_dead_dates_temp5,fix_dead_dates_temp6)

fix_dead_dates_all1 <- llply(.data=fix_dead_dates_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(fix_dead_dates_all0)

fix_dead_dates_all1 <- fix_dead_dates_all1[order(-sapply(fix_dead_dates_all1,"[[","size"),
                                             -sapply(fix_dead_dates_all1,"[[","ncol"),
                                             -sapply(fix_dead_dates_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))

Merge_IDs[,c("pull")] <- NULL

invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=fix_dead_dates_all1, .fun = function(x,ids,merge_data){
  
  # x <- fix_dead_dates_all1[[1]]
  # x <- fix_dead_dates_all1[[2]]
  # x <- fix_dead_dates_all1[[3]]
  # x <- fix_dead_dates_all1[[4]]
  # x <- fix_dead_dates_all1[[5]]
  # x <- fix_dead_dates_all1[[6]]
  
  # ids <- c("pull_trim","pull_trim2")
  # merge_data <- "Merge_IDs"
  
  require(data.table)
  
  #data_name_temp <- x[[which(names(x)==c("data_in"))]]
  #data_temp <- get(x[[which(names(x)==c("data_in"))]])
  
  cat(x[[which(names(x)==c("data_in"))]], "\n")
  
  files_temp1 <- x[[which(names(x)==c("data_in"))]]
  files_temp1_trim <- data.frame(file=x[[which(names(x)==c("data_in"))]],row_str=NA,col_str=NA,file_str=NA,stringsAsFactors=FALSE)
  #files_temp1_trim[,"row_str"] <- paste(files_temp1_trim[,"file"],"[",",","'pull_trim2'","]","==","'",x[,"pull_trim2"],"'",sep="")
  files_temp1_trim[,"row_str"] <- ""
  #files_temp1_trim[,"col_str"] <- paste("!(colnames(",files_temp1_trim[,"file"],") %in% c('Dead_Date'))",sep="")
  files_temp1_trim[,"col_str"] <- ""
  files_temp1_trim[,"file_str"] <- paste(files_temp1_trim[,"file"],"[",files_temp1_trim[,"row_str"],",",files_temp1_trim[,"col_str"],"]",sep="")
  
  rm(files_temp1)
  
  files_temp2 <- merge_data
  files_temp2_trim <- data.frame(file=files_temp2[!is.na(files_temp2)],row_str=NA,col_str=NA,file_str=NA,stringsAsFactors=FALSE)
  #files_temp2_trim[,"row_str"] <- paste(files_temp2_trim[,"file"],"[",",","'pull_trim2'","]","==","'",x[,"pull_trim2"],"'",sep="")
  files_temp2_trim[,"row_str"] <- ""
  #files_temp2_trim[,"col_str"] <- paste("!(colnames(",files_temp2_trim[,"file"],") %in% c('pull'))",sep="")
  files_temp2_trim[,"col_str"] <- ""
  files_temp2_trim[,"file_str"] <- paste(files_temp2_trim[,"file"],"[",files_temp2_trim[,"row_str"],",",files_temp2_trim[,"col_str"],"]",sep="")
  
  rm(files_temp2)
  
  files_temp_all <- rbind(files_temp2_trim,files_temp1_trim)
  
  rm(files_temp2_trim,files_temp1_trim)
  
  #merge_temp <- eval(parse(text=files_temp_all[1,"file_str"]))
  
  #for (i in 1:nrow(files_temp_all)) {print(paste(files_temp_all[i,"file"],": nrow = ",nrow(eval(parse(text=files_temp_all[i,"file_str"]))),sep=""))}
  
  #for (i in 2:nrow(files_temp_all)) {
  
  # i <- 2
  
  #data_temp <- data.frame(file_str=rbind(files_temp_all[1,"file_str"],files_temp_all[2,"file_str"]),stringsAsFactors=FALSE)
  
  #Merge_IDs[,!(colnames(Merge_IDs) %in% c('pull'))]
  
  if("Dead_Date" %in% colnames(get(files_temp_all[2,"file"]))){eval(parse(text=paste("set(",files_temp_all[2,"file"],",j=which(colnames(",files_temp_all[2,"file"],") %in% c('Dead_Date')), value=NULL)",sep="")))} 
  
  common_ids1a <- adply(.data=files_temp_all, .margins=1, .fun = function(x){
    
    # x <- files_temp_all[1,]
    # x <- files_temp_all[2,]
    
    temp_order <- data.frame(cols=eval(parse(text=paste("colnames(",x[,"file"],")[",x[,"col_str"],"]",sep=""))),order=NA,stringsAsFactors=FALSE)
    temp_order[,"order"] <- seq(1,nrow(temp_order))
    return(temp_order)
    
  }, .expand = FALSE, .progress = "none")
  
  common_ids1b <- ddply(.data=common_ids1a[,!(colnames(common_ids1a) %in% c("X1"))], .variables=c("cols"),.fun = function(x){
    
    return(data.frame(freq=nrow(x),avg_order=mean(x[,"order"]),stringsAsFactors=FALSE))
    
  }, .progress = "none")
  
  rm(common_ids1a)
  
  common_ids1b <- common_ids1b[order(common_ids1b[,"avg_order"],common_ids1b[,"freq"],common_ids1b[,"cols"]),]
  row.names(common_ids1b) <- seq(nrow(common_ids1b))
  
  common_ids1 <- common_ids1b[common_ids1b[,"freq"]==nrow(files_temp_all),] 
  
  rm(common_ids1b)
  
  #setkeyv(eval(parse(text=paste("",files_temp_all[1,"file"],"",sep=""))),common_ids1[,'cols'])
  #setkeyv(eval(parse(text=paste("",files_temp_all[2,"file"],"",sep=""))),common_ids1[,'cols'])
  
  eval(parse(text=paste("setkeyv(",files_temp_all[1,"file"],",common_ids1[,'cols'])",sep=""))) 
  eval(parse(text=paste("setkeyv(",files_temp_all[2,"file"],",common_ids1[,'cols'])",sep="")))
  
  #key(Merge_IDs)
  #key(EurekahedgeHF_Fund_Details)
  #key(EurekahedgeHF_Stats_noreturns)
  
  #eval(parse(text=paste(files_temp_all[2,"file"],"[,c('Dead_Date')] <- NULL",sep="")))
  #eval(parse(text=paste("set(",files_temp_all[2,"file"],",j=which(colnames(",files_temp_all[2,"file"],") %in% c('Dead_Date')), value=NULL)",sep="")))
  
  merge_temp <- merge(eval(parse(text=paste("",files_temp_all[1,"file"],"",sep=""))),
                      eval(parse(text=paste("",files_temp_all[2,"file"],"",sep=""))),
                      by.x=common_ids1[,"cols"], by.y=common_ids1[,"cols"], 
                      all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))
  
  
  rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
  rm(files_temp_all,common_ids1)
  invisible(gc(verbose = FALSE, reset = TRUE))
  
  order_ids <- c("pull_trim","pull_trim2",
                 "Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason",
                 "date","yr","month","bad_min","bad_max")
  
  order_ids_trim <- order_ids[order_ids %in% colnames(merge_temp)]
  
  rm(order_ids)
  
  setcolorder(merge_temp, c(order_ids_trim,colnames(merge_temp)[!(colnames(merge_temp) %in% c(order_ids_trim))]))
  setorderv(merge_temp, c("Fund_ID","date","pull_trim","pull_trim"),c(1,1,1,1))
  
  setkeyv(merge_temp,NULL)
  
  #assign(x[[which(names(x)==c("data_out"))]], merge_temp, envir = .GlobalEnv)
  write.csv(merge_temp, file=paste(final_folder_expand3_path,"//",x[[which(names(x)==c("data_out"))]],".csv",sep=""),row.names=FALSE)
  
  rm(merge_temp,order_ids_trim)
  
  invisible(gc(verbose = FALSE, reset = TRUE))
  
},ids=c("pull_trim","pull_trim2"), merge_data="Merge_IDs", .progress = "text")

rm2(Merge_IDs,fix_dead_dates_all1)
invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))


# ###############################################################################
# cat("SECTION: IMPORT OTHER DATA", "\n")
# ###############################################################################

invisible(gc(verbose = FALSE, reset = TRUE))

file.copy(paste(final_folder_expand2_path,"//",final_folder_expand2_files[final_folder_expand2_files[,"import"] %in% c(3),"file_name"],".csv",sep=""), 
          paste(final_folder_expand3_path,"//",gsub("_part2","_part3",final_folder_expand2_files[final_folder_expand2_files[,"import"] %in% c(3),"file_name"]),".csv",sep=""),overwrite=TRUE)

file.copy(paste(final_folder_path,"//",final_folder_files[final_folder_files[,"import"] %in% c(2),"file_name"],".csv",sep=""), 
          paste(final_folder_expand3_path,"//",gsub("_part2","_part3",final_folder_files[final_folder_files[,"import"] %in% c(2),"file_name"]),".csv",sep=""),overwrite=TRUE)


