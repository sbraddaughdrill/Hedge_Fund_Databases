# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Revisors.R
# Version: 1.0
# Date:    02.25.2015
# Purpose: Find Hedge Funds that Revise
#
###############################################################################

###############################################################################
cat("SECTION: INITIAL SETUP","\n")
###############################################################################

# Clear workspace
rm(list=ls(all=T))
rm(list=ls(all.names=T))
invisible(gc(verbose=F,reset=T))

# Limit History to not exceed 500 lines
Sys.setenv(R_HISTSIZE=500)

repo <- c("http://cran.us.r-project.org")
options(repos=structure(repo))
options(install.packages.check.source=F)

# String as factors is F -- used for read.csv
options(StringsAsFactors=F)

# Default maxprint option
options(max.print=500)
# options(max.print=99999)

# Memory limit
#memory.limit(size=8183)

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
  
  #input_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  input_directory <- normalizePath("F:/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  output_directory <- normalizePath("F:/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  #function_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)    
  function_directory <- normalizePath("F:/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)  
  
} else if (Location == 2) {
  
  input_directory <- normalizePath("C:/Users/bdaughdr/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  output_directory <- normalizePath("C:/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  function_directory <- normalizePath("C:/Users/bdaughdr/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)   
  
} else if (Location == 3) {
  
  input_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  output_directory <- normalizePath("C:/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  function_directory <- normalizePath("C:/Users/S.Brad/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)
  
} else if (Location == 4) {
  
  input_directory <- normalizePath("H:/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  #output_directory <- normalizePath("C:/Users/bdaughdr/Documents/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  output_directory <- normalizePath("H:/Research/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  #function_directory <- normalizePath("//tsclient/C/Users/S.Brad/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)
  function_directory <- normalizePath("//tsclient/F/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)
  
} else if (Location == 5) {
  
  input_directory <- normalizePath("H:/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  #output_directory <- normalizePath("C:/Users/bdaughdr/Documents/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  output_directory <- normalizePath("H:/Research/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  function_directory <- normalizePath("//tsclient/C/Users/bdaughdr/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)
  
} else if (Location == 6) {
  
  input_directory <- normalizePath("H:/Research/Hedge_Fund_Databases/Data",winslash="\\",mustWork=T)
  #output_directory <- normalizePath("C:/Users/bdaughdr/Documents/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  output_directory <- normalizePath("H:/Research/Import_Data/Data/Eurekahedge",winslash="\\",mustWork=T)
  #function_directory <- normalizePath("//tsclient/C/Users/S.Brad/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)
  function_directory <- normalizePath("//tsclient/F/Dropbox/Research_Methods/R",winslash="\\",mustWork=T)
  
} else {
  
  cat("ERROR ASSIGNING DIRECTORIES","\n")
  
}
rm(Location)


###############################################################################
cat("SECTION: FUNCTIONS","\n")
###############################################################################

#source(file=paste(function_directory,"functions_db.R",sep="\\"),echo=F)
#source(file=paste(function_directory,"functions_statistics.R",sep="\\"),echo=F)
#source(file=paste(function_directory,"functions_text_analysis.R",sep="\\"),echo=F)
#source(file=paste(function_directory,"functions_text_parse.R",sep="\\"),echo=F)
#source(file=paste(function_directory,"functions_utilities.R",sep="\\"),echo=F)

source(file=paste(function_directory,"functions_load_unload.R",sep="\\"),echo=F)

###############################################################################
# LIBRARIES;
cat("SECTION: LIBRARIES","\n")
###############################################################################

#LoadFromPackage("gdata",trim,bindData)
#LoadFromPackage("plyr",adply)
#test <- data.frame(one=c(1,2,3,4),two=c(1,2,3,4),stringsAsFactors=F)
#tets_out <- adply(.data=test,.margins=1,.fun=function(x){
#  cat(x[,"one"]+x[,"two"],"\n")
#},.expand=T)



#Load External Packages

LoadFromSource(file=paste(function_directory,"functions_utilities.R",sep="\\"),load_external_packages,list_installed_packages)

#"cwhmisc","DataCombine","formatR","gtools","Hmisc","lubridate","memisc",
#"mitools","plm","splitstackshape","stringi","stringr","taRifx","tm","zoo"

external_packages <- c("data.table","gdata","limma","plyr","reshape2","zoo")
invisible(unlist(sapply(external_packages,load_external_packages,repo_str=repo,simplify=F,USE.NAMES=F)))
installed_packages <- list_installed_packages(external_packages)

rm(external_packages,installed_packages,repo)

Unload(load_external_packages,list_installed_packages)


###############################################################################
cat("SECTION: SETUP","\n")
###############################################################################

identifier <- "Fund_ID"
Ret_col_name <- "Monthly_Ret"
AUM_col_name <- "AUM"


#final_folder_path <- paste(output_directory,"Final",sep="//",collapse="//")  
final_folder_path <- paste(output_directory,"Final_Expand3",sep="//",collapse="//")  

melt_folder_path <- paste(output_directory,"NAV_AUM_melt",sep="//",collapse="//")  


LoadFromSource(file=paste(function_directory,"functions_utilities.R",sep="\\"),create_directory)

revision_folder_path <- paste(output_directory,"Revision",sep="//",collapse="//")  
create_directory(revision_folder_path,remove=1)

Unload(create_directory)


###############################################################################
cat("SECTION: IMPORT DEAD DATES","\n")
###############################################################################

#Dead_cols_keep <- c("pull",identifier,"Dead_Date","date","min_date","max_date","bad_min","bad_max")
Dead_cols_keep <- c("pull",identifier,"Dead_Date","date")

#Dead_dates <- data.frame(read.csv(file=paste(final_folder_path,"\\","EurekahedgeHF_NAV_AUM_Ret",".csv",sep=""),header=T,na.strings="NA",stringsAsFactors=F),stringsAsFactors=F)[Dead_cols_keep]
Dead_dates <- data.frame(read.columns(file=paste(final_folder_path,"\\","EurekahedgeHF_NAV_AUM_Ret",".csv",sep=""),required.col=Dead_cols_keep,sep=",",na.strings="NA",stringsAsFactors=F))

for(i in which(sapply(Dead_dates,class)=="character"))
{
  # Dead_dates[[i]] <- trim(Dead_dates[[i]])
  Dead_dates[[i]] <- gsub("^\\s+|\\s+$","",Dead_dates[[i]],perl=T)
}
rm(i)
for (i in 1:ncol(Dead_dates))
{
  Dead_dates[,i] <- unknownToNA(Dead_dates[,i],unknown=unknowns_strings,force=T)
  Dead_dates[,i] <- ifelse(is.na(Dead_dates[,i]),NA,Dead_dates[,i])
} 
rm(i)

rm(Dead_cols_keep)

#test <- Dead_dates
Dead_dates[,"Dead_Date"] <- as.Date(Dead_dates[,"Dead_Date"],format="%Y-%m-%d")
Dead_dates[,"date"] <- as.Date(Dead_dates[,"date"],format="%Y-%m-%d")
#Dead_dates[,"min_date"] <- as.Date(Dead_dates[,"min_date"],format="%Y-%m-%d")
#Dead_dates[,"max_date"] <- as.Date(Dead_dates[,"max_date"],format="%Y-%m-%d")

Dead_dates_cols_id <- c(identifier,"date")
Dead_dates_pulls_u <- unique(Dead_dates[,"pull"])
Dead_dates_pulls_u_keep <- Dead_dates_pulls_u[!grepl("201207",Dead_dates_pulls_u)]
Dead_dates_pulls_u_drop <- Dead_dates_pulls_u[grepl("201207",Dead_dates_pulls_u)]

Dead_dates_trim <- data.frame(Dead_dates[Dead_dates[,"pull"] %in% Dead_dates_pulls_u_keep,],alive_flag=NA,stringsAsFactors=F)
Dead_dates_trim[,"alive_flag"] <- ifelse(is.na(Dead_dates_trim[,"Dead_Date"]),1,Dead_dates_trim[,"alive_flag"])
Dead_dates_trim[,"alive_flag"] <- ifelse((!is.na(Dead_dates_trim[,"Dead_Date"]) & Dead_dates_trim[,"Dead_Date"]>=Dead_dates_trim[,"date"]),2,Dead_dates_trim[,"alive_flag"])
Dead_dates_trim[,"alive_flag"] <- ifelse((!is.na(Dead_dates_trim[,"Dead_Date"]) & Dead_dates_trim[,"Dead_Date"]<Dead_dates_trim[,"date"]),0,Dead_dates_trim[,"alive_flag"])

#Dead_dates_trim_dead_date <- unique(Dead_dates_trim[,c("pull",identifier,"Dead_Date")])
Dead_dates_trim_dead_date <- unique(Dead_dates_trim[,c(identifier,"Dead_Date")])

#Dead_dates_trim_dead_date_count <- count(Dead_dates_trim_dead_date[!is.na(Dead_dates_trim_dead_date[,"Dead_Date"]),],identifier)
Dead_dates_trim_dead_date_count <- ddply(.data=Dead_dates_trim_dead_date[!is.na(Dead_dates_trim_dead_date[,"Dead_Date"]),],.variables=identifier,.fun=function(x){
  return(data.frame(x,freq=nrow(x),stringsAsFactors=F))
},.progress="none")

Dead_dates_trim_dead_date_count <- Dead_dates_trim_dead_date_count[order(Dead_dates_trim_dead_date_count[,identifier],
                                                                         Dead_dates_trim_dead_date_count[,"Dead_Date"]),] 
row.names(Dead_dates_trim_dead_date_count) <- seq(nrow(Dead_dates_trim_dead_date_count))

Dead_dates_trim_dead_date_good <- ddply(.data=Dead_dates_trim_dead_date_count,.variables=identifier,.fun=function(x){return(tail(x,1))},.progress="none")

rm(Dead_dates_trim_dead_date,Dead_dates_trim_dead_date_count)
rm(Dead_dates,Dead_dates_trim)
rm(Dead_dates_cols_id,Dead_dates_pulls_u,Dead_dates_pulls_u_keep,Dead_dates_pulls_u_drop)


###############################################################################
cat("SECTION: IMPORT ADDED DATES","\n")
###############################################################################

#Date_cols_keep <- c("pull_trim","pull_trim2",identifier,"Date_Added","date","bad_min","bad_max")
Date_cols_keep <- c("pull_trim","pull_trim2",identifier,"Date_Added","date")

#Added_dates <- data.frame(read.csv(file=paste(final_folder_path,"\\","EurekahedgeHF_Stats_noreturns_part3",".csv",sep=""),header=T,na.strings="NA",stringsAsFactors=F),stringsAsFactors=F)[Date_cols_keep]
Added_dates <- data.frame(read.columns(file=paste(final_folder_path,"\\","EurekahedgeHF_Stats_noreturns_part3",".csv",sep=""),required.col=Date_cols_keep,sep=",",na.strings="NA",stringsAsFactors=F))


for(i in which(sapply(Added_dates,class)=="character"))
{
  # Added_dates[[i]] <- trim(Added_dates[[i]])
  Added_dates[[i]] <- gsub("^\\s+|\\s+$","",Added_dates[[i]],perl=T)
}
rm(i)
for (i in 1:ncol(Added_dates))
{
  Added_dates[,i] <- unknownToNA(Added_dates[,i],unknown=unknowns_strings,force=T)
  Added_dates[,i] <- ifelse(is.na(Added_dates[,i]),NA,Added_dates[,i])
} 
rm(i)

rm(Date_cols_keep)

Added_dates[,"Date_Added"] <- as.yearmon(Added_dates[,"Date_Added"],format="%b %Y")
#Added_dates[,"Date_Added"] <- as.Date(Added_dates[,"Date_Added"],format="%m-%d-%Y")
Added_dates[,"Date_Added"] <- as.Date(Added_dates[,"Date_Added"],format="%Y-%m-%d")
Added_dates[,"date"] <- as.Date(Added_dates[,"date"],format="%Y-%m-%d")

Added_dates_cols_id <- c(identifier,"date")
Added_dates_pulls_u <- unique(Added_dates[,"pull_trim2"])
Added_dates_pulls_u_keep <- Added_dates_pulls_u[!grepl("201207",Added_dates_pulls_u)]
Added_dates_pulls_u_drop <- Added_dates_pulls_u[grepl("201207",Added_dates_pulls_u)]

Added_dates_trim <- data.frame(Added_dates[Added_dates[,"pull_trim2"] %in% Added_dates_pulls_u_keep,],alive_flag=NA,stringsAsFactors=F)
Added_dates_trim[,"alive_flag"] <- ifelse(is.na(Added_dates_trim[,"Date_Added"]),1,Added_dates_trim[,"alive_flag"])
Added_dates_trim[,"alive_flag"] <- ifelse((!is.na(Added_dates_trim[,"Date_Added"]) & Added_dates_trim[,"Date_Added"]>=Added_dates_trim[,"date"]),2,Added_dates_trim[,"alive_flag"])
Added_dates_trim[,"alive_flag"] <- ifelse((!is.na(Added_dates_trim[,"Date_Added"]) & Added_dates_trim[,"Date_Added"]<Added_dates_trim[,"date"]),0,Added_dates_trim[,"alive_flag"])

#Added_dates_trim_date_added <- unique(Added_dates_trim[,c("pull",identifier,"Date_Added")])
Added_dates_trim_date_added <- unique(Added_dates_trim[,c(identifier,"Date_Added")])

#Added_dates_trim_date_added_count <- count(Added_dates_trim_date_added[!is.na(Added_dates_trim_date_added[,"Date_Added"]),],identifier)
Added_dates_trim_date_added_count <- ddply(.data=Added_dates_trim_date_added[!is.na(Added_dates_trim_date_added[,"Date_Added"]),],.variables=identifier,.fun=function(x){
  return(data.frame(x,freq=nrow(x),stringsAsFactors=F))
},.progress="none")

Added_dates_trim_date_added_count <- Added_dates_trim_date_added_count[order(Added_dates_trim_date_added_count[,identifier],
                                                                             Added_dates_trim_date_added_count[,"Date_Added"]),] 
row.names(Added_dates_trim_date_added_count) <- seq(nrow(Added_dates_trim_date_added_count))

Added_dates_trim_date_added_good <- ddply(.data=Added_dates_trim_date_added_count,.variables=identifier,.fun=function(x){return(tail(x,1))},.progress="none")

rm(Added_dates_trim_date_added,Added_dates_trim_date_added_count)
rm(Added_dates,Added_dates_trim)
rm(Added_dates_cols_id,Added_dates_pulls_u,Added_dates_pulls_u_keep,Added_dates_pulls_u_drop)


###############################################################################
cat("SECTION: MEERGE DATE ADDED AND DEAD DATE","\n")
###############################################################################

Dates_comb0 <- data.frame(temp_col=sort(unique(c(Added_dates_trim_date_added_good[,identifier],Dead_dates_trim_dead_date_good[,identifier]))),stringsAsFactors=F)
colnames(Dates_comb0)[match("temp_col",names(Dates_comb0))] <- identifier

Dates_comb1 <- merge(Dates_comb0,Added_dates_trim_date_added_good[,colnames(Added_dates_trim_date_added_good)[!colnames(Added_dates_trim_date_added_good) %in% "freq"]],
                     by.x=identifier,by.y=identifier,all.x=T,all.y=F,sort=F,suffixes=c(".x",".y"))

Dates_comb2 <- merge(Dates_comb1,Dead_dates_trim_dead_date_good[,colnames(Dead_dates_trim_dead_date_good)[!colnames(Dead_dates_trim_dead_date_good) %in% "freq"]],
                     by.x=identifier,by.y=identifier,all.x=T,all.y=F,sort=F,suffixes=c(".x",".y"))

Dates_comb_full <- Dates_comb2

Dates_comb_full <- Dates_comb_full[order(Dates_comb_full[,identifier],Dates_comb_full[,"Date_Added"],Dates_comb_full[,"Dead_Date"]),] 
row.names(Dates_comb_full) <- seq(nrow(Dates_comb_full))

rm(Dates_comb0,Dates_comb1,Dates_comb2)
rm(Added_dates_trim_date_added_good,Dead_dates_trim_dead_date_good)


###############################################################################
cat("SECTION: IMPORT RET/AUM DATA","\n")
###############################################################################

NAV_AUM_cols_keep <- c("pull",identifier,"date",Ret_col_name,AUM_col_name)

NAV_AUM_input_files <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_NAV_AUM_files",".csv",sep=""),header=T,na.strings="NA",stringsAsFactors=F),stringsAsFactors=F)

NAV_AUM_input_files_trim <- NAV_AUM_input_files[!grepl("201207",NAV_AUM_input_files[,"file_name"]),]

NAV_AUM_concatentate0 <- alply(.data=NAV_AUM_input_files_trim,.margins=1,.fun=function(x,directory_in){
  
  # x <- NAV_AUM_input_files_trim[7,]
  # x <- NAV_AUM_input_files_trim[8,]
  # directory_in <- melt_folder_path
  
  #input <- data.frame(read.csv(file=paste(directory_in,"\\",x[,"pull"],".csv",sep=""),header=T,na.strings="NA",stringsAsFactors=F),stringsAsFactors=F)[NAV_AUM_cols_keep]
  input <- data.frame(read.columns(file=paste(directory_in,"\\",x[,"pull"],".csv",sep=""),required.col=NAV_AUM_cols_keep,sep=",",na.strings="NA",stringsAsFactors=F),
                      #Dead_Date=NA,Date_Added=NA,
                      #Ret_lag_temp=NA,AUM_lag_temp=NA,
                      #first_source=NA,last_source=NA,
                      #added_flag=NA,dead_flag=NA,first_flag=NA,revision_flag=NA,
                      stringsAsFactors=F)
  #colnames(NAV_AUM_concatentate)[match("Ret_lag_temp",names(NAV_AUM_concatentate))] <- paste(Ret_col_name,"lag1",sep="_")
  #colnames(NAV_AUM_concatentate)[match("AUM_lag_temp",names(NAV_AUM_concatentate))] <- paste(AUM_col_name,"lag1",sep="_")
  
  return(input)
  
},directory_in=melt_folder_path,.expand=T,.progress="text")

rm(NAV_AUM_input_files,NAV_AUM_input_files_trim)
invisible(gc(verbose=F,reset=T))

NAV_AUM_concatentate <- as.data.table(rbind.fill(NAV_AUM_concatentate0))
setkeyv(NAV_AUM_concatentate,NULL)

rm(NAV_AUM_concatentate0)
invisible(gc(verbose=F,reset=T))

# for(i in which(sapply(NAV_AUM_concatentate,class)=="character"))
# {
#   # NAV_AUM_concatentate[[i]] <- trim(NAV_AUM_concatentate[[i]])
#   NAV_AUM_concatentate[[i]] <- gsub("^\\s+|\\s+$","",NAV_AUM_concatentate[[i]],perl=T)
# }
# rm(i)
# for (i in 1:ncol(NAV_AUM_concatentate))
# {
#   NAV_AUM_concatentate[,i] <- unknownToNA(NAV_AUM_concatentate[,i],unknown=unknowns_strings,force=T)
#   NAV_AUM_concatentate[,i] <- ifelse(is.na(NAV_AUM_concatentate[,i]),NA,NAV_AUM_concatentate[,i])
# } 
# rm(i)

for (k in which(sapply(NAV_AUM_concatentate,class)=="character")) 
{
  set(NAV_AUM_concatentate,i=NULL,j=k,value=gsub("^\\s+|\\s+$","",NAV_AUM_concatentate[[k]],perl=T))
}
rm(k)
for (k in 1:ncol(NAV_AUM_concatentate)) 
{
  set(NAV_AUM_concatentate,i=NULL,j=k,value=unknownToNA(NAV_AUM_concatentate[[k]],unknown=unknowns_strings,force=T))
  set(NAV_AUM_concatentate,i=NULL,j=k,value=ifelse(is.na(NAV_AUM_concatentate[[k]]),NA,NAV_AUM_concatentate[[k]]))
}
rm(k)

#NAV_AUM_concatentate <- NAV_AUM_concatentate[,c("pull","Fund_ID","date",Ret_col_name,AUM_col_name)]

#NAV_AUM_concatentate[,"date"] <- as.yearmon(NAV_AUM_concatentate[,"date"],format="%b %Y")
#NAV_AUM_concatentate[,"date"] <- as.Date(NAV_AUM_concatentate[,"date"],format="%m-%d-%Y")

set(NAV_AUM_concatentate,i=NULL,j=which(colnames(NAV_AUM_concatentate)==c("date")),
    value=as.yearmon(NAV_AUM_concatentate[[which(colnames(NAV_AUM_concatentate)==c("date"))]],format="%b %Y"))
set(NAV_AUM_concatentate,i=NULL,j=which(colnames(NAV_AUM_concatentate)==c("date")),
    value=as.Date(NAV_AUM_concatentate[[which(colnames(NAV_AUM_concatentate)==c("date"))]],format="%m-%d-%Y"))


# Find lag pull and exclusion dates

NAV_AUM_pulls <- copy(NAV_AUM_concatentate)
NAV_AUM_pulls <- unique(NAV_AUM_pulls[,c(identifier,Ret_col_name,AUM_col_name):=list(NULL,NULL,NULL),by=NULL])

setnames(NAV_AUM_pulls,"pull","pull_trim2")
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c("pull_trim2")),
    value=gsub("[^[:digit:]]","",NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c("pull_trim2"))]]))
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c("pull_trim2")),
    value=as.integer(NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c("pull_trim2"))]]))
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c("pull_trim2")),
    value=paste(NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c("pull_trim2"))]],"01",sep=""))
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c("pull_trim2")),
    value=as.Date(NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c("pull_trim2"))]],format="%Y%m%d"))

setkeyv(NAV_AUM_pulls,c("date","pull_trim2"))
setorderv(NAV_AUM_pulls,c("date","pull_trim2"),c(1,1))
setkeyv(NAV_AUM_pulls,c("pull_trim2"))

NAV_AUM_pulls[,c("min_date","max_date"):=list(min(date,na.rm=T),max(date,na.rm=T)),by=c("pull_trim2")]
NAV_AUM_pulls <- unique(NAV_AUM_pulls[,c("date"):=list(NULL),by=NULL])

NAV_AUM_pulls[,c(paste("pull_trim2","lag1",sep="_")):=list(c(NA,head(pull_trim2,-1))),by=NULL]
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c(paste("pull_trim2","lag1",sep="_"))),
    value=as.Date(NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c(paste("pull_trim2","lag1",sep="_")))]]))

NAV_AUM_pulls[,c(paste("min_date","lag1",sep="_")):=list(c(NA,head(min_date,-1))),by=NULL]
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c(paste("min_date","lag1",sep="_"))),
    value=as.Date(NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c(paste("min_date","lag1",sep="_")))]]))

NAV_AUM_pulls[,c(paste("max_date","lag1",sep="_")):=list(c(NA,head(max_date,-1))),by=NULL]
set(NAV_AUM_pulls,i=NULL,j=which(colnames(NAV_AUM_pulls)==c(paste("max_date","lag1",sep="_"))),
    value=as.Date(NAV_AUM_pulls[[which(colnames(NAV_AUM_pulls)==c(paste("max_date","lag1",sep="_")))]]))

NAV_AUM_pulls[,c("min_date",paste("min_date","lag1",sep="_")):=list(NULL,NULL),by=NULL]

NAV_AUM_pulls <- adply(.data=as.data.frame(NAV_AUM_pulls,stringsAsFactors=F),.margins=1,.fun=function(x){
  
  # x <- as.data.frame(NAV_AUM_pulls,stringsAsFactors=F)[1,]
  # x <- as.data.frame(NAV_AUM_pulls,stringsAsFactors=F)[2,]
  
  out <- data.frame(x,temp_col=NA,stringsAsFactors=F)
  out[,"temp_col"] <- ifelse(is.na(out[,paste("max_date","lag1",sep="_")]),NA,seq(out[,paste("max_date","lag1",sep="_")],length=2,by="-12 months")[2])
  out[,"temp_col"] <- as.Date(out[,"temp_col"])
  return(out)
})

colnames(NAV_AUM_pulls)[match("temp_col",names(NAV_AUM_pulls))] <- paste("max_date","lag1","prior_12",sep="_")

NAV_AUM_pulls <- as.data.table(NAV_AUM_pulls)
setkeyv(NAV_AUM_pulls,c("pull_trim2"))
setorderv(NAV_AUM_pulls,c("pull_trim2"),c(1))

NAV_AUM_pulls_trim <- copy(NAV_AUM_pulls)
#NAV_AUM_pulls_trim[,c("max_date",paste("pull_trim2","lag1",sep="_"),paste("max_date","lag1",sep="_")):=list(NULL,NULL,NULL),by=NULL]
NAV_AUM_pulls_trim[,c("max_date",paste("max_date","lag1",sep="_")):=list(NULL,NULL),by=NULL]


## Find first and last sources for each fund

NAV_AUM_concatentate_trim_both <- copy(NAV_AUM_concatentate)[!(is.na(Monthly_Ret) & is.na(AUM))]
#NAV_AUM_concatentate_trim_both <- copy(NAV_AUM_concatentate)[(!is.na(Monthly_Ret) | !is.na(AUM))]

setnames(NAV_AUM_concatentate_trim_both,"pull","pull_trim2")
set(NAV_AUM_concatentate_trim_both,i=NULL,j=which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2")),
    value=gsub("[^[:digit:]]","",NAV_AUM_concatentate_trim_both[[which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2"))]]))
set(NAV_AUM_concatentate_trim_both,i=NULL,j=which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2")),
    value=as.integer(NAV_AUM_concatentate_trim_both[[which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2"))]]))
set(NAV_AUM_concatentate_trim_both,i=NULL,j=which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2")),
    value=paste(NAV_AUM_concatentate_trim_both[[which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2"))]],"01",sep=""))
set(NAV_AUM_concatentate_trim_both,i=NULL,j=which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2")),
    value=as.Date(NAV_AUM_concatentate_trim_both[[which(colnames(NAV_AUM_concatentate_trim_both)==c("pull_trim2"))]],format="%Y%m%d"))

setkeyv(NAV_AUM_concatentate_trim_both,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_trim_both,c(identifier,"date","pull_trim2"),c(1,1,1))
setkeyv(NAV_AUM_concatentate_trim_both,c(identifier,"date"))
#NAV_AUM_concatentate_trim_both[,c("first_source","last_source"):=list(min(pull_trim2,na.rm=T),max(pull_trim2,na.rm=T)),by=c(identifier,"date")]


#NAV_AUM_concatentate_trim_ret_drop_col <- colnames(NAV_AUM_concatentate_trim_both)[!(colnames(NAV_AUM_concatentate_trim_both) %in% c("AUM","first_source","last_source"))]
#NAV_AUM_concatentate_trim_ret <- as.data.table(as.data.frame(NAV_AUM_concatentate_trim_both,stringsAsFactors=F)[,NAV_AUM_concatentate_trim_ret_drop_col])
NAV_AUM_concatentate_trim_ret <- copy(NAV_AUM_concatentate_trim_both)
#NAV_AUM_concatentate_trim_ret[,c("AUM","first_source","last_source"):=list(NULL,NULL,NULL),by=NULL]
NAV_AUM_concatentate_trim_ret[,c("AUM"):=list(NULL),by=NULL]
setkeyv(NAV_AUM_concatentate_trim_ret,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_trim_ret,c(identifier,"date","pull_trim2"),c(1,1,1))
setkeyv(NAV_AUM_concatentate_trim_ret,c(identifier,"date"))
NAV_AUM_concatentate_trim_ret <- NAV_AUM_concatentate_trim_ret[!is.na(Monthly_Ret),]
NAV_AUM_concatentate_trim_ret[,c("first_source","last_source"):=list(min(pull_trim2,na.rm=T),max(pull_trim2,na.rm=T)),by=c(identifier,"date")]
NAV_AUM_concatentate_trim_ret[,c("pull_trim2","Monthly_Ret"):=list(NULL,NULL),by=NULL]
NAV_AUM_concatentate_trim_ret <- unique(NAV_AUM_concatentate_trim_ret)

#AV_AUM_concatentate_trim_aum_drop_col <- colnames(NAV_AUM_concatentate_trim_both)[!(colnames(NAV_AUM_concatentate_trim_both) %in% c("Monthly_Ret","first_source","last_source"))]
#NAV_AUM_concatentate_trim_aum <- as.data.table(as.data.frame(NAV_AUM_concatentate_trim_both,stringsAsFactors=F)[,NAV_AUM_concatentate_trim_aum_drop_col])
NAV_AUM_concatentate_trim_aum <- copy(NAV_AUM_concatentate_trim_both)
#NAV_AUM_concatentate_trim_aum[,c("Monthly_Ret","first_source","last_source"):=list(NULL,NULL,NULL),by=NULL]
NAV_AUM_concatentate_trim_aum[,c("Monthly_Ret"):=list(NULL),by=NULL]
setkeyv(NAV_AUM_concatentate_trim_aum,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_trim_aum,c(identifier,"date","pull_trim2"),c(1,1,1))
setkeyv(NAV_AUM_concatentate_trim_aum,c(identifier,"date"))
NAV_AUM_concatentate_trim_aum <- NAV_AUM_concatentate_trim_aum[!is.na(AUM),]
NAV_AUM_concatentate_trim_aum[,c("first_source","last_source"):=list(min(pull_trim2,na.rm=T),max(pull_trim2,na.rm=T)),by=c(identifier,"date")]
NAV_AUM_concatentate_trim_aum[,c("pull_trim2","AUM"):=list(NULL,NULL),by=NULL]
NAV_AUM_concatentate_trim_aum <- unique(NAV_AUM_concatentate_trim_aum)


# Expand Dates and Pulls

#NAV_AUM_concatentate_cast <- dcast(data=NAV_AUM_concatentate_trim[,c("pull",identifier,"date",Ret_col_name)],Fund_ID + date ~ pull,margins=NULL,subset=NULL,fill=NULL,drop=F,value.var="Monthly_Ret")
NAV_AUM_concatentate_cast <- dcast(data=as.data.frame(NAV_AUM_concatentate,stringsAsFactors=F)[,c("pull",identifier,"date",Ret_col_name)],Fund_ID + date ~ pull,margins=NULL,subset=NULL,fill=NULL,drop=F,value.var="Monthly_Ret")

rm(NAV_AUM_concatentate)
invisible(gc(verbose=F,reset=T))

NAV_AUM_concatentate_melt <- data.frame(melt(NAV_AUM_concatentate_cast,id=c(identifier,"date"),na.rm=F),AUM=NA,stringsAsFactors=F)

rm(NAV_AUM_concatentate_cast)
invisible(gc(verbose=F,reset=T))

colnames(NAV_AUM_concatentate_melt)[match("variable",names(NAV_AUM_concatentate_melt))] <- "pull"
colnames(NAV_AUM_concatentate_melt)[match("value",names(NAV_AUM_concatentate_melt))] <- Ret_col_name

NAV_AUM_concatentate_melt[sapply(NAV_AUM_concatentate_melt,is.factor)] <- lapply(NAV_AUM_concatentate_melt[sapply(NAV_AUM_concatentate_melt,is.factor)], as.character)
#NAV_AUM_concatentate_melt <- NAV_AUM_concatentate_melt[,colnames(NAV_AUM_concatentate_trim)]

NAV_AUM_concatentate_melt <- as.data.table(NAV_AUM_concatentate_melt[,colnames(NAV_AUM_concatentate_melt)[!colnames(NAV_AUM_concatentate_melt) %in% AUM_col_name]])

setnames(NAV_AUM_concatentate_melt,"pull","pull_trim2")
set(NAV_AUM_concatentate_melt,i=NULL,j=which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2")),
    value=gsub("[^[:digit:]]","",NAV_AUM_concatentate_melt[[which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2"))]]))
set(NAV_AUM_concatentate_melt,i=NULL,j=which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2")),
    value=as.integer(NAV_AUM_concatentate_melt[[which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2"))]]))
set(NAV_AUM_concatentate_melt,i=NULL,j=which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2")),
    value=paste(NAV_AUM_concatentate_melt[[which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2"))]],"01",sep=""))
set(NAV_AUM_concatentate_melt,i=NULL,j=which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2")),
    value=as.Date(NAV_AUM_concatentate_melt[[which(colnames(NAV_AUM_concatentate_melt)==c("pull_trim2"))]],format="%Y%m%d"))


## Merge

setcolorder(NAV_AUM_concatentate_melt,c(c(identifier,"date","pull_trim2"),
                                             colnames(NAV_AUM_concatentate_melt)[!(colnames(NAV_AUM_concatentate_melt) %in% c(identifier,"date","pull_trim2"))]))
NAV_AUM_concatentate_melt[,c(Ret_col_name):=list(NULL),by=NULL]
setkeyv(NAV_AUM_concatentate_melt,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_melt,c(identifier,"date","pull_trim2"),c(1,1,1))
invisible(gc(verbose=F,reset=T))

setcolorder(NAV_AUM_concatentate_trim_both,c(c(identifier,"date","pull_trim2"),
                                             colnames(NAV_AUM_concatentate_trim_both)[!(colnames(NAV_AUM_concatentate_trim_both) %in% c(identifier,"date","pull_trim2"))]))
setkeyv(NAV_AUM_concatentate_trim_both,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_trim_both,c(identifier,"date","pull_trim2"),c(1,1,1))
invisible(gc(verbose=F,reset=T))

NAV_AUM_concatentate_expand <- merge(NAV_AUM_concatentate_melt,NAV_AUM_concatentate_trim_both,
                                     by.x=c(identifier,"date","pull_trim2"),by.y=c(identifier,"date","pull_trim2"),
                                     all.x=T,all.y=F,sort=F,suffixes=c(".x",".y"))

rm(NAV_AUM_concatentate_melt,NAV_AUM_concatentate_trim_both)
invisible(gc(verbose=F,reset=T))

setkeyv(NAV_AUM_concatentate_expand,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_expand,c(identifier,"date","pull_trim2"),c(1,1,1))
setkeyv(NAV_AUM_concatentate_expand,c(identifier,"date"))

NAV_AUM_concatentate_expand[,c("ret_all_na","aum_all_na"):=list(mean(Monthly_Ret,na.rm=T),mean(AUM,na.rm=T)), by=c(identifier,"date")]

NAV_AUM_concatentate_expand_trim <- NAV_AUM_concatentate_expand[!(is.na(ret_all_na) & is.na(aum_all_na))]
NAV_AUM_concatentate_expand_trim[,c("ret_all_na","aum_all_na"):=list(NULL,NULL), by=NULL]

rm(NAV_AUM_concatentate_expand)
invisible(gc(verbose=F,reset=T))

## Merge in dead dates and date added

setkeyv(NAV_AUM_concatentate_expand_trim,c(identifier))
setorderv(NAV_AUM_concatentate_expand_trim,c(identifier),c(1))

Dates_comb_full <- as.data.table(Dates_comb_full)
setkeyv(Dates_comb_full,c(identifier))
setorderv(Dates_comb_full,c(identifier),c(1))

NAV_AUM_concatentate_expand2 <- merge(NAV_AUM_concatentate_expand_trim,Dates_comb_full,
                                     by.x=c(identifier),by.y=c(identifier),
                                     all.x=T,all.y=F,sort=F,suffixes=c(".x",".y"))

rm(NAV_AUM_concatentate_expand_trim,Dates_comb_full)
invisible(gc(verbose=F,reset=T))

setkeyv(NAV_AUM_concatentate_expand2,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_expand2,c(identifier,"date","pull_trim2"),c(1,1,1))


# Merge in lag pull and exclusion dates

setkeyv(NAV_AUM_concatentate_expand2,c("pull_trim2"))
setorderv(NAV_AUM_concatentate_expand2,c("pull_trim2"),c(1))

setkeyv(NAV_AUM_pulls_trim,c("pull_trim2"))
setorderv(NAV_AUM_pulls_trim,c("pull_trim2"),c(1))

NAV_AUM_concatentate_expand3 <- merge(NAV_AUM_concatentate_expand2,NAV_AUM_pulls_trim,
                                      by.x=c("pull_trim2"),by.y=c("pull_trim2"),
                                      all.x=T,all.y=F,sort=F,suffixes=c(".x",".y"))

rm(NAV_AUM_concatentate_expand2,NAV_AUM_pulls_trim)
invisible(gc(verbose=F,reset=T))

setkeyv(NAV_AUM_concatentate_expand3,c(identifier,"date","pull_trim2"))
setorderv(NAV_AUM_concatentate_expand3,c(identifier,"date","pull_trim2"),c(1,1,1))


###############################################################################
cat("SECTION: FIND REVISIONS","\n")
###############################################################################

Ret0 <- copy(NAV_AUM_concatentate_expand3)
AUM0 <- copy(NAV_AUM_concatentate_expand3)

rm(NAV_AUM_concatentate_expand3)
invisible(gc(verbose=F,reset=T))

## RET

Ret0[,c(AUM_col_name):=list(NULL),by=NULL]

setkeyv(Ret0,c(identifier,"date","pull_trim2"))
setorderv(Ret0,c(identifier,"date","pull_trim2"),c(1,1,1))
setkeyv(Ret0,c(identifier,"date"))

setkeyv(NAV_AUM_concatentate_trim_ret,c(identifier,"date"))
setorderv(NAV_AUM_concatentate_trim_ret,c(identifier,"date"),c(1,1))

Ret0_sources <- merge(Ret0,NAV_AUM_concatentate_trim_ret,
                      by.x=c(identifier,"date"),by.y=c(identifier,"date"),
                      all.x=T,all.y=F,sort=F,suffixes=c(".x",".y"))

Ret0_sources[,c(paste(Ret_col_name,"lag1",sep="_")):=list(c(NA,head(Monthly_Ret,-1))),by=c(identifier,"date")]
#Ret0_sources[,c(paste(Ret_col_name,"lead1",sep="_")):=list(c(tail(Monthly_Ret,-1),NA)),by=c(identifier,"date")]

setkeyv(Ret0_sources,c(identifier,"date","pull_trim2"))
setorderv(Ret0_sources,c(identifier,"date","pull_trim2"),c(1,1,1))

#Ret0_dead_good <- copy(Ret0_sources)
#Ret0_dead_good <- Ret0_dead_good[(is.na(Dead_Date) | date<=Dead_Date)]

#Ret0_dead_bad <- copy(Ret0_sources)
#Ret0_dead_bad <- Ret0_dead_bad[(!is.na(Dead_Date) & date>Dead_Date)]

Ret0_sources_trim0 <- copy(Ret0_sources)
Ret0_sources_trim1 <- Ret0_sources_trim0[(is.na(Dead_Date) | date<=Dead_Date)]
Ret0_sources_trim2 <- Ret0_sources_trim1[(pull_trim2>=first_source)]
Ret0_sources_trim3 <- Ret0_sources_trim2[(pull_trim2<=last_source)]
#Ret0_sources_trim4 <- Ret0_sources_trim3[!(max_date_lag1_prior_12 <= date < pull_trim2)]
Ret0_sources_trim4 <- Ret0_sources_trim3[date < max_date_lag1_prior_12]
Ret0_sources_trim5 <- Ret0_sources_trim4[(pull_trim2>first_source)]

Ret0_sources_flags <- copy(Ret0_sources_trim5)

rm(Ret0_sources_trim0,Ret0_sources_trim1,Ret0_sources_trim2,Ret0_sources_trim3,Ret0_sources_trim4,Ret0_sources_trim5)
invisible(gc(verbose=F,reset=T))

Ret0_sources_flags[,c("Date_Added","Dead_Date",paste("pull_trim2","lag1",sep="_"),"max_date_lag1_prior_12","first_source","last_source"):=list(NULL,NULL,NULL,NULL,NULL,NULL),by=NULL]
Ret0_sources_flags[,c("Status","Addition_Flag","Deletion_Flag","Revision_Flag"):=list(NA_character_,NA_integer_,NA_integer_,NA_integer_),by=NULL]
Ret0_sources_flags[,c("Revision_Size","Revision_Size_Abs","Revision_1BP","Revision_10BP","Revision_50BP","Revision_100BP"):=list(NA_real_,NA_real_,NA_integer_,NA_integer_,NA_integer_,NA_integer_),by=NULL]

val_idx <- which(colnames(Ret0_sources_flags)==c(Ret_col_name))
val_lag_idx <- which(colnames(Ret0_sources_flags)==c(paste(Ret_col_name,"lag1",sep="_")))
sflag_idx <- which(colnames(Ret0_sources_flags)==c("Status"))
aflag_idx <- which(colnames(Ret0_sources_flags)==c("Addition_Flag"))
dflag_idx <- which(colnames(Ret0_sources_flags)==c("Deletion_Flag"))
rflag_idx <- which(colnames(Ret0_sources_flags)==c("Revision_Flag"))
rs_idx <- which(colnames(Ret0_sources_flags)==c("Revision_Size"))
rsa_idx <- which(colnames(Ret0_sources_flags)==c("Revision_Size_Abs"))
r1bp_idx <- which(colnames(Ret0_sources_flags)==c("Revision_1BP"))
r10bp_idx <- which(colnames(Ret0_sources_flags)==c("Revision_10BP"))
r50bp_idx <- which(colnames(Ret0_sources_flags)==c("Revision_50BP"))
r100bp_idx <- which(colnames(Ret0_sources_flags)==c("Revision_100BP"))

set(Ret0_sources_flags,i=NULL,j=sflag_idx,value=ifelse((is.na(Ret0_sources_flags[[val_idx]]) & is.na(Ret0_sources_flags[[val_lag_idx]])),NA,
                                                       ifelse((!is.na(Ret0_sources_flags[[val_idx]]) & is.na(Ret0_sources_flags[[val_lag_idx]])),"Addition",
                                                              ifelse((is.na(Ret0_sources_flags[[val_idx]]) & !is.na(Ret0_sources_flags[[val_lag_idx]])),"Deletion",
                                                                     ifelse((Ret0_sources_flags[[val_idx]]!=Ret0_sources_flags[[val_lag_idx]]),"Revision","Nonrevision")))))

set(Ret0_sources_flags,i=NULL,j=aflag_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,ifelse(Ret0_sources_flags[[sflag_idx]]=="Addition",1,0)))
set(Ret0_sources_flags,i=NULL,j=dflag_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,ifelse(Ret0_sources_flags[[sflag_idx]]=="Deletion",1,0)))
set(Ret0_sources_flags,i=NULL,j=rflag_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,ifelse(Ret0_sources_flags[[sflag_idx]]=="Revision",1,0)))

set(Ret0_sources_flags,i=NULL,j=rs_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,
                                                    ifelse(Ret0_sources_flags[[rflag_idx]]==0,NA,Ret0_sources_flags[[val_idx]]-Ret0_sources_flags[[val_lag_idx]])))
set(Ret0_sources_flags,i=NULL,j=rsa_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,
                                                    ifelse(Ret0_sources_flags[[rflag_idx]]==0,NA,abs(Ret0_sources_flags[[rs_idx]]))))

## Returns are in percentage form
set(Ret0_sources_flags,i=NULL,j=r1bp_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,
                                                      ifelse(Ret0_sources_flags[[rflag_idx]]==0,NA,
                                                             ifelse(Ret0_sources_flags[[rsa_idx]]>=0.01,1,0))))
set(Ret0_sources_flags,i=NULL,j=r10bp_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,
                                                      ifelse(Ret0_sources_flags[[rflag_idx]]==0,NA,
                                                             ifelse(Ret0_sources_flags[[rsa_idx]]>=0.10,1,0))))
set(Ret0_sources_flags,i=NULL,j=r50bp_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,
                                                      ifelse(Ret0_sources_flags[[rflag_idx]]==0,NA,
                                                             ifelse(Ret0_sources_flags[[rsa_idx]]>=0.50,1,0))))
set(Ret0_sources_flags,i=NULL,j=r100bp_idx,value=ifelse(is.na(Ret0_sources_flags[[sflag_idx]]),NA,
                                                      ifelse(Ret0_sources_flags[[rflag_idx]]==0,NA,
                                                             ifelse(Ret0_sources_flags[[rsa_idx]]>=1.0,1,0))))


rm(val_idx,val_lag_idx,sflag_idx,aflag_idx,dflag_idx,rflag_idx)
rm(rs_idx,rsa_idx,r1bp_idx,r10bp_idx,r50bp_idx,r100bp_idx)
invisible(gc(verbose=F,reset=T))

Ret0_sources_flags_sum <- copy(Ret0_sources_flags)
Ret0_sources_flags_sum[,c(Ret_col_name,paste(Ret_col_name,"lag1",sep="_"),"Status"):=list(NULL,NULL,NULL),by=NULL]

# setkeyv(Ret0_sources_flags_sum,c(identifier,"date"))
# setorderv(Ret0_sources_flags_sum,c(identifier,"date"),c(1,1))
# Ret0_sources_flags_sum[,c("Addition_Flag_total","Deletion_Flag_total","Revision_Flag_total"):=list(sum(Addition_Flag,na.rm=T),sum(Deletion_Flag,na.rm=T),sum(Revision_Flag,na.rm=T)),by=c(identifier,"date")]
# Ret0_sources_flags_sum[,c("pull_trim2","Addition_Flag","Deletion_Flag","Revision_Flag"):=list(NULL,NULL,NULL,NULL),by=NULL]

# setkeyv(Ret0_sources_flags_sum,c(identifier,"pull_trim2"))
# setorderv(Ret0_sources_flags_sum,c(identifier,"pull_trim2"),c(1,1))
# Ret0_sources_flags_sum[,c("Addition_Flag_total","Deletion_Flag_total","Revision_Flag_total"):=list(sum(Addition_Flag,na.rm=T),sum(Deletion_Flag,na.rm=T),sum(Revision_Flag,na.rm=T)),by=c(identifier,"pull_trim2")]
# Ret0_sources_flags_sum[,c("date","Addition_Flag","Deletion_Flag","Revision_Flag"):=list(NULL,NULL,NULL,NULL),by=NULL]

setkeyv(Ret0_sources_flags_sum,c(identifier))
setorderv(Ret0_sources_flags_sum,c(identifier),c(1))
Ret0_sources_flags_sum[,c("Addition_Flag_total","Deletion_Flag_total","Revision_Flag_total"):=list(sum(Addition_Flag,na.rm=T),sum(Deletion_Flag,na.rm=T),sum(Revision_Flag,na.rm=T)),by=c(identifier)]
Ret0_sources_flags_sum[,c("Revision_1BP_Total","Revision_10BP_Total","Revision_50BP_Total","Revision_100BP_Total"):=list(sum(Revision_1BP,na.rm=T),sum(Revision_10BP,na.rm=T),sum(Revision_50BP,na.rm=T),sum(Revision_100BP,na.rm=T)),by=c(identifier)]

Ret0_sources_flags_sum[,c("date","pull_trim2","Addition_Flag","Deletion_Flag","Revision_Flag"):=list(NULL,NULL,NULL,NULL,NULL),by=NULL]
Ret0_sources_flags_sum[,c("Revision_Size","Revision_Size_Abs"):=list(NULL,NULL),by=NULL]
Ret0_sources_flags_sum[,c("Revision_1BP","Revision_10BP","Revision_50BP","Revision_100BP"):=list(NULL,NULL,NULL,NULL),by=NULL]

Ret0_sources_flags_sum <- unique(Ret0_sources_flags_sum)
# Ret0_sources_flags_sum[1:50,]


Ret0_sources_flags_sum[,c("Addition_DV","Deletion_DV","Revision_DV"):=list(NA_integer_,NA_integer_,NA_integer_),by=NULL]
Ret0_sources_flags_sum[,c("Revision_1BP_DV","Revision_10BP_DV","Revision_50BP_DV","Revision_100BP_DV"):=list(NA_integer_,NA_integer_,NA_integer_,NA_integer_),by=NULL]

set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Addition_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Addition_Flag_total"))]]),NA,
                                                                                                           ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Addition_Flag_total"))]]>0,1,0)))
set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Deletion_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Deletion_Flag_total"))]]),NA,
                                                                                                           ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Deletion_Flag_total"))]]>0,1,0)))
set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Revision_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_Flag_total"))]]),NA,
                                                                                                           ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_Flag_total"))]]>0,1,0)))

set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Revision_1BP_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_1BP_Total"))]]),NA,
                                                                                                               ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_1BP_Total"))]]>0,1,0)))
set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Revision_10BP_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_10BP_Total"))]]),NA,
                                                                                                                ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_10BP_Total"))]]>0,1,0)))
set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Revision_50BP_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_50BP_Total"))]]),NA,
                                                                                                                ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_50BP_Total"))]]>0,1,0)))
set(Ret0_sources_flags_sum,i=NULL,j=which(colnames(Ret0_sources_flags_sum)==c("Revision_100BP_DV")),value=ifelse(is.na(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_100BP_Total"))]]),NA,
                                                                                                                 ifelse(Ret0_sources_flags_sum[[which(colnames(Ret0_sources_flags_sum)==c("Revision_100BP_Total"))]]>0,1,0)))


Ret0_sources_flags_sum <- as.data.frame(Ret0_sources_flags_sum,stringsAsFactors=F)

# mean(Ret0_sources_flags_sum[,"Addition_DV"],na.rm=T)
# mean(Ret0_sources_flags_sum[,"Deletion_DV"],na.rm=T)
# mean(Ret0_sources_flags_sum[,"Revision_DV"],na.rm=T)
# mean(Ret0_sources_flags_sum[,"Revision_1BP_DV"],na.rm=T)
# mean(Ret0_sources_flags_sum[,"Revision_10BP_DV"],na.rm=T)
# mean(Ret0_sources_flags_sum[,"Revision_50BP_DV"],na.rm=T)
# mean(Ret0_sources_flags_sum[,"Revision_100BP_DV"],na.rm=T)


write.csv(Ret0_sources,file=paste(revision_folder_path,"//","Ret0_sources",".csv",sep=""),row.names=FALSE)
write.csv(Ret0_sources_flags,file=paste(revision_folder_path,"//","Ret0_sources_flags",".csv",sep=""),row.names=FALSE)
write.csv(Ret0_sources_flags_sum,file=paste(revision_folder_path,"//","Ret0_sources_flags_sum",".csv",sep=""),row.names=FALSE)




rm(Ret0,NAV_AUM_concatentate_trim_ret)
rm(AUM0,NAV_AUM_concatentate_trim_aum)
rm(NAV_AUM_pulls)
invisible(gc(verbose=F,reset=T))






