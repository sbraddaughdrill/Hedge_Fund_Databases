# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Stats_Expand_Returns.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Expand yearly and monthly returns from Eurekahedge Data
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

#Check to see if nav_aum_ret folder exists.  If not, create it.
nav_aum_ret_folder_path <- paste(output_directory, "NAV_AUM_Ret", sep = "//", collapse = "//")  
create_directory(nav_aum_ret_folder_path,remove=1)

Stats <- data.frame(read.csv(file=paste(nav_aum_ret_folder_path,"\\","EurekahedgeHF_Stats_with_Rets",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                    stringsAsFactors=FALSE)

for(i in which(sapply(Stats,class)=="character"))
{
  Stats[[i]] = trim(Stats[[i]])
}
rm2(i)
for (i in 1:ncol(Stats))
{
  Stats[,i] <- unknownToNA(Stats[,i], unknown=unknowns_strings,force=TRUE)
  Stats[,i] <- ifelse(is.na(Stats[,i]),NA, Stats[,i])
} 
rm2(i)

Stats  <- Stats[order(Stats[,"pull"],Stats[,"Fund_ID"],Stats[,"Fund_Name"]),]
row.names(Stats) <- seq(nrow(Stats))


###############################################################################
cat("SECTION: SEPERATE RETURN AND NONRETURN COLUMNS", "\n")
###############################################################################

Stats_all_cols <- colnames(Stats)

Stats_id_cols <- c("pull","Fund_ID","Fund_Name","Date_Added","Flagship","Closed","Limited","Dead","Dead_Date","Dead_Reason")
Stats_nonid_cols <- Stats_all_cols[!(Stats_all_cols %in% c(Stats_id_cols))]

Stats_yearlyret_cols <- Stats_nonid_cols[grep("yearlyret", Stats_nonid_cols)]
Stats_nonyearlyret_cols <- Stats_nonid_cols[!(Stats_nonid_cols %in% c(Stats_yearlyret_cols))]

Stats_monthlyret_cols <- Stats_nonyearlyret_cols[grep("monthlyret", Stats_nonyearlyret_cols)]
Stats_other_cols <- Stats_nonyearlyret_cols[!(Stats_nonyearlyret_cols %in% c(Stats_monthlyret_cols))]

rm2(Stats_all_cols,Stats_nonid_cols,Stats_nonyearlyret_cols)


#Seperate nonreturn columns

Stats_nonreturn <- Stats[,c(Stats_id_cols,Stats_other_cols)]

Stats_nonreturn_inception_col <- Stats_other_cols[grep("Inception", Stats_other_cols)]
Stats_nonreturn_noninception_col <- Stats_other_cols[!(Stats_other_cols %in% c(Stats_nonreturn_inception_col))]

Stats_nonreturn_ret_col <- Stats_nonreturn_noninception_col[grep("Return", Stats_nonreturn_noninception_col)]
Stats_nonreturn_nonret_col <- Stats_nonreturn_noninception_col[!(Stats_nonreturn_noninception_col %in% c(Stats_nonreturn_ret_col))]

Stats_nonreturn_ratio_col <- Stats_nonreturn_nonret_col[grep("Ratio", Stats_nonreturn_nonret_col)]
Stats_nonreturn_nonratio_col <- Stats_nonreturn_nonret_col[!(Stats_nonreturn_nonret_col %in% c(Stats_nonreturn_ratio_col))]

Stats_nonreturn <- Stats_nonreturn[,c(Stats_id_cols,Stats_nonreturn_inception_col,Stats_nonreturn_ret_col,
                                      Stats_nonreturn_nonratio_col,Stats_nonreturn_ratio_col)]

rm2(Stats_nonreturn_inception_col,Stats_nonreturn_noninception_col)
rm2(Stats_nonreturn_ret_col,Stats_nonreturn_nonret_col)
rm2(Stats_nonreturn_ratio_col,Stats_nonreturn_nonratio_col)

Stats_nonreturn  <- Stats_nonreturn[order(Stats_nonreturn[,"Fund_ID"],Stats_nonreturn[,"pull"]),]
row.names(Stats_nonreturn) <- seq(nrow(Stats_nonreturn))


#Seperate yearly return columns

Stats_yearly_return <- Stats[,c("pull","Fund_ID",Stats_yearlyret_cols)]

Stats_yearly_return  <- Stats_yearly_return[order(Stats_yearly_return[,"Fund_ID"],Stats_yearly_return[,"pull"]),]
row.names(Stats_yearly_return) <- seq(nrow(Stats_yearly_return))


#Seperate monthly return columns

Stats_monthly_return <- Stats[,c("pull","Fund_ID",Stats_monthlyret_cols)]

Stats_monthly_return  <- Stats_monthly_return[order(Stats_monthly_return[,"Fund_ID"],Stats_monthly_return[,"pull"]),]
row.names(Stats_monthly_return) <- seq(nrow(Stats_monthly_return))

rm2(Stats)


###############################################################################
cat("SECTION: SEPERATE MONTHLY RETURNS", "\n")
###############################################################################

#aa <- Stats_monthly_return[Stats_monthly_return[,"Fund_ID"] %in% c(5000,5002,5003,5004),]

#Melt Returns

Stats_monthly_return_melt <- data.frame(melt(Stats_monthly_return, id=c("pull","Fund_ID"),na.rm=FALSE,value.name = "value"),
                                        yr=NA, month=NA, stringsAsFactors=FALSE)

rm2(Stats_monthly_return)

for(j in 1:ncol(Stats_monthly_return_melt))
{
  #j <- 1
  Stats_monthly_return_melt[,j] = trim(Stats_monthly_return_melt[,j])
}
rm2(j)

colnames(Stats_monthly_return_melt)[match("variable",names(Stats_monthly_return_melt))] <- "date"
colnames(Stats_monthly_return_melt)[match("value",names(Stats_monthly_return_melt))] <- "Monthly_Ret2"

Stats_monthly_return_melt <- Stats_monthly_return_melt[order(Stats_monthly_return_melt[,"Fund_ID"],
                                                             Stats_monthly_return_melt[,"pull"],
                                                             Stats_monthly_return_melt[,"date"]),]
row.names(Stats_monthly_return_melt) <- seq(nrow(Stats_monthly_return_melt))


Stats_monthly_return_melt[,"date"] <- gsub(pattern="Return_", replacement="", x=Stats_monthly_return_melt[,"date"])
Stats_monthly_return_melt[,"date"] <- gsub(pattern="\\_", replacement="-", x=Stats_monthly_return_melt[,"date"])
Stats_monthly_return_melt[,"date"] <- gsub(pattern=" ", replacement="", x=Stats_monthly_return_melt[,"date"])
Stats_monthly_return_melt[,"date"] <- gsub(pattern=" ", replacement="", x=Stats_monthly_return_melt[,"date"])
Stats_monthly_return_melt[,"date"] <- gsub(pattern="_monthlyret", replacement="", x=Stats_monthly_return_melt[,"date"])
Stats_monthly_return_melt[,"date"] <- gsub(pattern="-monthlyret", replacement="", x=Stats_monthly_return_melt[,"date"])

for(j in 1:ncol(Stats_monthly_return_melt))
{
  #j <- 1
  Stats_monthly_return_melt[,j] = trim(Stats_monthly_return_melt[,j])
  
}
rm2(j)

#Find Year and Months

Stats_monthly_return_melt[,"yr"] <- Stats_monthly_return_melt[,"date"]
Stats_monthly_return_melt[,"yr"] <- substr(Stats_monthly_return_melt[,"yr"], 1, 4)
Stats_monthly_return_melt[,"yr"] <- as.integer(Stats_monthly_return_melt[,"yr"])

Stats_monthly_return_melt[,"month"] <- Stats_monthly_return_melt[,"date"]
Stats_monthly_return_melt[,"month"] <- substr(Stats_monthly_return_melt[,"month"], 6, 7)
Stats_monthly_return_melt[,"month"] <- as.integer(Stats_monthly_return_melt[,"month"])

#Stats_monthly_return_melt[,"date"] <- format(as.Date(paste(Stats_monthly_return_melt[,"yr"],Stats_monthly_return_melt[,"month"],"01",sep="-"),
#                                                    format="%y-%m-%d"),"%Y-%m-%d")

Stats_monthly_return_melt[,"date"] <- paste(Stats_monthly_return_melt[,"yr"],Stats_monthly_return_melt[,"month"],"01",sep="-")
Stats_monthly_return_melt[,"date"] <- as.Date(Stats_monthly_return_melt[,"date"],format="%Y-%m-%d")

Stats_monthly_return_melt[,"yr"] <- year(Stats_monthly_return_melt[,"date"])
Stats_monthly_return_melt[,"month"] <- month(Stats_monthly_return_melt[,"date"])

Stats_monthly_return_melt[,"Monthly_Ret2"] <- as.numeric(Stats_monthly_return_melt[,"Monthly_Ret2"])

Stats_monthly_return_melt <- Stats_monthly_return_melt[rowSums(is.na(Stats_monthly_return_melt[,1:ncol(Stats_monthly_return_melt)]))<ncol(Stats_monthly_return_melt),]

Stats_monthly_return_melt[,"Fund_ID"] <- as.integer(Stats_monthly_return_melt[,"Fund_ID"])

Stats_monthly_return_melt <- Stats_monthly_return_melt[order(Stats_monthly_return_melt[,"Fund_ID"],
                                                             Stats_monthly_return_melt[,"yr"],
                                                             Stats_monthly_return_melt[,"month"],
                                                             Stats_monthly_return_melt[,"pull"]),]
row.names(Stats_monthly_return_melt) <- seq(nrow(Stats_monthly_return_melt))

#Only keep one return per date

Stats_monthly_return_melt_u0 <- Stats_monthly_return_melt[!is.na(Stats_monthly_return_melt[,"Monthly_Ret2"]),]

rm2(Stats_monthly_return_melt)

Stats_monthly_return_melt_u0 <- Stats_monthly_return_melt_u0[order(Stats_monthly_return_melt_u0[,"Fund_ID"],
                                                                   Stats_monthly_return_melt_u0[,"yr"],
                                                                   Stats_monthly_return_melt_u0[,"month"]),]
row.names(Stats_monthly_return_melt_u0) <- seq(nrow(Stats_monthly_return_melt_u0))

Stats_monthly_return_melt_u <- ddply(.data=Stats_monthly_return_melt_u0, .variables=c("Fund_ID","date"),function(x){ 
  return(tail(x,1))}, .progress = "text")

#Stats_monthly_return_melt_u0 <- data.table(Stats_monthly_return_melt_u0)
#setkeyv(Stats_monthly_return_melt_u0, c("Fund_ID","date"))
#Stats_monthly_return_melt_u <- Stats_monthly_return_melt_u0[, tail(.SD, 1), by=c("Fund_ID","date")]

rm2(Stats_monthly_return_melt_u0)

Stats_monthly_return_melt_u <- Stats_monthly_return_melt_u[order(Stats_monthly_return_melt_u[,"Fund_ID"],
                                                                 Stats_monthly_return_melt_u[,"yr"],
                                                                 Stats_monthly_return_melt_u[,"month"]),]
row.names(Stats_monthly_return_melt_u) <- seq(nrow(Stats_monthly_return_melt_u))

Stats_monthly_return_melt_u <- Stats_monthly_return_melt_u[,c(colnames(Stats_monthly_return_melt_u)[!(colnames(Stats_monthly_return_melt_u) %in% "Monthly_Ret2")],
                                                              "Monthly_Ret2")]


###############################################################################
cat("SECTION: SEPERATE YEARLY RETURNS", "\n")
###############################################################################

#aa <- Stats_yearly_return[Stats_yearly_return[,"Fund_ID"] %in% c(5000,5002,5003,5004),]

#Melt Returns

Stats_yearly_return_melt <- data.frame(melt(Stats_yearly_return, id=c("pull","Fund_ID"),na.rm=FALSE,value.name = "value"),
                                       yr=NA, stringsAsFactors=FALSE)

rm2(Stats_yearly_return)

for(j in 1:ncol(Stats_yearly_return_melt))
{
  #j <- 1
  Stats_yearly_return_melt[,j] = trim(Stats_yearly_return_melt[,j])
}
rm2(j)

colnames(Stats_yearly_return_melt)[match("variable",names(Stats_yearly_return_melt))] <- "date"
colnames(Stats_yearly_return_melt)[match("value",names(Stats_yearly_return_melt))] <- "Yearly_Ret2"

Stats_yearly_return_melt <- Stats_yearly_return_melt[order(Stats_yearly_return_melt[,"Fund_ID"],
                                                           Stats_yearly_return_melt[,"pull"],
                                                           Stats_yearly_return_melt[,"date"]),]
row.names(Stats_yearly_return_melt) <- seq(nrow(Stats_yearly_return_melt))


Stats_yearly_return_melt[,"date"] <- gsub(pattern="Return_", replacement="", x=Stats_yearly_return_melt[,"date"])
Stats_yearly_return_melt[,"date"] <- gsub(pattern="\\_", replacement="-", x=Stats_yearly_return_melt[,"date"])
Stats_yearly_return_melt[,"date"] <- gsub(pattern=" ", replacement="", x=Stats_yearly_return_melt[,"date"])
Stats_yearly_return_melt[,"date"] <- gsub(pattern=" ", replacement="", x=Stats_yearly_return_melt[,"date"])
Stats_yearly_return_melt[,"date"] <- gsub(pattern="_yearlyret", replacement="", x=Stats_yearly_return_melt[,"date"])
Stats_yearly_return_melt[,"date"] <- gsub(pattern="-yearlyret", replacement="", x=Stats_yearly_return_melt[,"date"])

for(j in 1:ncol(Stats_yearly_return_melt))
{
  #j <- 1
  Stats_yearly_return_melt[,j] = trim(Stats_yearly_return_melt[,j])
  
}
rm2(j)

#Find Year and Months

Stats_yearly_return_melt[,"yr"] <- Stats_yearly_return_melt[,"date"]
Stats_yearly_return_melt[,"yr"] <- substr(Stats_yearly_return_melt[,"yr"], 1, 4)
Stats_yearly_return_melt[,"yr"] <- as.integer(Stats_yearly_return_melt[,"yr"])

Stats_yearly_return_melt[,"Yearly_Ret2"] <- as.numeric(Stats_yearly_return_melt[,"Yearly_Ret2"])

Stats_yearly_return_melt <- Stats_yearly_return_melt[rowSums(is.na(Stats_yearly_return_melt[,1:ncol(Stats_yearly_return_melt)]))<ncol(Stats_yearly_return_melt),]

Stats_yearly_return_melt[,"Fund_ID"] <- as.integer(Stats_yearly_return_melt[,"Fund_ID"])

Stats_yearly_return_melt <- Stats_yearly_return_melt[order(Stats_yearly_return_melt[,"Fund_ID"],
                                                           Stats_yearly_return_melt[,"yr"],
                                                           Stats_yearly_return_melt[,"pull"]),]
row.names(Stats_yearly_return_melt) <- seq(nrow(Stats_yearly_return_melt))

#Only keep one return per date

Stats_yearly_return_melt_u0 <- Stats_yearly_return_melt[!is.na(Stats_yearly_return_melt[,"Yearly_Ret2"]),]

rm2(Stats_yearly_return_melt)

Stats_yearly_return_melt_u0 <- Stats_yearly_return_melt_u0[order(Stats_yearly_return_melt_u0[,"Fund_ID"],
                                                                 Stats_yearly_return_melt_u0[,"yr"]),]
row.names(Stats_yearly_return_melt_u0) <- seq(nrow(Stats_yearly_return_melt_u0))

Stats_yearly_return_melt_u <- ddply(.data=Stats_yearly_return_melt_u0, .variables=c("Fund_ID","date"),function(x){ 
  return(tail(x,1))}, .progress = "text")

#Stats_yearly_return_melt_u0 <- data.table(Stats_yearly_return_melt_u0)
#setkeyv(Stats_yearly_return_melt_u0, c("Fund_ID","date"))
#Stats_yearly_return_melt_u <- Stats_yearly_return_melt_u0[, tail(.SD, 1), by=c("Fund_ID","date")]

rm2(Stats_yearly_return_melt_u0)

Stats_yearly_return_melt_u <- Stats_yearly_return_melt_u[order(Stats_yearly_return_melt_u[,"Fund_ID"],
                                                               Stats_yearly_return_melt_u[,"yr"]),]
row.names(Stats_yearly_return_melt_u) <- seq(nrow(Stats_yearly_return_melt_u))

Stats_yearly_return_melt_u <- Stats_yearly_return_melt_u[,c(colnames(Stats_yearly_return_melt_u)[!(colnames(Stats_yearly_return_melt_u) %in% "Yearly_Ret2")],
                                                            "Yearly_Ret2")]


###############################################################################
cat("SECTION: IMPORT NAV & AUM", "\n")
###############################################################################

NAV_AUM <- data.frame(read.csv(file=paste(nav_aum_ret_folder_path,"\\","EurekahedgeHF_NAV_AUM_merge",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      stringsAsFactors=FALSE)

#NAV_AUM[,"pull"] <- as.Date(paste("01",NAV_AUM[,"pull"],sep="-"),format="%d-%m-%Y")
#NAV_AUM[,"pull"] <- as.POSIXct(paste("01",NAV_AUM[,"pull"],sep="-"),format="%d-%m-%Y")


###############################################################################
cat("SECTION: MERGE MONTHLY_RET2 AND YEARLY_RET2", "\n")
###############################################################################

Stats_monthly_yearly_return_yr_min <- min(c(Stats_yearly_return_melt_u[,"yr"],Stats_monthly_return_melt_u[,"yr"]))
Stats_monthly_yearly_return_yr_max <- max(c(Stats_yearly_return_melt_u[,"yr"],Stats_monthly_return_melt_u[,"yr"]))
Stats_monthly_yearly_return_yr_seq <- seq(Stats_monthly_yearly_return_yr_min,Stats_monthly_yearly_return_yr_max,by=1)

Stats_monthly_yearly_return_id <- unique(c(Stats_yearly_return_melt_u[,"Fund_ID"],Stats_monthly_return_melt_u[,"Fund_ID"]))

Stats_monthly_yearly_return0 <- expand.grid(Fund_ID=Stats_monthly_yearly_return_id, 
                                            yr=Stats_monthly_yearly_return_yr_seq,
                                            month=seq(1,12,by=1))
Stats_monthly_yearly_return0 <- as.data.frame(Stats_monthly_yearly_return0,stringsAsFactors=FALSE)
Stats_monthly_yearly_return0[,"yr"] <- as.integer(Stats_monthly_yearly_return0[,"yr"])
Stats_monthly_yearly_return0[,"month"] <- as.integer(Stats_monthly_yearly_return0[,"month"])

rm2(Stats_monthly_yearly_return_yr_min,Stats_monthly_yearly_return_yr_max,Stats_monthly_yearly_return_yr_seq)
rm2(Stats_monthly_yearly_return_id)

Stats_monthly_yearly_return0 <- Stats_monthly_yearly_return0[order(Stats_monthly_yearly_return0[,"Fund_ID"],
                                                                   Stats_monthly_yearly_return0[,"yr"],
                                                                   Stats_monthly_yearly_return0[,"month"]),]
row.names(Stats_monthly_yearly_return0) <- seq(nrow(Stats_monthly_yearly_return0))

Stats_monthly_yearly_return1 <- merge(Stats_monthly_yearly_return0,
                                      Stats_yearly_return_melt_u[,!names(Stats_yearly_return_melt_u) %in% c("pull","date")], 
                                      by.x=c("Fund_ID","yr"), by.y=c("Fund_ID","yr"), 
                                      all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(Stats_monthly_yearly_return0,Stats_yearly_return_melt_u)

Stats_monthly_yearly_return <- merge(Stats_monthly_yearly_return1,
                                     Stats_monthly_return_melt_u[,!names(Stats_monthly_return_melt_u) %in% c("pull","date")], 
                                     by.x=c("Fund_ID","yr","month"), by.y=c("Fund_ID","yr","month"), 
                                     all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(Stats_monthly_yearly_return1,Stats_monthly_return_melt_u)

Stats_monthly_yearly_return <- Stats_monthly_yearly_return[,c("Fund_ID","yr","month",
                                                              colnames(Stats_monthly_yearly_return)[!(colnames(Stats_monthly_yearly_return) %in% c("Fund_ID","yr","month"))])]

Stats_monthly_yearly_return <- Stats_monthly_yearly_return[order(Stats_monthly_yearly_return[,"Fund_ID"],
                                                                 Stats_monthly_yearly_return[,"yr"],
                                                                 Stats_monthly_yearly_return[,"month"]),]
row.names(Stats_monthly_yearly_return) <- seq(nrow(Stats_monthly_yearly_return))

NAV_AUM_rets <- merge(NAV_AUM,Stats_monthly_yearly_return, 
                      by.x=c("Fund_ID","yr","month"), by.y=c("Fund_ID","yr","month"), 
                      all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(NAV_AUM,Stats_monthly_yearly_return)

NAV_AUM_rets_id_cols <- c("pull_trim","pull","Fund_ID","date","yr","month","min_date","max_date","bad_min","bad_max")
NAV_AUM_rets_other_cols <- c("Monthly_Ret","Monthly_Ret2","Yearly_Ret2","AUM")

NAV_AUM_rets <- NAV_AUM_rets[,c(NAV_AUM_rets_id_cols,NAV_AUM_rets_other_cols)]

NAV_AUM_rets <- NAV_AUM_rets[order(NAV_AUM_rets[,"Fund_ID"],
                                   NAV_AUM_rets[,"yr"],
                                   NAV_AUM_rets[,"month"]),]
row.names(NAV_AUM_rets) <- seq(nrow(NAV_AUM_rets))

rm2(NAV_AUM_rets_id_cols,NAV_AUM_rets_other_cols)


###############################################################################
cat("SECTION: REMOVE RETURNS AFTER FUND DEATH", "\n")
###############################################################################

NAV_AUM_yr_min <- min(NAV_AUM_rets[,"yr"])
NAV_AUM_yr_max <- max(NAV_AUM_rets[,"yr"])
NAV_AUM_yr_seq <- seq(NAV_AUM_yr_min,NAV_AUM_yr_max,by=1)

dead_dates_flag <- data.frame(Stats_nonreturn[,c("pull","Fund_ID","Dead_Date")],dead_flag=NA,stringsAsFactors=FALSE)
dead_dates_flag[,"dead_flag"] <- ifelse(is.na(dead_dates_flag[,"Dead_Date"]),0,1)

dead_dates_flag[,"Dead_Date"] <- as.yearmon(dead_dates_flag[,"Dead_Date"],format="%b %Y")
#dead_dates_flag[,"Dead_Date"] <- as.POSIXct(dead_dates_flag[,"Dead_Date"],format="%Y-%m-%d")
dead_dates_flag[,"Dead_Date"] <- as.Date(dead_dates_flag[,"Dead_Date"],format="%Y-%m-%d")

#Find funds with multiple death dates

dead_dates_no <- dead_dates_flag[dead_dates_flag[,"dead_flag"]==0,]
dead_dates_yes <- dead_dates_flag[dead_dates_flag[,"dead_flag"]==1,]

rm2(dead_dates_flag)

dead_dates_u <- unique(dead_dates_yes[,!(colnames(dead_dates_yes) %in% c("pull"))])
dead_dates_u <- dead_dates_u[order(dead_dates_u[,"Fund_ID"],dead_dates_u[,"Dead_Date"]),]
row.names(dead_dates_u) <- seq(nrow(dead_dates_u))

dead_dates_u_count <- count(dead_dates_u[,"Fund_ID"])
colnames(dead_dates_u_count)[match("x",names(dead_dates_u_count))] <- "Fund_ID"

rm2(dead_dates_u)

dead_dates_u_count_single <- dead_dates_u_count[dead_dates_u_count[,"freq"]==1,]
dead_dates_u_single <- dead_dates_yes[(dead_dates_yes[,"Fund_ID"] %in% dead_dates_u_count_single[,"Fund_ID"]),]

dead_dates_u_count_multiple <- dead_dates_u_count[dead_dates_u_count[,"freq"]>=2,]
dead_dates_u_multiple <- dead_dates_yes[(dead_dates_yes[,"Fund_ID"] %in% dead_dates_u_count_multiple[,"Fund_ID"]),]

rm2(dead_dates_u_count,dead_dates_u_count_single,dead_dates_u_count_multiple)


#Take the latest death from funds with multiple death funds

colnames(dead_dates_u_multiple)[match("Dead_Date",names(dead_dates_u_multiple))] <- "Dead_Date_org"

dead_dates_u_multiple_replace0 <- ddply(.data=dead_dates_u_multiple, .variables="Fund_ID",function(x){
  
  return(max(x[,"Dead_Date_org"],na.rm=TRUE))}, .progress = "none")

colnames(dead_dates_u_multiple_replace0) <- c("Fund_ID","Dead_Date")

dead_dates_u_multiple_replace <- merge(dead_dates_u_multiple, 
                                       dead_dates_u_multiple_replace0, 
                                       by.x=c("Fund_ID"), by.y=c("Fund_ID"), 
                                       all.x=TRUE, all.y=FALSE, sort=FALSE,pull=c(".x",".y"))

rm2(dead_dates_u_multiple,dead_dates_u_multiple_replace0)

dead_dates_u_multiple_replace[,"Dead_Date"] <- as.character(dead_dates_u_multiple_replace[,"Dead_Date"])
dead_dates_u_multiple_replace[,"Dead_Date_org"] <- as.character(dead_dates_u_multiple_replace[,"Dead_Date_org"])

dead_dates_u_multiple_replace[,"Dead_Date"] <- ifelse(is.na(dead_dates_u_multiple_replace[,"Dead_Date_org"]),
                                                      NA,dead_dates_u_multiple_replace[,"Dead_Date"])

#dead_dates_u_multiple_replace[,"Dead_Date"] <- as.POSIXct(dead_dates_u_multiple_replace[,"Dead_Date"],format="%Y-%m-%d")
dead_dates_u_multiple_replace[,"Dead_Date"] <- as.Date(dead_dates_u_multiple_replace[,"Dead_Date"],format="%Y-%m-%d")

#Concatenate Death.Dates

dead_dates_fixed0 <- rbind(dead_dates_u_single,dead_dates_u_multiple_replace[,!(colnames(dead_dates_u_multiple_replace) %in% c("Dead_Date_org"))])

rm2(dead_dates_u_single,dead_dates_u_multiple_replace)

dead_dates_fixed0 <- dead_dates_fixed0[order(dead_dates_fixed0[,"Fund_ID"],dead_dates_fixed0[,"pull"],dead_dates_fixed0[,"Dead_Date"]),]
row.names(dead_dates_fixed0) <- seq(nrow(dead_dates_fixed0))

dead_dates_fixed <- rbind(dead_dates_fixed0,dead_dates_no)

rm2(dead_dates_fixed0,dead_dates_no,dead_dates_yes)

dead_dates_fixed <- dead_dates_fixed[order(dead_dates_fixed[,"Fund_ID"],
                                           dead_dates_fixed[,"pull"],
                                           dead_dates_fixed[,"Dead_Date"]),]
row.names(dead_dates_fixed) <- seq(nrow(dead_dates_fixed))

###TEST
# dead_dates[dead_dates[,"Fund_ID"]==5248,]
# dead_dates_fixed[dead_dates_fixed[,"Fund_ID"]==5248,]
# dead_dates[dead_dates[,"Fund_ID"]==26603,]
# dead_dates_fixed[dead_dates_fixed[,"Fund_ID"]==26603,]
# dead_dates_u0[dead_dates_u0[,"Fund_ID"]==5248,max(dead_dates_u0[,"Dead_Date"],na.rm=TRUE)]
# dead_dates_u0_dt <- data.table(dead_dates_u0)
# setkeyv(dead_dates_u0_dt, c("Fund_ID"))
# dead_dates_u0_no_multiples <- dead_dates_u0_dt[, max(.SD), by=c("Fund_ID")]
# dead_dates_u0_no_multiples <- as.data.frame(dead_dates_u0_no_multiples,stringsAsFactors=FALSE)
# rm2(dead_dates_u0,dead_dates_u0_dt)
###END TEST


#Expand death dates

dead_dates_fixed_u0 <- unique(dead_dates_fixed[,c("Fund_ID","Dead_Date")])

rm2(dead_dates_fixed)

dead_dates_expand0 <- data.frame(dead_dates_fixed_u0[!is.na(dead_dates_fixed_u0[,c("Dead_Date")]),c("Fund_ID","Dead_Date")],
                                 date=NA,yr=NA,month=NA,flag=NA,stringsAsFactors=FALSE)

dead_dates_expand0  <- dead_dates_expand0[order(dead_dates_expand0[,"Fund_ID"],dead_dates_expand0[,"Dead_Date"]),]
row.names(dead_dates_expand0) <- seq(nrow(dead_dates_expand0))

rm2(dead_dates_fixed_u0)

dead_dates_expand <- coredata(dead_dates_expand0)[rep(seq(nrow(dead_dates_expand0)),length(NAV_AUM_yr_seq)*12),]

dead_dates_expand <- dead_dates_expand[order(dead_dates_expand[,"Fund_ID"]),]
row.names(dead_dates_expand) <- seq(nrow(dead_dates_expand))

dead_dates_yr_month_expand <- expand.grid(yr=NAV_AUM_yr_seq, month=seq(1,12,by=1))
dead_dates_yr_month_expand <- dead_dates_yr_month_expand[order(dead_dates_yr_month_expand[,"yr"],dead_dates_yr_month_expand[,"month"]),]
row.names(dead_dates_yr_month_expand) <- seq(nrow(dead_dates_yr_month_expand))

rm2(NAV_AUM_yr_max,NAV_AUM_yr_min,NAV_AUM_yr_seq)

dead_dates_expand[,c("yr","month")] <- coredata(dead_dates_yr_month_expand)[rep(seq(nrow(dead_dates_yr_month_expand)),nrow(dead_dates_expand0)),]

rm2(dead_dates_yr_month_expand)

dead_dates_expand[,c("date")] <- as.Date(paste(dead_dates_expand[,"yr"],dead_dates_expand[,"month"],"1",sep="-"),format="%Y-%m-%d")
#dead_dates_expand[,c("date")] <- as.POSIXct(paste(dead_dates_expand[,"yr"],dead_dates_expand[,"month"],"1",sep="-"),format="%Y-%m-%d")

dead_dates_expand[,c("flag")] <- ifelse(dead_dates_expand[,"Dead_Date"]>=dead_dates_expand[,"date"],1,0)

dead_dates_expand <- dead_dates_expand[order(dead_dates_expand[,"Fund_ID"],dead_dates_expand[,"yr"],dead_dates_expand[,"month"]),]
row.names(dead_dates_expand) <- seq(nrow(dead_dates_expand))

rm2(dead_dates_expand0)

NAV_AUM_trim0 <- merge(NAV_AUM_rets, dead_dates_expand[,!(colnames(dead_dates_expand) %in% c("date"))], 
                       by.x=c("Fund_ID","yr","month"), by.y=c("Fund_ID","yr","month"), 
                       all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

NAV_AUM_trim0 <- NAV_AUM_trim0[,c(colnames(NAV_AUM_rets),
                                  colnames(NAV_AUM_trim0)[!(colnames(NAV_AUM_trim0) %in% colnames(NAV_AUM_rets))])]

NAV_AUM_trim0 <- NAV_AUM_trim0[order(NAV_AUM_trim0[,"Fund_ID"],
                                     NAV_AUM_trim0[,"yr"],
                                     NAV_AUM_trim0[,"month"]),]
row.names(NAV_AUM_trim0) <- seq(nrow(NAV_AUM_trim0))

rm2(NAV_AUM_rets,dead_dates_expand)

NAV_AUM_trim0[,"flag"] <- ifelse(is.na(NAV_AUM_trim0[,"flag"]),2, NAV_AUM_trim0[,"flag"])

NAV_AUM_trim <- NAV_AUM_trim0[NAV_AUM_trim0[,"flag"] %in% c(1,2),!(colnames(NAV_AUM_trim0) %in% "flag")]

rm2(NAV_AUM_trim0)

NAV_AUM_trim <- NAV_AUM_trim[,c("pull_trim","pull","Fund_ID","Dead_Date","date",
                                colnames(NAV_AUM_trim)[!(colnames(NAV_AUM_trim) %in% c("pull_trim","pull","Fund_ID","Dead_Date","date"))])]

NAV_AUM_trim <- NAV_AUM_trim[order(NAV_AUM_trim[,"Fund_ID"],
                                   NAV_AUM_trim[,"date"],
                                   NAV_AUM_trim[,"pull_trim"],
                                   NAV_AUM_trim[,"pull"]),]
row.names(NAV_AUM_trim) <- seq(nrow(NAV_AUM_trim))


###TEST
# aa <- count(NAV_AUM_trim0[,"Fund_ID"], vars = NULL, wt_var = NULL)
# colnames(aa)[match("freq",names(aa))] <- "a"
# 
# bb <- count(NAV_AUM[,"Fund_ID"], vars = NULL, wt_var = NULL)
# colnames(bb)[match("freq",names(bb))] <- "b"
# 
# cc <- merge(aa, bb, 
#             by.x=c("x"), by.y=c("x"), 
#             all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))
# 
# dd <- data.frame(cc,flag=NA,stringsAsFactors=FALSE)
# dd[,c("flag")] <- ifelse(dd[,"a"]==dd[,"b"],0,1)
# dd[dd[,"flag"]!=0,]
# 
# ee1 <- NAV_AUM_trim0[NAV_AUM_trim0[,"Fund_ID"]==26603,]
# ee2 <- NAV_AUM[NAV_AUM[,"Fund_ID"]==26603,]
# ff <- dead_dates_expand[dead_dates_expand[,"Fund_ID"]==26603,!(colnames(dead_dates_expand) %in% c("Dead_Date","date"))]
#
# gg <- NAV_AUM_trim[NAV_AUM_trim[,"Fund_ID"]==5004,]
###END TEST


###############################################################################
cat("SECTION: OUTPUT FILES", "\n")
###############################################################################

#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

write.csv(Stats_nonreturn, file=paste(final_folder_path,"//","EurekahedgeHF_Stats_noreturns",".csv",sep=""),row.names=FALSE)
write.csv(NAV_AUM_trim, file=paste(final_folder_path,"//","EurekahedgeHF_NAV_AUM_Ret",".csv",sep=""),row.names=FALSE)
write.csv(unique(NAV_AUM_trim[,c("pull_trim","pull","Fund_ID","Dead_Date")]), file=paste(final_folder_path,"//","EurekahedgeHF_Dead_Dates_fixed",".csv",sep=""),row.names=FALSE)

#aaa <- unique(NAV_AUM_trim[,c("pull_trim,"pull","Fund_ID","Dead_Date")])

rm2(final_folder_path,nav_aum_ret_folder_path)

rm2(Stats_nonreturn,NAV_AUM_trim)
rm2(Stats_id_cols,Stats_yearlyret_cols,Stats_monthlyret_cols,Stats_other_cols)


