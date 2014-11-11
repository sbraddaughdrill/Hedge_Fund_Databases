# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_NAV_AUM_merge.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Merge NAV & AUM Eurekahedge Data
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
cat("SECTION: IMPORT NAV & AUM", "\n")
###############################################################################

NAV_AUM_pull <- c("04-2014","02-2013","02-2012","02-2011","02-2010","02-2009","02-2008","02-2007")
NAV_AUM_yr <- c("2014","2013","2012","2011","2010","2009","2008","2007")

NAV_AUM_files <- c("EurekahedgeHF_EXCEL_2014Apr_NAV_AUM","EurekahedgeHF_EXCEL_2013Feb_NAV_AUM",
                   "EurekahedgeHF_EXCEL_2012Feb_NAV_AUM","EurekahedgeHF_EXCEL_2011Feb_NAV_AUM",
                   "EurekahedgeHF_EXCEL_2010Feb_NAV_AUM","EurekahedgeHF_EXCEL_2009Feb_NAV_AUM",
                   "EurekahedgeHF_EXCEL_2008Feb_NAV_AUM","EurekahedgeHF_EXCEL_2007Feb_NAV_AUM")

NAV_AUM_input <- data.frame(matrix(NA, ncol=3, nrow=length(NAV_AUM_files), dimnames=list(c(), c("pull","yr","file"))), 
                            stringsAsFactors=FALSE)

NAV_AUM_input[,"pull"] <- NAV_AUM_pull
NAV_AUM_input[,"yr"] <- NAV_AUM_yr
NAV_AUM_input[,"file"] <- NAV_AUM_files

rm2(NAV_AUM_pull,NAV_AUM_yr,NAV_AUM_files)


###############################################################################
cat("SECTION: CONCATENATE NAV & AUM", "\n")
###############################################################################

melt_folder_path <- paste(output_directory, "NAV_AUM_melt", sep = "//", collapse = "//")  

NAV_AUM_id_cols <- c("Ret_AUM","Fund.ID","Fund.Name")

NAV_AUM_concatentate0 <- alply(.data=NAV_AUM_input, .margins=1, .fun = function(x,directory_in){
  
  # x <- NAV_AUM_input[7,]
  # directory_in <- melt_folder_path
  
  input <- data.frame(read.csv(file=paste(directory_in,"\\",x[,"file"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      yr=NA,month=NA,
                      stringsAsFactors=FALSE)
  #bad <- input[is.na(input[,"Fund.ID"]),]
  
  return(input)
  
}, directory_in=melt_folder_path, .expand = TRUE, .progress = "text")

NAV_AUM_concatentate <- rbind.fill(NAV_AUM_concatentate0)

rm2(NAV_AUM_concatentate0)

for(i in which(sapply(NAV_AUM_concatentate,class)=="character"))
{
  NAV_AUM_concatentate[[i]] = trim(NAV_AUM_concatentate[[i]])
}
rm2(i)
for (i in 1:ncol(NAV_AUM_concatentate))
{
  NAV_AUM_concatentate[,i] <- unknownToNA(NAV_AUM_concatentate[,i], unknown=unknowns_strings,force=TRUE)
  NAV_AUM_concatentate[,i] <- ifelse(is.na(NAV_AUM_concatentate[,i]),NA, NAV_AUM_concatentate[,i])
} 
rm2(i)

rm2(NAV_AUM_input,NAV_AUM_id_cols)


#Names_concatentate <- NAV_AUM_concatentate[,c("pull","Fund.ID","Fund.Name")]
#Names_concatentate_u <- unique(Names_concatentate)
#row.names(Names_concatentate_u) <- seq(nrow(Names_concatentate_u))
#write.csv(Names_concatentate_u,file=paste(output_directory,"//","Names_concatentate_u",".csv",sep=""),row.names=FALSE)

NAV_AUM_concatentate <- NAV_AUM_concatentate[,c("pull","Fund.ID","Fund.Name","date","yr","month","Monthly_Ret","AUM")]

#NAV_AUM_concatentate <- NAV_AUM_concatentate[rowSums(is.na(NAV_AUM_concatentate[,1:ncol(NAV_AUM_concatentate)]))<ncol(NAV_AUM_concatentate),]

NAV_AUM_concatentate[,"date"] <- as.yearmon(NAV_AUM_concatentate[,"date"],format="%b %Y")
NAV_AUM_concatentate[,"date"] <- as.Date(NAV_AUM_concatentate[,"date"],format="%m-%d-%Y")

NAV_AUM_concatentate[,"yr"] <- year(NAV_AUM_concatentate[,"date"])
NAV_AUM_concatentate[,"month"] <- month(NAV_AUM_concatentate[,"date"])

NAV_AUM_concatentate <- NAV_AUM_concatentate[order(NAV_AUM_concatentate[,"pull"],NAV_AUM_concatentate[,"Fund.ID"],
                                                   NAV_AUM_concatentate[,"Fund.Name"], NAV_AUM_concatentate[,"date"]),]

#NAV_AUM_concatentate <- NAV_AUM_concatentate[order(NAV_AUM_concatentate[,"pull"],NAV_AUM_concatentate[,"Fund.ID"],
#                                                    NAV_AUM_concatentate[,"Fund.Name"],NAV_AUM_concatentate[,"yr"],
#                                                    NAV_AUM_concatentate[,"month"]),]

row.names(NAV_AUM_concatentate) <- seq(nrow(NAV_AUM_concatentate))

write.csv(NAV_AUM_concatentate[,c("pull","Fund.ID","Fund.Name","date","Monthly_Ret","AUM")], 
          file=paste(output_directory,"//","EurekahedgeHF_NAV_AUM_melt",".csv",sep=""),row.names=FALSE)



















###############################################################################
cat("SECTION: IMPORT STATS", "\n")
###############################################################################

Import_Stats_yr <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Stats","EurekahedgeHF_EXCEL_2013Feb_Data_Stats",
                     "EurekahedgeHF_EXCEL_2012Feb_Data_Stats","EurekahedgeHF_EXCEL_2011Feb_Data_Stats",
                     "EurekahedgeHF_EXCEL_2010Feb_Data_Stats","EurekahedgeHF_EXCEL_2009Feb_Data_Stats",
                     "EurekahedgeHF_EXCEL_2008Feb_Data_Stats","EurekahedgeHF_EXCEL_2007Feb_Data_Stats")

Import_Stats_files <- c("EurekahedgeHF_EXCEL_2014Apr_Data_Stats","EurekahedgeHF_EXCEL_2013Feb_Data_Stats",
                        "EurekahedgeHF_EXCEL_2012Feb_Data_Stats","EurekahedgeHF_EXCEL_2011Feb_Data_Stats",
                        "EurekahedgeHF_EXCEL_2010Feb_Data_Stats","EurekahedgeHF_EXCEL_2009Feb_Data_Stats",
                        "EurekahedgeHF_EXCEL_2008Feb_Data_Stats","EurekahedgeHF_EXCEL_2007Feb_Data_Stats")

Import_Stats_input <- data.frame(matrix(NA, ncol=2, nrow=length(Import_Stats_files), dimnames=list(c(), c("yr","file"))), 
                                 stringsAsFactors=FALSE)









###############################################################################
cat("SECTION: CLEAN EurekahedgeHF_Excel_aca_Instruments_Traded", "\n")
###############################################################################

EurekahedgeHF_Excel_aca_Instruments_Traded <- read.csv(file=paste(output_directory,files[3,1],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)
for(i in which(sapply(EurekahedgeHF_Excel_aca_Instruments_Traded,class)=="character"))
{
  EurekahedgeHF_Excel_aca_Instruments_Traded[[i]] = trim(EurekahedgeHF_Excel_aca_Instruments_Traded[[i]])
}
for (i in 1:ncol(EurekahedgeHF_Excel_aca_Instruments_Traded))
{
  EurekahedgeHF_Excel_aca_Instruments_Traded[,i] <- unknownToNA(EurekahedgeHF_Excel_aca_Instruments_Traded[,i], unknown=c("",".","n/a","na","NA",NA,"null","NULL",NULL,"nan","NaN",NaN,
                                                                                                                          NA_integer_,"NA_integer_",NA_complex_,"NA_complex_",
                                                                                                                          NA_character_,"NA_character_",NA_real_,"NA_real_"),force=TRUE)
  EurekahedgeHF_Excel_aca_Instruments_Traded[,i] <- ifelse(is.na(EurekahedgeHF_Excel_aca_Instruments_Traded[,i]),NA, EurekahedgeHF_Excel_aca_Instruments_Traded[,i])
} 

EurekahedgeHF_Excel_aca_Instruments_Traded  <- EurekahedgeHF_Excel_aca_Instruments_Traded[order(EurekahedgeHF_Excel_aca_Instruments_Traded[,"Fund.ID"],
                                                                                                EurekahedgeHF_Excel_aca_Instruments_Traded[,"Fund.Name"],
                                                                                                EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"]),]

row.names(EurekahedgeHF_Excel_aca_Instruments_Traded) <- seq(nrow(EurekahedgeHF_Excel_aca_Instruments_Traded))

EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"] <-  gsub(pattern=" ", replacement=".", x=EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"])
EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"] <-  gsub(pattern="\\.{2,}", replacement="\\.", x=EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"])
EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"] <-  gsub(pattern="-", replacement="_", x=EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"])

unique_instruments <-  unique(EurekahedgeHF_Excel_aca_Instruments_Traded[,"Instrument.Traded"])

#Instruments Traded
Instruments_Traded <- data.frame(EurekahedgeHF_Excel_aca_Instruments_Traded, 
                                 matrix(0, ncol=length(unique_instruments), nrow=nrow(EurekahedgeHF_Excel_aca_Instruments_Traded), dimnames=list(c(), paste("Instrument.Traded",unique_instruments,sep="_"))), 
                                 stringsAsFactors=FALSE)

for (i in 1:length(unique_instruments))
{
  
  Instruments_Traded[,paste("Instrument.Traded",unique_instruments[i],sep="_")] <- 
    ifelse(Instruments_Traded[,"Instrument.Traded"]==unique_instruments[i], 1, Instruments_Traded[,paste("Instrument.Traded",unique_instruments[i],sep="_")])
  
} 

Instruments_Traded_comb <- aggregate(Instruments_Traded[,(ncol(Instruments_Traded)-length(unique_instruments)+1):ncol(Instruments_Traded)], by=list(Instruments_Traded[,"Fund.ID"]), FUN=sum, na.rm=TRUE)
colnames(Instruments_Traded_comb)[1] <- "Fund.ID"

rm2(Instruments_Traded)

Instruments_Traded_comb <- Instruments_Traded_comb[,sort(colnames(Instruments_Traded_comb))]
Instruments_Traded_comb <- Instruments_Traded_comb[,c("Fund.ID",colnames(Instruments_Traded_comb)[-which(colnames(Instruments_Traded_comb) %in% "Fund.ID")])]

#Exposure
Exposure <- data.frame(EurekahedgeHF_Excel_aca_Instruments_Traded, 
                       matrix(NA, ncol=length(unique_instruments), nrow=nrow(EurekahedgeHF_Excel_aca_Instruments_Traded), dimnames=list(c(), paste("Exposure",unique_instruments,sep="_"))), 
                       stringsAsFactors=FALSE)

rm2(EurekahedgeHF_Excel_aca_Instruments_Traded)

for (i in 1:length(unique_instruments))
{
  
  Exposure[,paste("Exposure",unique_instruments[i],sep="_")] <- 
    ifelse(Exposure[,"Instrument.Traded"]==unique_instruments[i], Exposure[,"Exposure"], Exposure[,paste("Exposure",unique_instruments[i],sep="_")])
  
} 

Exposure_comb <- dcast(Exposure, Fund.ID ~ Instrument.Traded, value.var = 'Exposure')
colnames(Exposure_comb)[2:ncol(Exposure_comb)] <- paste("Exposure",colnames(Exposure_comb)[2:ncol(Exposure_comb)],sep="_")
colnames(Exposure_comb)[1] <- "Fund.ID"

rm2(Exposure)

Exposure_comb <- Exposure_comb[,sort(colnames(Exposure_comb))]
Exposure_comb <- Exposure_comb[,c("Fund.ID",colnames(Exposure_comb)[-which(colnames(Exposure_comb) %in% "Fund.ID")])]


EurekahedgeHF_Excel_aca_Instruments_Traded_merge <- merge(Instruments_Traded_comb, Exposure_comb, by.x=c("Fund.ID"), by.y=c("Fund.ID"), 
                                                          all.x=TRUE, all.y=FALSE, sort=TRUE, suffixes=c(".x",".y"))

rm(Instruments_Traded_comb,Exposure_comb)

EurekahedgeHF_Excel_aca_Instruments_Traded_merge <- EurekahedgeHF_Excel_aca_Instruments_Traded_merge[rowSums(is.na(EurekahedgeHF_Excel_aca_Instruments_Traded_merge[,1:ncol(EurekahedgeHF_Excel_aca_Instruments_Traded_merge)]))<ncol(EurekahedgeHF_Excel_aca_Instruments_Traded_merge),]

EurekahedgeHF_Excel_aca_Instruments_Traded_merge <- EurekahedgeHF_Excel_aca_Instruments_Traded_merge[order(EurekahedgeHF_Excel_aca_Instruments_Traded_merge[,"Fund.ID"]),]

row.names(EurekahedgeHF_Excel_aca_Instruments_Traded_merge) <- seq(nrow(EurekahedgeHF_Excel_aca_Instruments_Traded_merge))


###############################################################################
cat("SECTION: CLEAN EurekahedgeHF_Excel_aca", "\n")
###############################################################################

EurekahedgeHF_Excel_aca <- read.csv(file=paste(output_directory,files[1,1],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE)
for(i in which(sapply(EurekahedgeHF_Excel_aca,class)=="character"))
{
  EurekahedgeHF_Excel_aca[[i]] = trim(EurekahedgeHF_Excel_aca[[i]])
}
for (i in 1:ncol(EurekahedgeHF_Excel_aca))
{
  EurekahedgeHF_Excel_aca[,i] <- unknownToNA(EurekahedgeHF_Excel_aca[,i], unknown=c("",".","n/a","na","NA",NA,"null","NULL",NULL,"nan","NaN",NaN,
                                                                                    NA_integer_,"NA_integer_",NA_complex_,"NA_complex_",
                                                                                    NA_character_,"NA_character_",NA_real_,"NA_real_"),force=TRUE)
  EurekahedgeHF_Excel_aca[,i] <- ifelse(is.na(EurekahedgeHF_Excel_aca[,i]),NA, EurekahedgeHF_Excel_aca[,i])
} 

EurekahedgeHF_Excel_aca_cols <- colnames(EurekahedgeHF_Excel_aca)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="\\.{2,}", replacement="\\.", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="[[:punct:]]?$", replacement="", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="X2011.Returns", replacement="Returns.2011", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="X2012.Returns", replacement="Returns.2012", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="VaR.90", replacement="VaR.90pct", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="VaR.95", replacement="VaR.95pct", x=EurekahedgeHF_Excel_aca_cols)
EurekahedgeHF_Excel_aca_cols <- gsub(pattern="VaR.99", replacement="VaR.99pct", x=EurekahedgeHF_Excel_aca_cols)
colnames(EurekahedgeHF_Excel_aca) <- EurekahedgeHF_Excel_aca_cols

EurekahedgeHF_Excel_aca[,"Date.Added"] <- gsub(pattern="/", replacement="-", x=EurekahedgeHF_Excel_aca[,"Date.Added"])
EurekahedgeHF_Excel_aca[,"Date.Added"] <- as.Date(EurekahedgeHF_Excel_aca[,"Date.Added"],format="%m-%d-%Y")

EurekahedgeHF_Excel_aca[,"Dead.Date"] <- gsub(pattern="/", replacement="-", x=EurekahedgeHF_Excel_aca[,"Dead.Date"])
EurekahedgeHF_Excel_aca[,"Dead.Date"] <- as.Date(EurekahedgeHF_Excel_aca[,"Dead.Date"],format="%m-%d-%Y")

EurekahedgeHF_Excel_aca[,"Inception.Date"] <- gsub(pattern="/", replacement="-", x=EurekahedgeHF_Excel_aca[,"Inception.Date"])
EurekahedgeHF_Excel_aca[,"Inception.Date"] <- as.Date(EurekahedgeHF_Excel_aca[,"Inception.Date"],format="%m-%d-%Y")

EurekahedgeHF_Excel_aca <- EurekahedgeHF_Excel_aca[rowSums(is.na(EurekahedgeHF_Excel_aca[,1:ncol(EurekahedgeHF_Excel_aca)]))<ncol(EurekahedgeHF_Excel_aca),]

EurekahedgeHF_Excel_aca <- EurekahedgeHF_Excel_aca[order(EurekahedgeHF_Excel_aca[,"Fund.ID"]),]

row.names(EurekahedgeHF_Excel_aca) <- seq(nrow(EurekahedgeHF_Excel_aca))

write.csv(EurekahedgeHF_Excel_aca, file=paste(output_directory,file="EurekahedgeHF_Excel_aca",".csv",sep=""),row.names=FALSE)

rm2(EurekahedgeHF_Excel_aca_cols)


###############################################################################
cat("SECTION: REMOVE MONTHLY RETURNS FROM EurekahedgeHF_Excel_aca", "\n")
###############################################################################

monthly_ret_cols <- c("Jun.12.Returns","May.12.Returns","Apr.12.Returns")

EurekahedgeHF_Excel_aca_monthly_ret_temp <- EurekahedgeHF_Excel_aca[,c("Fund.ID",monthly_ret_cols)]

EurekahedgeHF_Excel_aca_monthly_ret_temp  <- EurekahedgeHF_Excel_aca_monthly_ret_temp[order(EurekahedgeHF_Excel_aca_monthly_ret_temp[,"Fund.ID"]),]

row.names(EurekahedgeHF_Excel_aca_monthly_ret_temp) <- seq(nrow(EurekahedgeHF_Excel_aca_monthly_ret_temp))

# for(i in 1:length(monthly_ret_cols))
# {
#   #i <- 1
#   #i <- 2
#   #i <- 3
#   
#   #temp <-  melt(EurekahedgeHF_Excel_aca_monthly_ret_temp[,c("Fund.ID",monthly_ret_cols[i])], id=c("Fund.ID"), na.rm=FALSE)
#   temp <- data.frame(melt(EurekahedgeHF_Excel_aca_monthly_ret_temp[,c("Fund.ID",monthly_ret_cols[i])], id=c("Fund.ID"), na.rm=FALSE), 
#                      yr=NA, 
#                      month=NA, stringsAsFactors=FALSE)
#   
#   for(j in 1:ncol(temp))
#   {
#     #j <- 1
#     
#     temp[,j] = trim(temp[,j])
#     
#   }
#   
#   temp[,"Fund.ID"] <- as.integer(temp[,"Fund.ID"])
#   
#   colnames(temp)[match("variable",names(temp))] <- "date"
#   colnames(temp)[match("value",names(temp))] <- "Monthly_Ret2"
#   
#   temp[,"date"] <- gsub(pattern="X", replacement="", x=temp[,"date"])
#   temp[,"date"] <- gsub(pattern="\\.", replacement="-", x=temp[,"date"])
#   temp[,"date"] <- gsub(pattern=" ", replacement="", x=temp[,"date"])
#   temp[,"date"] <- gsub(pattern=" ", replacement="", x=temp[,"date"])
#   temp[,"date"] <- gsub(pattern="-Returns", replacement="", x=temp[,"date"])
#   
#   for(j in 1:ncol(temp))
#   {
#     #j <- 1
#     
#     temp[,j] = trim(temp[,j])
#     
#   }
#   
#   temp_month <- temp[,"date"]
#   temp_month <- substr(temp_month, 1, 3)
#   temp_month <- match(tolower(temp_month), tolower(month.abb))
#   
#   temp_yr <- temp[,"date"]
#   temp_yr <- substr(temp_yr, 5, 6)
#   temp_yr <- as.integer(temp_yr)
#   
#   temp_dt <- format(as.Date(paste(temp_yr,temp_month,"01",sep="-"),format="%y-%m-%d"),"%Y-%m-%d")
#   temp_dt <- as.Date(temp_dt,format="%Y-%m-%d")
#   
#   temp[,"date"] <- temp_dt
#   
#   temp[,"Monthly_Ret2"] <- as.numeric(temp[,"Monthly_Ret2"])
#   
#   cat("LOOP: ",i, "\n")
#   
#   if(i==1)
#   {
#     EurekahedgeHF_Excel_aca_monthly_ret <- temp
#     
#   } else
#   {
#     EurekahedgeHF_Excel_aca_monthly_ret <- rbind(EurekahedgeHF_Excel_aca_monthly_ret,temp)
#   }
#   
# }
# 
# EurekahedgeHF_Excel_aca_monthly_ret[,"yr"] <- year(EurekahedgeHF_Excel_aca_monthly_ret[,"date"])
# EurekahedgeHF_Excel_aca_monthly_ret[,"month"] <- month(EurekahedgeHF_Excel_aca_monthly_ret[,"date"])

EurekahedgeHF_Excel_aca_monthly_ret <- ldply(.data=monthly_ret_cols, .fun = function(x,data){
  
  temp <- data.frame(melt(data[,c("Fund.ID",x)], id=c("Fund.ID"), na.rm=FALSE),  yr=NA, month=NA, stringsAsFactors=FALSE)
  
  return(temp)
  
}, data=EurekahedgeHF_Excel_aca_monthly_ret_temp, 
.progress = "none", .inform = FALSE,.parallel = FALSE, .paropts = NULL, .id = NA)

rm(EurekahedgeHF_Excel_aca_monthly_ret_temp)

for(j in 1:ncol(EurekahedgeHF_Excel_aca_monthly_ret))
{
  #j <- 1
  
  EurekahedgeHF_Excel_aca_monthly_ret[,j] = trim(EurekahedgeHF_Excel_aca_monthly_ret[,j])
  
}

EurekahedgeHF_Excel_aca_monthly_ret[,"Fund.ID"] <- as.integer(EurekahedgeHF_Excel_aca_monthly_ret[,"Fund.ID"])

colnames(EurekahedgeHF_Excel_aca_monthly_ret)[match("variable",names(EurekahedgeHF_Excel_aca_monthly_ret))] <- "date"
colnames(EurekahedgeHF_Excel_aca_monthly_ret)[match("value",names(EurekahedgeHF_Excel_aca_monthly_ret))] <- "Monthly_Ret2"

EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- gsub(pattern="X", replacement="", x=EurekahedgeHF_Excel_aca_monthly_ret[,"date"])
EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- gsub(pattern="\\.", replacement="-", x=EurekahedgeHF_Excel_aca_monthly_ret[,"date"])
EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- gsub(pattern=" ", replacement="", x=EurekahedgeHF_Excel_aca_monthly_ret[,"date"])
EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- gsub(pattern=" ", replacement="", x=EurekahedgeHF_Excel_aca_monthly_ret[,"date"])
EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- gsub(pattern="-Returns", replacement="", x=EurekahedgeHF_Excel_aca_monthly_ret[,"date"])

for(j in 1:ncol(EurekahedgeHF_Excel_aca_monthly_ret))
{
  #j <- 1
  
  EurekahedgeHF_Excel_aca_monthly_ret[,j] = trim(EurekahedgeHF_Excel_aca_monthly_ret[,j])
  
}

EurekahedgeHF_Excel_aca_monthly_ret[,"month"] <- EurekahedgeHF_Excel_aca_monthly_ret[,"date"]
EurekahedgeHF_Excel_aca_monthly_ret[,"month"] <- substr(EurekahedgeHF_Excel_aca_monthly_ret[,"month"], 1, 3)
EurekahedgeHF_Excel_aca_monthly_ret[,"month"] <- match(tolower(EurekahedgeHF_Excel_aca_monthly_ret[,"month"]), tolower(month.abb))

EurekahedgeHF_Excel_aca_monthly_ret[,"yr"] <- EurekahedgeHF_Excel_aca_monthly_ret[,"date"]
EurekahedgeHF_Excel_aca_monthly_ret[,"yr"] <- substr(EurekahedgeHF_Excel_aca_monthly_ret[,"yr"], 5, 6)
EurekahedgeHF_Excel_aca_monthly_ret[,"yr"] <- as.integer(EurekahedgeHF_Excel_aca_monthly_ret[,"yr"])

EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- format(as.Date(paste(EurekahedgeHF_Excel_aca_monthly_ret[,"yr"],
                                                                     EurekahedgeHF_Excel_aca_monthly_ret[,"month"],
                                                                     "01",sep="-"),
                                                               format="%y-%m-%d"),"%Y-%m-%d")
EurekahedgeHF_Excel_aca_monthly_ret[,"date"] <- as.Date(EurekahedgeHF_Excel_aca_monthly_ret[,"date"],format="%Y-%m-%d")

EurekahedgeHF_Excel_aca_monthly_ret[,"Monthly_Ret2"] <- as.numeric(EurekahedgeHF_Excel_aca_monthly_ret[,"Monthly_Ret2"])

EurekahedgeHF_Excel_aca_monthly_ret[,"yr"] <- year(EurekahedgeHF_Excel_aca_monthly_ret[,"date"])
EurekahedgeHF_Excel_aca_monthly_ret[,"month"] <- month(EurekahedgeHF_Excel_aca_monthly_ret[,"date"])

EurekahedgeHF_Excel_aca_monthly_ret <- EurekahedgeHF_Excel_aca_monthly_ret[rowSums(is.na(EurekahedgeHF_Excel_aca_monthly_ret[,1:ncol(EurekahedgeHF_Excel_aca_monthly_ret)]))<ncol(EurekahedgeHF_Excel_aca_monthly_ret),]

EurekahedgeHF_Excel_aca_monthly_ret <- EurekahedgeHF_Excel_aca_monthly_ret[order(EurekahedgeHF_Excel_aca_monthly_ret[,"Fund.ID"],
                                                                                 EurekahedgeHF_Excel_aca_monthly_ret[,"date"],
                                                                                 EurekahedgeHF_Excel_aca_monthly_ret[,"yr"],
                                                                                 EurekahedgeHF_Excel_aca_monthly_ret[,"month"]),]

row.names(EurekahedgeHF_Excel_aca_monthly_ret) <- seq(nrow(EurekahedgeHF_Excel_aca_monthly_ret))


###############################################################################
cat("SECTION: REMOVE YEARLY RETURNS FROM EurekahedgeHF_Excel_aca", "\n")
###############################################################################

yearly_ret_cols <- c("Returns.2011","Returns.2012")

EurekahedgeHF_Excel_aca_yearly_ret_temp <- EurekahedgeHF_Excel_aca[,c("Fund.ID",yearly_ret_cols)]

EurekahedgeHF_Excel_aca_yearly_ret_temp  <- EurekahedgeHF_Excel_aca_yearly_ret_temp[order(EurekahedgeHF_Excel_aca_yearly_ret_temp[,"Fund.ID"]),]

row.names(EurekahedgeHF_Excel_aca_yearly_ret_temp) <- seq(nrow(EurekahedgeHF_Excel_aca_yearly_ret_temp))

# for(i in 1:length(yearly_ret_cols))
# {
#   #i <- 1
#   #i <- 2
#   #i <- 3
#   
#   temp <-  melt(EurekahedgeHF_Excel_aca_yearly_ret_temp[,c("Fund.ID",yearly_ret_cols[i])], id=c("Fund.ID"), na.rm=FALSE)
#   
#   for(j in 1:ncol(temp))
#   {
#     #j <- 1
#     
#     temp[,j] = trim(temp[,j])
#     
#   }
#   
#   temp[,"Fund.ID"] <- as.integer(temp[,"Fund.ID"])
#   
#   colnames(temp)[match("variable",names(temp))] <- "yr"
#   colnames(temp)[match("value",names(temp))] <- "Yearly_Ret"
#   
#   temp[,"yr"] <- gsub(pattern="X", replacement="", x=temp[,"yr"])
#   temp[,"yr"] <- gsub(pattern="\\.", replacement="-", x=temp[,"yr"])
#   temp[,"yr"] <- gsub(pattern=" ", replacement="", x=temp[,"yr"])
#   temp[,"yr"] <- gsub(pattern=" ", replacement="", x=temp[,"yr"])
#   temp[,"yr"] <- gsub(pattern="Returns-", replacement="", x=temp[,"yr"])
#   
#   for(j in 1:ncol(temp))
#   {
#     #j <- 1
#     
#     temp[,j] = trim(temp[,j])
#     
#   }
#   
#   temp[,"yr"] <- as.integer(temp[,"yr"])
#   temp[,"Yearly_Ret"] <- as.numeric(temp[,"Yearly_Ret"])
#   
#   cat("LOOP: ",i, "\n")
#   
#   if(i==1)
#   {
#     EurekahedgeHF_Excel_aca_yearly_ret <- temp
#     
#   } else
#   {
#     EurekahedgeHF_Excel_aca_yearly_ret <- rbind(EurekahedgeHF_Excel_aca_yearly_ret,temp)
#   }
#   
# }

EurekahedgeHF_Excel_aca_yearly_ret <- ldply(.data=yearly_ret_cols, .fun = function(x,data){
  
  temp <- data.frame(melt(data[,c("Fund.ID",x)], id=c("Fund.ID"), na.rm=FALSE), stringsAsFactors=FALSE)
  
  return(temp)
  
}, data=EurekahedgeHF_Excel_aca_yearly_ret_temp, 
.progress = "none", .inform = FALSE,.parallel = FALSE, .paropts = NULL, .id = NA)

rm(EurekahedgeHF_Excel_aca_yearly_ret_temp)

for(j in 1:ncol(EurekahedgeHF_Excel_aca_yearly_ret))
{
  #j <- 1
  
  EurekahedgeHF_Excel_aca_yearly_ret[,j] = trim(EurekahedgeHF_Excel_aca_yearly_ret[,j])
  
}

EurekahedgeHF_Excel_aca_yearly_ret[,"Fund.ID"] <- as.integer(EurekahedgeHF_Excel_aca_yearly_ret[,"Fund.ID"])

colnames(EurekahedgeHF_Excel_aca_yearly_ret)[match("variable",names(EurekahedgeHF_Excel_aca_yearly_ret))] <- "yr"
colnames(EurekahedgeHF_Excel_aca_yearly_ret)[match("value",names(EurekahedgeHF_Excel_aca_yearly_ret))] <- "Yearly_Ret"

EurekahedgeHF_Excel_aca_yearly_ret[,"yr"] <- gsub(pattern="X", replacement="", x=EurekahedgeHF_Excel_aca_yearly_ret[,"yr"])
EurekahedgeHF_Excel_aca_yearly_ret[,"yr"] <- gsub(pattern="\\.", replacement="-", x=EurekahedgeHF_Excel_aca_yearly_ret[,"yr"])
EurekahedgeHF_Excel_aca_yearly_ret[,"yr"] <- gsub(pattern=" ", replacement="", x=EurekahedgeHF_Excel_aca_yearly_ret[,"yr"])
EurekahedgeHF_Excel_aca_yearly_ret[,"yr"] <- gsub(pattern=" ", replacement="", x=EurekahedgeHF_Excel_aca_yearly_ret[,"yr"])
EurekahedgeHF_Excel_aca_yearly_ret[,"yr"] <- gsub(pattern="Returns-", replacement="", x=EurekahedgeHF_Excel_aca_yearly_ret[,"yr"])

for(j in 1:ncol(EurekahedgeHF_Excel_aca_yearly_ret))
{
  #j <- 1
  
  EurekahedgeHF_Excel_aca_yearly_ret[,j] = trim(EurekahedgeHF_Excel_aca_yearly_ret[,j])
  
}

EurekahedgeHF_Excel_aca_yearly_ret[,"yr"] <- as.integer(EurekahedgeHF_Excel_aca_yearly_ret[,"yr"])
EurekahedgeHF_Excel_aca_yearly_ret[,"Yearly_Ret"] <- as.numeric(EurekahedgeHF_Excel_aca_yearly_ret[,"Yearly_Ret"])

EurekahedgeHF_Excel_aca_yearly_ret <- EurekahedgeHF_Excel_aca_yearly_ret[rowSums(is.na(EurekahedgeHF_Excel_aca_yearly_ret[,1:ncol(EurekahedgeHF_Excel_aca_yearly_ret)]))<ncol(EurekahedgeHF_Excel_aca_yearly_ret),]

EurekahedgeHF_Excel_aca_yearly_ret <- EurekahedgeHF_Excel_aca_yearly_ret[order(EurekahedgeHF_Excel_aca_yearly_ret[,"Fund.ID"],
                                                                               EurekahedgeHF_Excel_aca_yearly_ret[,"yr"]),]

row.names(EurekahedgeHF_Excel_aca_yearly_ret) <- seq(nrow(EurekahedgeHF_Excel_aca_yearly_ret))


###############################################################################
cat("SECTION: MERGE DATA", "\n")
###############################################################################

#EurekahedgeHF_Excel_aca_NAV_AUM_melt[,!names(EurekahedgeHF_Excel_aca_NAV_AUM_melt) %in% c("date")]
EurekahedgeHF_Excel_aca_full0 <- merge(EurekahedgeHF_Excel_aca_NAV_AUM_melt,
                                       EurekahedgeHF_Excel_aca_monthly_ret[,!names(EurekahedgeHF_Excel_aca_monthly_ret) %in% c("date")], 
                                       by.x=c("Fund.ID","yr","month"), by.y=c("Fund.ID","yr","month"), 
                                       all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(EurekahedgeHF_Excel_aca_NAV_AUM_melt,EurekahedgeHF_Excel_aca_monthly_ret)

EurekahedgeHF_Excel_aca_full1 <- merge(EurekahedgeHF_Excel_aca_full0, 
                                       EurekahedgeHF_Excel_aca_yearly_ret, 
                                       by.x=c("Fund.ID","yr"), by.y=c("Fund.ID","yr"), 
                                       all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(EurekahedgeHF_Excel_aca_full0,EurekahedgeHF_Excel_aca_yearly_ret)

EurekahedgeHF_Excel_aca_full2 <- merge(EurekahedgeHF_Excel_aca[,!names(EurekahedgeHF_Excel_aca) %in% c(monthly_ret_cols,yearly_ret_cols)],
                                       EurekahedgeHF_Excel_aca_full1, 
                                       by.x=c("Fund.ID"), by.y=c("Fund.ID"), 
                                       all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(EurekahedgeHF_Excel_aca,EurekahedgeHF_Excel_aca_full1)

EurekahedgeHF_Excel_aca_full3 <- merge(EurekahedgeHF_Excel_aca_full2, EurekahedgeHF_Excel_aca_Instruments_Traded_merge, 
                                       by.x=c("Fund.ID"), by.y=c("Fund.ID"), 
                                       all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(EurekahedgeHF_Excel_aca_full2,EurekahedgeHF_Excel_aca_Instruments_Traded_merge)

EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_full3[!is.na(EurekahedgeHF_Excel_aca_full3[,"Strategy"]),]
EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[EurekahedgeHF_Excel_aca_merge[,"Base.Currency"]=="USD",]
EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[EurekahedgeHF_Excel_aca_merge[,"Minimum.Investment.Currency"]=="USD",]

rm2(EurekahedgeHF_Excel_aca_full3)

EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[rowSums(is.na(EurekahedgeHF_Excel_aca_merge[,1:ncol(EurekahedgeHF_Excel_aca_merge)]))<ncol(EurekahedgeHF_Excel_aca_merge),]

EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[order(EurekahedgeHF_Excel_aca_merge[,"Fund.ID"],
                                                                     EurekahedgeHF_Excel_aca_merge[,"Fund.Name"],
                                                                     EurekahedgeHF_Excel_aca_merge[,"date"],
                                                                     EurekahedgeHF_Excel_aca_merge[,"yr"],                                                                  
                                                                     EurekahedgeHF_Excel_aca_merge[,"month"]),]

row.names(EurekahedgeHF_Excel_aca_merge) <- seq(nrow(EurekahedgeHF_Excel_aca_merge))


###############################################################################
cat("REORDER COLUMNS", "\n")
###############################################################################

#EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[!is.na(EurekahedgeHF_Excel_aca_merge[,"Strategy"]),]
#EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[EurekahedgeHF_Excel_aca_merge[,"Base.Currency"]=="USD",]

starting_cols <- c("Fund.ID","Fund.Name","Date.Added","Flagship","Closed","Limited","Dead","Dead.Date","Dead.Reason",
                   "Eurekahedge.ID","ISIN","SEDOL","Valoren","CUSIP","Bloomberg","Reuters",
                   "date","yr","month","AUM","Yearly_Ret","Monthly_Ret","Monthly_Ret2")

all_cols <- colnames(EurekahedgeHF_Excel_aca_merge)

other_cols <- all_cols[-which(all_cols %in% starting_cols)]

EurekahedgeHF_Excel_aca_merge <- EurekahedgeHF_Excel_aca_merge[,c(starting_cols,other_cols)]

write.csv(EurekahedgeHF_Excel_aca_merge, file=paste(output_directory,file="EurekahedgeHF_Excel_aca_merge",".csv",sep=""),row.names=FALSE)

EurekahedgeHF_Excel_aca_merge_trim <- data.frame(EurekahedgeHF_Excel_aca_merge[,c("Fund.ID","Fund.Name","date","yr","month","Strategy")],
                                                 stringsAsFactors=FALSE)

rm2(EurekahedgeHF_Excel_aca_merge)

EurekahedgeHF_Excel_aca_merge_trim <- unique(EurekahedgeHF_Excel_aca_merge_trim[,c("Fund.ID","Fund.Name","yr","Strategy")])

write.csv(EurekahedgeHF_Excel_aca_merge_trim, file=paste(output_directory,file="EurekahedgeHF_Excel_aca_merge_trim",".csv",sep=""),row.names=FALSE)

rm2(EurekahedgeHF_Excel_aca_merge_trim)
