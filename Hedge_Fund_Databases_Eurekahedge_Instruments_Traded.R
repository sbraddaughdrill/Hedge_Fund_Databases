# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Instruments_Traded.R
# Version: 1.0
# Date:    11.10.2014
# Purpose: Combine Instruments Traded Eurekahedge Data
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
cat("SECTION: IMPORT INSTRUMENTS", "\n")
###############################################################################

Instruments_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Instruments_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)

# Instruments_input <- data.frame(pull=NA,file_name=NA,read.csv(file=paste(output_directory,"\\","EurekahedgeHF_Instruments_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
#                             stringsAsFactors=FALSE)
# Instruments_input[,"file_name"] <- Instruments_input[,"file_clean"]
# Instruments_input[,"file_name"] <- gsub("\\\\","/",Instruments_input[,"file_name"])
# Instruments_input[,"file_name"] <- gsub("//","/",Instruments_input[,"file_name"])
# Instruments_input[,"file_name"] <- gsub("//","/",Instruments_input[,"file_name"])
# 
# Instruments_input[,"file_name"] <- encodeString(Instruments_input[,"file_name"])
# 
# #Instruments_input[,"pull"] <- regexpr("/[^/]*$", Instruments_input[,"file_name"])
# Instruments_input[,"pull"] <- sapply(gregexpr("\\/", Instruments_input[,"file_name"]), tail, 1)
# Instruments_input[,"file_name"] <- substr(Instruments_input[,"file_name"],Instruments_input[,"pull"]+1,nchar(Instruments_input[,"file_name"]))
# 
# Instruments_input[,"pull"] <- Instruments_input[,"file_name"]
# Instruments_input[,"pull"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=Instruments_input[,"pull"])

#Instruments_pull <- c("04-2014","02-2013","02-2012")
#Instruments_yr <- c(2014,2013,2012)

#Instruments_files <- c("EurekahedgeHF_EXCEL_2014Apr_Instruments_Traded","EurekahedgeHF_EXCEL_2013Feb_Instruments_Traded",
#                       "EurekahedgeHF_EXCEL_2012Feb_Instruments_Traded")

#Instruments_input <- data.frame(matrix(NA, ncol=3, nrow=length(Instruments_files), dimnames=list(c(), c("pull","yr","file"))), 
#                                stringsAsFactors=FALSE)

#Instruments_input[,"pull"] <- Instruments_pull
#Instruments_input[,"yr"] <- Instruments_yr
#Instruments_input[,"file"] <- Instruments_files

Instruments_concatenate0 <- alply(.data=Instruments_input, .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- Instruments_input[1,]
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
  
  #assign(x[,"file"], input, envir = .GlobalEnv)
  
  input  <- input[order(input[,"Fund_ID"],input[,"Fund_Name"],input[,"Instrument_Traded"]),]
  row.names(input) <- seq(nrow(input))
  
  gc()
  
  return(input)
  
}, directory_in=output_directory, unknowns=unknowns_strings, .expand = TRUE, .progress = "text")

Instruments_concatenate <- rbind.fill(Instruments_concatenate0)

#rm2(Instruments_pull,Instruments_yr,Instruments_files)
rm2(Instruments_concatenate0)

for(i in which(sapply(Instruments_concatenate,class)=="character"))
{
  Instruments_concatenate[[i]] = trim(Instruments_concatenate[[i]])
}
rm2(i)
for (i in 1:ncol(Instruments_concatenate))
{
  Instruments_concatenate[,i] <- unknownToNA(Instruments_concatenate[,i], unknown=unknowns_strings,force=TRUE)
  Instruments_concatenate[,i] <- ifelse(is.na(Instruments_concatenate[,i]),NA, Instruments_concatenate[,i])
} 
rm2(i)

Instruments_concatenate  <- Instruments_concatenate[order(Instruments_concatenate[,"Fund_ID"],
                                                          Instruments_concatenate[,"Fund_Name"],
                                                          Instruments_concatenate[,"Instrument_Traded"]),]
row.names(Instruments_concatenate) <- seq(nrow(Instruments_concatenate))

rm2(Instruments_input)


###############################################################################
cat("SECTION: CLEAN INSTRUMENTS", "\n")
###############################################################################

Instruments_clean <- Instruments_concatenate

Instruments_clean[,"Instrument_Traded"] <-  gsub(pattern=" ", replacement=".", x=Instruments_clean[,"Instrument_Traded"])
Instruments_clean[,"Instrument_Traded"] <-  gsub(pattern="\\.{2,}", replacement="\\.", x=Instruments_clean[,"Instrument_Traded"])
Instruments_clean[,"Instrument_Traded"] <-  gsub(pattern="-", replacement="_", x=Instruments_clean[,"Instrument_Traded"])

Instruments_u <-  unique(Instruments_clean[!is.na(Instruments_clean[,"Instrument_Traded"]),"Instrument_Traded"])

rm2(Instruments_concatenate)


###############################################################################
cat("SECTION: EXPAND IDS", "\n")
###############################################################################

Instruments_names <- data.frame(unique(Instruments_clean[,c("pull","Fund_ID","Fund_Name")]),stringsAsFactors=FALSE)

Instruments_ids <- data.frame(unique(Instruments_clean[,c("pull","Fund_ID")]),temp_col=NA,stringsAsFactors=FALSE)
colnames(Instruments_ids)[match("temp_col",names(Instruments_ids))] <- "Instrument_Traded_type"

Instruments_ids_expand <- coredata(Instruments_ids)[rep(seq(nrow(Instruments_ids)),length(Instruments_u)),]
Instruments_ids_expand  <- Instruments_ids_expand[order(Instruments_ids_expand[,"Fund_ID"],Instruments_ids_expand[,"pull"]),]
row.names(Instruments_ids_expand) <- seq(nrow(Instruments_ids_expand))
Instruments_ids_expand[,"Instrument_Traded_type"] <- rep(Instruments_u,nrow(Instruments_ids))

rm2(Instruments_ids,Instruments_u)

Instruments_expand <- Instruments_clean[!is.na(Instruments_clean[,"Instrument_Traded"]),]
colnames(Instruments_expand)[match("Instrument_Traded",names(Instruments_expand))] <- "Instrument_Traded_type"
colnames(Instruments_expand)[match("Fund_Name",names(Instruments_expand))] <- "Instrument_Traded"
Instruments_expand[,"Instrument_Traded"] <- 1

rm2(Instruments_clean)

Instruments_ids_merge <- merge(Instruments_ids_expand, Instruments_expand, 
                               by.x=c("pull","Fund_ID","Instrument_Traded_type"), 
                               by.y=c("pull","Fund_ID","Instrument_Traded_type"), 
                               all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

rm2(Instruments_ids_expand,Instruments_expand)

Instruments_ids_merge  <- Instruments_ids_merge[order(Instruments_ids_merge[,"pull"],
                                                      Instruments_ids_merge[,"Fund_ID"],
                                                      Instruments_ids_merge[,"Instrument_Traded_type"]),]
row.names(Instruments_ids_merge) <- seq(nrow(Instruments_ids_merge))

Instruments_ids_merge[,"Instrument_Traded"] <- ifelse(is.na(Instruments_ids_merge[,"Instrument_Traded"]), 0, Instruments_ids_merge[,"Instrument_Traded"])

Instruments_cast <- reshape(Instruments_ids_merge, direction = "wide", idvar=c("pull", "Fund_ID"), timevar="Instrument_Traded_type")

Instruments_cast  <- Instruments_cast[order(Instruments_cast[,"Fund_ID"],Instruments_cast[,"pull"]),]
row.names(Instruments_cast) <- seq(nrow(Instruments_cast))

rm2(Instruments_ids_merge)


###############################################################################
cat("SECTION: MERGE NAMES IDS", "\n")
###############################################################################

Instruments_final <-  merge(Instruments_names, Instruments_cast, 
                            by.x=c("pull","Fund_ID"), 
                            by.y=c("pull","Fund_ID"), 
                            all.x=TRUE, all.y=FALSE, sort=FALSE,suffixes=c(".x",".y"))

Instruments_final  <- Instruments_final[order(Instruments_final[,"Fund_ID"],
                                              Instruments_final[,"pull"]),]
row.names(Instruments_final) <- seq(nrow(Instruments_final))

rm2(Instruments_names,Instruments_cast)


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

#Check to see if final folder exists.  If not, create it.
final_folder_path <- paste(output_directory, "Final", sep = "//", collapse = "//")  
create_directory(final_folder_path,remove=1)

write.csv(Instruments_final, file=paste(final_folder_path,"//","EurekahedgeHF_Instruments_Traded",".csv",sep=""),row.names=FALSE)

rm2(Instruments_final,final_folder_path)

