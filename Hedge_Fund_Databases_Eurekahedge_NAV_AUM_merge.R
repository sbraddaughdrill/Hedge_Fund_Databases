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

#source(file=paste(function_directory,"functions_db.R",sep="\\"),echo=FALSE)
source(file=paste(function_directory,"functions_statistics.R",sep="\\"),echo=FALSE)
#source(file=paste(function_directory,"functions_text_analysis.R",sep="\\"),echo=FALSE)
#source(file=paste(function_directory,"functions_text_parse.R",sep="\\"),echo=FALSE)
source(file=paste(function_directory,"functions_utilities.R",sep="\\"),echo=FALSE)


###############################################################################
# LIBRARIES;
cat("SECTION: LIBRARIES", "\n")
###############################################################################

#Load External Packages

external_packages <- c("compare","cwhmisc","data.table","DataCombine","fastmatch","foreign","formatR","gdata","gtools",
                       "Hmisc","lubridate","pbapply","plyr","reshape2",
                       "stringr","zoo")
invisible(unlist(sapply(external_packages,load_external_packages, repo_str=repo, simplify=FALSE, USE.NAMES=FALSE)))
installed_packages <- list_installed_packages(external_packages)

rm(external_packages,installed_packages,repo)


###############################################################################
cat("SECTION: IMPORT NAV & AUM", "\n")
###############################################################################

NAV_AUM_input <- data.frame(read.csv(file=paste(output_directory,"\\","EurekahedgeHF_NAV_AUM_files",".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),stringsAsFactors=FALSE)


###############################################################################
cat("SECTION: CONCATENATE NAV & AUM", "\n")
###############################################################################

melt_folder_path <- paste(output_directory, "NAV_AUM_melt", sep = "//", collapse = "//")  

NAV_AUM_id_cols <- c("RetAUM","Fund_ID","Fund_Name")

NAV_AUM_concatentate0 <- alply(.data=NAV_AUM_input, .margins=1, .fun = function(x,directory_in){
  
  # x <- NAV_AUM_input[7,]
  # directory_in <- melt_folder_path
  
  input <- data.frame(read.csv(file=paste(directory_in,"\\",x[,"pull"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      #read.csv(file=paste(x[,"file_clean"],sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE),
                      #yr=NA,month=NA,
                      stringsAsFactors=FALSE)
  #bad <- input[is.na(input[,"Fund.ID"]),]
  
  return(input)
  
}, directory_in=melt_folder_path, .expand = TRUE, .progress = "text")

NAV_AUM_concatentate <- rbind.fill(NAV_AUM_concatentate0)

#rm2(NAV_AUM_pull,NAV_AUM_yr,NAV_AUM_files)
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

NAV_AUM_concatentate <- NAV_AUM_concatentate[,c("pull","Fund_ID","Fund_Name","date","Monthly_Ret","AUM")]

#NAV_AUM_concatentate <- NAV_AUM_concatentate[rowSums(is.na(NAV_AUM_concatentate[,1:ncol(NAV_AUM_concatentate)]))<ncol(NAV_AUM_concatentate),]

NAV_AUM_concatentate[,"date"] <- as.yearmon(NAV_AUM_concatentate[,"date"],format="%b %Y")
NAV_AUM_concatentate[,"date"] <- as.Date(NAV_AUM_concatentate[,"date"],format="%m-%d-%Y")

#NAV_AUM_concatentate[,"yr"] <- year(NAV_AUM_concatentate[,"date"])
#NAV_AUM_concatentate[,"month"] <- month(NAV_AUM_concatentate[,"date"])

invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: FIND PULL-FUND COMBOS", "\n")
###############################################################################

NAV_AUM_concatentate <- NAV_AUM_concatentate[order(NAV_AUM_concatentate[,"Fund_ID"],NAV_AUM_concatentate[,"Fund_Name"],
                                                   NAV_AUM_concatentate[,"date"], NAV_AUM_concatentate[,"pull"]),]

row.names(NAV_AUM_concatentate) <- seq(nrow(NAV_AUM_concatentate))

invisible(gc(verbose = FALSE, reset = TRUE))

#aa <- NAV_AUM_concatentate[1:10000,]

# NAV_AUM_dates <- ddply(.data=NAV_AUM_concatentate[,(colnames(NAV_AUM_concatentate) %in% c("pull","date"))], .variables=c("pull"), .fun = function(x){
#   
#   return(data.frame(min_date=min(x[,"date"]),max_date=max(x[,"date"]),stringsAsFactors=FALSE)) 
# },.progress = "text")

NAV_AUM_dates_dt <- data.table(NAV_AUM_concatentate[,(colnames(NAV_AUM_concatentate) %in% c("pull","date"))])
setkeyv(NAV_AUM_dates_dt,"pull")
NAV_AUM_dates <- NAV_AUM_dates_dt[,j=list(min_date=min(date),max_date=max(date)), by = c("pull")]

rm2(NAV_AUM_dates_dt)

NAV_AUM_dates <- as.data.frame(NAV_AUM_dates,stringsAsFactors=FALSE)

NAV_AUM_dates <- NAV_AUM_dates[order(NAV_AUM_dates[,"pull"]),]
row.names(NAV_AUM_dates) <- seq(nrow(NAV_AUM_dates))

invisible(gc(verbose = FALSE, reset = TRUE))

bb_pull_id0 <- data.table(NAV_AUM_concatentate[,(colnames(NAV_AUM_concatentate) %in% c("pull","date"))])
setkeyv(bb_pull_id0,NULL)
bb_pull_id0 <- as.data.frame(unique(bb_pull_id0),pull_trim=NA,stringsAsFactors=FALSE)

bb_pull_id0[,"pull_trim"] <- bb_pull_id0[,"pull"]
bb_pull_id0[,"pull_trim"] <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="",bb_pull_id0[,"pull_trim"])

# bb <- data.table(bb_pull_id0[,(colnames(bb_pull_id0) %in% c("pull_trim","Fund_ID","date","yr","month"))])
# setkeyv(bb, c("Fund_ID","date"))
# bb <- unique(bb,by=c("Fund_ID","date"),fromFirst=TRUE)

bb_pull_id1 <- data.table(bb_pull_id0[,!(colnames(bb_pull_id0) %in% c("pull"))])
setkeyv(bb_pull_id1, c("date"))
bb_pull_id1 <- unique(bb_pull_id1,by=c("date"),fromFirst=TRUE)
bb_pull_id1 <- as.data.frame(bb_pull_id1,stringsAsFactors=FALSE)

bb_pull_id2 <- merge(bb_pull_id1, unique(bb_pull_id0[,(colnames(bb_pull_id0) %in% c("pull_trim","pull"))]), 
                     by.x=c("pull_trim"), by.y=c("pull_trim"), 
                     all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))

rm2(bb_pull_id0,bb_pull_id1)

bb_pull_id2 <- bb_pull_id2[order(bb_pull_id2[,"date"]),]
row.names(bb_pull_id2) <- seq(nrow(bb_pull_id2))

bb_pull_id3 <- merge(bb_pull_id2, NAV_AUM_dates, 
                     by.x=c("pull"), by.y=c("pull"), 
                     all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))

rm2(bb_pull_id2,NAV_AUM_dates)

bb_pull_id3 <- bb_pull_id3[order(bb_pull_id3[,"date"]),]
row.names(bb_pull_id3) <- seq(nrow(bb_pull_id3))


bb_pull_id <- data.frame(bb_pull_id3,yr=NA,month=NA,bad_min=NA,bad_max=NA,stringsAsFactors=FALSE)

rm2(bb_pull_id3)

bb_pull_id[,"yr"] <- year(bb_pull_id[,"date"])
bb_pull_id[,"month"] <- month(bb_pull_id[,"date"])
bb_pull_id[,"bad_min"] <- ifelse(bb_pull_id[,"date"]<bb_pull_id[,"min_date"],1,0)
bb_pull_id[,"bad_max"] <- ifelse(bb_pull_id[,"date"]>bb_pull_id[,"max_date"],1,0)

bb_pull_id <- bb_pull_id[,c("pull_trim","pull","date","yr","month","min_date","max_date","bad_min","bad_max")]

invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: REMOVE DUPLICATES", "\n")
###############################################################################

#aa <- NAV_AUM_concatentate[1:10000,]

bb_seq_dt <- data.table(NAV_AUM_concatentate[,!(colnames(NAV_AUM_concatentate) %in% c("Fund_Name"))])

#rm2(NAV_AUM_concatentate)

bb_seq <- bb_seq_dt[, Order := sequence(.N), by = c("Fund_ID","date")]

rm2(bb_seq_dt)

#bb_seq <- as.data.frame(bb_seq,stringsAsFactors=FALSE)
#bb_seq <- bb_seq[order(bb_seq[,"Fund_ID"],bb_seq[,"date"],bb_seq[,"date"]),]
#row.names(bb_seq) <- seq(nrow(bb_seq))

#bb_seq <- ddply(.data=NAV_AUM_concatentate[,!(colnames(NAV_AUM_concatentate) %in% c("Fund_Name","yr","month"))],.variables=c("Fund_ID","date"), .fun = function(y){
#  return(data.frame(y,Order=seq(1,nrow(y)),stringsAsFactors=FALSE))
#  #return(data.frame(y,Order=seq(1,nrow(y)),last=c(rep.int(0, nrow(y)-1),1),stringsAsFactors=FALSE))
#},.progress = "text")

### Ret
# 
# #bb_ret_freq <- count(bb_seq[,(colnames(bb_seq) %in% c("Fund_ID","date","yr","month","Monthly_Ret"))],
# #                     c("Fund_ID","date","yr","month","Monthly_Ret"))
# 
# bb_ret_freq_dt <- bb_seq
# bb_ret_freq_dt[,c("AUM", "Order")] <- NULL
# bb_ret_freq <- bb_ret_freq_dt[, freq := .N, by = c("Fund_ID","date","Monthly_Ret")]
# bb_ret_freq[,c("pull")] <- NULL
# bb_ret_freq <- unique(bb_ret_freq)
# 
# rm2(bb_ret_freq_dt)
# 
# # bb_ret_order <- ddply(.data=bb_seq[,(colnames(bb_seq) %in% c("Fund_ID","date","yr","month","Monthly_Ret","Order"))],.variables=c("Fund_ID","date","yr","month","Monthly_Ret"), .fun = function(y){
# #   
# #   return(data.frame(avg_order=mean(y[,"Order"]),stringsAsFactors=FALSE))
# # },.progress = "text")
# 
# bb_ret_order_dt <- bb_seq
# bb_ret_order_dt[,c("pull","AUM")] <- NULL
# bb_ret_order <- bb_ret_order_dt[, avg_order := mean(Order), by = c("Fund_ID","date","Monthly_Ret")]
# bb_ret_order[,c("Order")] <- NULL
# bb_ret_order <- unique(bb_ret_order)
# 
# rm2(bb_ret_order_dt)
# 
# bb_ret_full <- merge(as.data.frame(bb_ret_freq,stringsAsFactors=FALSE),
#                      as.data.frame(bb_ret_order,stringsAsFactors=FALSE), 
#                      by.x=c("Fund_ID","date","Monthly_Ret"),
#                      by.y=c("Fund_ID","date","Monthly_Ret"), 
#                      all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))
# 
# bb_ret_full <- bb_ret_full[order(bb_ret_full[,"Fund_ID"],bb_ret_full[,"date"],-bb_ret_full[,"freq"],bb_ret_full[,"avg_order"]),]
# row.names(bb_ret_full) <- seq(nrow(bb_ret_full))

bb_ret_full_dt <- bb_seq
bb_ret_full_dt[,c("pull","AUM")] <- NULL
bb_ret_full <- bb_ret_full_dt[,j=list(freq=.N,avg_order=mean(Order)), by = c("Fund_ID","date","Monthly_Ret")]

rm2(bb_ret_full_dt)

bb_ret_full <- as.data.frame(bb_ret_full,stringsAsFactors=FALSE)

bb_ret_full <- bb_ret_full[order(bb_ret_full[,"Fund_ID"],bb_ret_full[,"date"],-bb_ret_full[,"freq"],bb_ret_full[,"avg_order"]),]
row.names(bb_ret_full) <- seq(nrow(bb_ret_full))

invisible(gc(verbose = FALSE, reset = TRUE))

# bb_ret_full_trim <- ddply(.data=bb_ret_full,.variables=c("Fund_ID","date","yr","month"), .fun = function(y){
#   
#   return(head(y,1))
# },.progress = "text")

bb_ret_full <- data.table(bb_ret_full)
setkeyv(bb_ret_full, c("Fund_ID","date"))
bb_ret_full_trim <- unique(bb_ret_full,by=c("Fund_ID","date"),fromFirst=TRUE)

rm2(bb_ret_full)

bb_ret_full_trim <- as.data.frame(bb_ret_full_trim,stringsAsFactors=FALSE)

invisible(gc(verbose = FALSE, reset = TRUE))


### AUM

bb_aum_full_dt <- bb_seq
bb_aum_full_dt[,c("pull","Monthly_Ret")] <- NULL
bb_aum_full <- bb_aum_full_dt[,j=list(freq=.N,avg_order=mean(Order)), by = c("Fund_ID","date","AUM")]

rm2(bb_aum_full_dt)

bb_aum_full <- as.data.frame(bb_aum_full,stringsAsFactors=FALSE)

bb_aum_full <- bb_aum_full[order(bb_aum_full[,"Fund_ID"],bb_aum_full[,"date"],-bb_aum_full[,"freq"],bb_aum_full[,"avg_order"]),]
row.names(bb_aum_full) <- seq(nrow(bb_aum_full))

invisible(gc(verbose = FALSE, reset = TRUE))

# bb_aum_full_trim <- ddply(.data=bb_aum_full,.variables=c("Fund_ID","date","yr","month"), .fun = function(y){
#   
#   aumurn(head(y,1))
# },.progress = "text")

bb_aum_full <- data.table(bb_aum_full)
setkeyv(bb_aum_full, c("Fund_ID","date"))
bb_aum_full_trim <- unique(bb_aum_full,by=c("Fund_ID","date"),fromFirst=TRUE)

rm2(bb_aum_full)

bb_aum_full_trim <- as.data.frame(bb_aum_full_trim,stringsAsFactors=FALSE)

invisible(gc(verbose = FALSE, reset = TRUE))

### Merge

rm2(bb_seq)

bb_ret_aum_full <- merge(bb_ret_full_trim[,!(colnames(bb_ret_full_trim) %in% c("freq","avg_order"))],
                         bb_aum_full_trim[,!(colnames(bb_aum_full_trim) %in% c("freq","avg_order"))], 
                         by.x=c("Fund_ID","date"),
                         by.y=c("Fund_ID","date"), 
                         all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))

rm2(bb_ret_full_trim,bb_aum_full_trim)
invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: OUTPUT DATA", "\n")
###############################################################################

NAV_AUM_collapse <- merge(bb_ret_aum_full,bb_pull_id, 
                         by.x=c("date"),by.y=c("date"), 
                         all.x=TRUE, all.y=FALSE, sort=FALSE, suffixes=c(".x",".y"))

rm2(bb_ret_aum_full,bb_pull_id)
invisible(gc(verbose = FALSE, reset = TRUE))

NAV_AUM_collapse <- NAV_AUM_collapse[order(NAV_AUM_collapse[,"Fund_ID"],
                                           NAV_AUM_collapse[,"date"],
                                           NAV_AUM_collapse[,"pull_trim"],
                                           NAV_AUM_collapse[,"pull"]),]
row.names(NAV_AUM_collapse) <- seq(nrow(NAV_AUM_collapse))

invisible(gc(verbose = FALSE, reset = TRUE))

NAV_AUM_collapse <- NAV_AUM_collapse[,c("pull_trim","pull","Fund_ID","Monthly_Ret","AUM",
                                        "date","yr","month","min_date","max_date","bad_min","bad_max")]

#Check to see if nav_aum_ret folder exists.  If not, create it.
nav_aum_ret_folder_path <- paste(output_directory, "NAV_AUM_Ret", sep = "//", collapse = "//")  
create_directory(nav_aum_ret_folder_path,remove=1)

write.csv(NAV_AUM_collapse, file=paste(nav_aum_ret_folder_path,"//","EurekahedgeHF_NAV_AUM_merge",".csv",sep=""),row.names=FALSE)

rm2(NAV_AUM_collapse,nav_aum_ret_folder_path)
