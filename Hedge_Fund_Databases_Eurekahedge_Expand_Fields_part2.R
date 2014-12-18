# TODO: Add comment
# 
# Author:  Brad
# File:    Hedge_Fund_Databases_Eurekahedge_Expand_Fields_part2.R
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

final_folder_expand1_files0 <- data.frame(files=list.files(path=final_folder_expand1_path),file_name=NA,import=NA,stringsAsFactors=FALSE)

#final_folder_expand1_files <- final_folder_expand1_files0[!grepl(".TXT|.txt", final_folder_expand1_files0[,"files"]),]
final_folder_expand1_files <- final_folder_expand1_files0[grepl(".CSV|.csv", final_folder_expand1_files0[,"files"]),]

rm2(final_folder_expand1_files0)

final_folder_expand1_files[,"file_name"] <- final_folder_expand1_files[,"files"]
final_folder_expand1_files[,"file_name"] <- gsub(pattern="(.CSV|.csv)", replacement="", x=final_folder_expand1_files[,"file_name"])

final_folder_expand1_files[,"import"] <- ifelse(grepl("(Stats|Fund_Detail|Fee_and_Redemption|Profile_Strategy|Identifier)",final_folder_expand1_files[,"file_name"]),1,final_folder_expand1_files[,"import"])
final_folder_expand1_files[,"import"] <- ifelse(grepl("(NAV_AUM_Ret|Instruments_Traded)",final_folder_expand1_files[,"file_name"]),2,final_folder_expand1_files[,"import"])
final_folder_expand1_files[,"import"] <- ifelse(grepl("(Other)",final_folder_expand1_files[,"file_name"]),3,final_folder_expand1_files[,"import"])
final_folder_expand1_files[,"import"] <- ifelse(is.na(final_folder_expand1_files[,"import"]),0,final_folder_expand1_files[,"import"])
final_folder_expand1_files[,"import"] <- ifelse(grepl("(_part2)",final_folder_expand1_files[,"file_name"]),0,final_folder_expand1_files[,"import"])

invisible(gc(verbose = FALSE, reset = TRUE))


###############################################################################
cat("SECTION: IMPORT FILES", "\n")
###############################################################################

a_ply(.data=final_folder_expand1_files[final_folder_expand1_files[,"import"] %in% c(1),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_expand1_files[1,]
  # x <- final_folder_expand1_files[3,]
  # x <- final_folder_expand1_files[4,] 
  
  # directory_in <- final_folder_expand2_path
  # unknowns <- unknowns_strings
  
  #temp_cols <- c("date","yr","month","bad_min","bad_max")
  
  input <- data.table(read.csv(file=paste(final_folder_expand1_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE))

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
  set(input,i=NULL,j=which(colnames(input)==c("droprow")),value=NULL)
  #input <- input[(droprow)][,droprow:=NULL][]
  
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

a_ply(.data=final_folder_expand1_files[final_folder_expand1_files[,"import"] %in% c(3),], .margins=1, .fun = function(x,directory_in,unknowns){
  
  # x <- final_folder_expand1_files[6,]
  
  # directory_in <- final_folder_expand2_path
  # unknowns <- unknowns_strings
  
  input <- data.table(read.csv(file=paste(final_folder_expand1_path,"//",x[,"file_name"],".csv",sep=""),header=TRUE,na.strings="NA",stringsAsFactors=FALSE))
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
  set(input,i=NULL,j=which(colnames(input)==c("droprow")),value=NULL)
  #input <- input[(droprow)][,droprow:=NULL][]
  
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


###############################################################################
cat("SECTION: SPLIT (EXPAND) COLUMNS", "\n")
###############################################################################

# EurekahedgeHF_Stats_noreturns_split_expand <- EurekahedgeHF_Stats_noreturns_yn_binary
# EurekahedgeHF_Fund_Details_split_expand <- EurekahedgeHF_Fund_Details_yn_binary
# EurekahedgeHF_Fee_and_Redemption_split_expand <- EurekahedgeHF_Fee_and_Redemption_yn_binary
# EurekahedgeHF_Profile_Strategy_split_expand <- EurekahedgeHF_Profile_Strategy_yn_binary
# EurekahedgeHF_Identifiers_split_expand <- EurekahedgeHF_Identifiers_yn_binary
# EurekahedgeHF_Other_split_expand <- EurekahedgeHF_Other_yn_binary
# #EurekahedgeHF_Instruments_Traded_split_expand <- EurekahedgeHF_Instruments_Traded_yn_binary
# 
# rm2(EurekahedgeHF_Stats_noreturns_yn_binary,EurekahedgeHF_Fund_Details_yn_binary,EurekahedgeHF_Fee_and_Redemption_yn_binary)
# rm2(EurekahedgeHF_Profile_Strategy_yn_binary,EurekahedgeHF_Identifiers_yn_binary,EurekahedgeHF_Other_yn_binary)
# #rm2(EurekahedgeHF_Instruments_Traded_split_expand)
# 
# split_expand_temp1 <- list(data=c("EurekahedgeHF_Other_split_expand"),
#                            col=c("Custodian","Principal_Prime_Broker_combcol","Secondary_Prime_Broker","Synthetic_Prime_Broker",
#                                  "Legal_Advisor_Offshore","Legal_Advisor_Onshore","Legal_Advisor"))
# 
# split_expand_all <- list(split_expand_temp1)
# 
# rm2(split_expand_temp1)
# invisible(gc(verbose = FALSE, reset = TRUE))


split_expand_temp1 <- list(data_in=c("EurekahedgeHF_Stats_noreturns_part1"),data_out=c("EurekahedgeHF_Stats_noreturns_part2"),
                           col=NULL)
split_expand_temp2 <- list(data_in=c("EurekahedgeHF_Fund_Details_part1"),data_out=c("EurekahedgeHF_Fund_Details_part2"),
                           col=NULL)
split_expand_temp3 <- list(data_in=c("EurekahedgeHF_Fee_and_Redemption_part1"),data_out=c("EurekahedgeHF_Fee_and_Redemption_part2"),
                           col=NULL)
split_expand_temp4 <- list(data_in=c("EurekahedgeHF_Profile_Strategy_part1"),data_out=c("EurekahedgeHF_Profile_Strategy_part2"),
                           col=NULL)
split_expand_temp5 <- list(data_in=c("EurekahedgeHF_Identifiers_part1"),data_out=c("EurekahedgeHF_Identifiers_part2"),
                           col=NULL)
split_expand_temp6 <- list(data_in=c("EurekahedgeHF_Other_part1"),data_out=c("EurekahedgeHF_Other_part2"),
                           col=c("Custodian","Principal_Prime_Broker_combcol","Secondary_Prime_Broker","Synthetic_Prime_Broker",
                                 "Legal_Advisor_Offshore","Legal_Advisor_Onshore","Legal_Advisor"))

split_expand_all0 <- list(split_expand_temp1,split_expand_temp2,split_expand_temp3,
                          split_expand_temp4,split_expand_temp5,split_expand_temp6)

rm2(split_expand_temp1,split_expand_temp2,split_expand_temp3,split_expand_temp4,split_expand_temp5,split_expand_temp6)
invisible(gc(verbose = FALSE, reset = TRUE))

split_expand_all1 <- llply(.data=split_expand_all0, .fun = function(x){
  
  return(xout <- c(x,ncol=ncol(get(x[[which(names(x)==c("data_in"))]])),
                   nrow=nrow(get(x[[which(names(x)==c("data_in"))]])),size=object.size(get(x[[which(names(x)==c("data_in"))]]))))
  
}, .progress = "text")

rm2(split_expand_all0)

split_expand_all1 <- split_expand_all1[order(-sapply(split_expand_all1,"[[","size"),
                                             -sapply(split_expand_all1,"[[","ncol"),
                                             -sapply(split_expand_all1,"[[","nrow"))]

invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))
invisible(gc(verbose = FALSE, reset = TRUE))

l_ply(.data=split_expand_all1, .fun = function(x){
  
  # x <- split_expand_all1[[1]]
  # x <- split_expand_all1[[2]]
  
  if (length(x[[which(names(x)==c("col"))]])==0) {
    
    #cat("NO", "\n")
    
    #assign(x[[which(names(x)==c("data_out"))]], get(x[[which(names(x)==c("data_in"))]]), envir = .GlobalEnv)
    write.csv(get(x[[which(names(x)==c("data_in"))]]), file=paste(final_folder_expand2_path,"//",x[[which(names(x)==c("data_out"))]],".csv",sep=""),row.names=FALSE)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
  } else {
    
    #cat("YES", "\n")
    
    data_cols <- ldply(.data=x[[which(names(x)==c("col"))]], .fun = function(y,data_temp){
      
      # y <- x[[which(names(x)==c("col"))]][1]
      
      countCharOccurrences <- function(char, s) {
        s2 <- gsub(char,"",s)
        return (nchar(s) - nchar(s2))
      }
      
      return(data.frame(order_org=NA,order_size=NA,col=y,max=max(countCharOccurrences(",",data_temp[[which(colnames(data_temp) %in% y)]])),stringsAsFactors=FALSE))
      
    }, data_temp=data.table(get(x[[which(names(x)==c("data_in"))]])), .progress = "none")

    data_cols[,"order_org"] <- seq(1,nrow(data_cols))
    data_cols  <- data_cols[order(-data_cols[,"max"],data_cols[,"col"]),]
    row.names(data_cols) <- seq(nrow(data_cols))
    #colnames(data_cols)[match("max",names(data_cols))] <- "order_size"
    data_cols[,"order_size"] <- seq(1,nrow(data_cols))
    
    invisible(gc(verbose = FALSE, reset = TRUE))
    
    data_count0 <- dlply(.data=data_cols,.variables=c("order_size"), .fun = function(y,data_temp){
      
      # y <- data_cols[data_cols[,"col"]=="Custodian",]
      y_trim <- y[,"col"]
      
      for (k in which(sapply(data_temp,class)=="character")) 
      {
        set(data_temp, i=NULL, j=k, value=gsub(" ,", ",", data_temp[[k]], perl=TRUE))
        set(data_temp, i=NULL, j=k, value=gsub(", ", ",", data_temp[[k]], perl=TRUE))
        set(data_temp, i=NULL, j=k, value=gsub(",+", ",", data_temp[[k]], perl=TRUE))
        set(data_temp, i=NULL, j=k, value=gsub(",", ", ", data_temp[[k]], perl=TRUE))
      }
      rm(k)
      
      temp4 <- data.frame(stri_split_fixed(data_temp[[which(colnames(data_temp) %in% y_trim)]], ",", simplify = TRUE),temp_count=NA,
                          stringsAsFactors=FALSE)
      
      temp4_expand_cols <- colnames(temp4)[!(colnames(temp4) %in% "temp_count")]
      temp4_expand_cols <- gsub(pattern="([[:alpha:]]|[[:punct:]])", replacement="", temp4_expand_cols, perl=TRUE)
      colnames(temp4) <- c(paste(y_trim,sprintf("%02d", seq(1,length(temp4_expand_cols))),sep=""),"temp_count")
      
      #temp4 <- as.data.table(stri_split_fixed(data_temp[[which(colnames(data_temp) %in% y_trim)]], ",", simplify = TRUE))
      #setnames(temp4,colnames(temp4),paste(y_trim,sprintf("%02d", seq(1,ncol(temp4))),sep=""))
      
      rm(temp4_expand_cols)
      invisible(gc(verbose = FALSE, reset = TRUE))
      
      temp4 <- as.data.table(temp4)
      
      cols <- colnames(temp4)
      for (l in cols) 
      {
        #l <- 1
        set(temp4, i=NULL, j=l, value=gsub(pattern=" {2,}", replacement=" ", temp4[[l]], perl=TRUE))
        set(temp4, i=NULL, j=l, value=gsub(pattern="^\\s+|\\s+$", replacement="", temp4[[l]], perl=TRUE))
      }
      rm(l,cols)
      
      temp4[,temp_count := rowSums(!is.na(.SD)),.SDcols=colnames(temp4)]
      
      invisible(gc(verbose = FALSE, reset = TRUE))
      
      cols <- "temp_count"
      for (l in cols) 
      {
        #l <- 1
        set(temp4, i=NULL, j=l, value=ifelse(temp4[[l]]==0,NA,temp4[[l]]))
      }
      rm(l,cols)
      
      invisible(gc(verbose = FALSE, reset = TRUE))
      
      #temp4[,temp_count:=NA_integer_] 
      
      setnames(temp4,"temp_count",paste(y_trim,"count",sep="_"))
      
      rm(y_trim)
      invisible(gc(verbose = FALSE, reset = TRUE))
      
      return(temp4)
      
    }, data_temp=data.table(get(x[[which(names(x)==c("data_in"))]])), .progress = "text")
    #}, data_temp=unique(set(data.table(get(x[[which(names(x)==c("data_in"))]])), j=which(!(colnames(data.table(get(x[[which(names(x)==c("data_in"))]]))) %in% c("Fund_ID",col_temp))), value=NULL)), .progress = "text") 
    
    rm(data_cols)
    invisible(gc(verbose = FALSE, reset = TRUE))
    
    df_all_outc <- data.frame(get(x[[which(names(x)==c("data_in"))]]),
                              matrix(NA, ncol=length(unlist(lapply(data_count0,colnames))), dimnames=list(c(), unlist(lapply(data_count0,colnames)))),
                              stringsAsFactors=FALSE)
    
    rm(list = x[[which(names(x)==c("data_in"))]],envir = .GlobalEnv)
    invisible(gc(verbose = FALSE, reset = TRUE))
    
    for (k in 1:length(data_count0)) 
    {
      for (l in 1:ncol(data_count0[[k]])) 
      {
        if(colnames(data_count0[[k]])[l] %in% colnames(df_all_outc)){set(df_all_outc, i=NULL, j=which(colnames(df_all_outc)==colnames(data_count0[[k]])[l]), value=data_count0[[k]][[l]])}
      
        invisible(gc(verbose = FALSE, reset = TRUE))
        
        progress_function(outer_loop_count=k, outer_loop_start_val=1, outer_loop_end_val=length(data_count0),inner_loop_count=l, inner_loop_start_val=1, inner_loop_end_val=ncol(data_count0[[k]]))
      }
      rm(l)
    }
    rm(k)
    
    rm(data_count0)
    invisible(gc(verbose = FALSE, reset = TRUE))
    
    #assign(x[[which(names(x)==c("data_out"))]], data.table(df_all_outc), envir = .GlobalEnv)
    write.csv(df_all_outc, file=paste(final_folder_expand2_path,"//",x[[which(names(x)==c("data_out"))]],".csv",sep=""),row.names=FALSE)
    
    rm(df_all_outc)

    invisible(gc(verbose = FALSE, reset = TRUE))
    
  }
  
}, .progress = "text")

rm2(split_expand_all1)


