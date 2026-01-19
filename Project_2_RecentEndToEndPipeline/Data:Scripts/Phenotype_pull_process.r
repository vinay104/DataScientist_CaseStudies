# Load required libraries 
library(dplyr)
library(arrow)
library(tidyverse)
library(paws)

aws_account = 'xxxx'
bucket = paste0("xxxx", aws_account, "xxxxx")
output_path = "xxxx/"

# AWS connection and pull the phenotype data 
Sys.setenv(
  AWS_SHARED_CREDENTIALS_FILE = "xxxx", ## path to AWS credentials
  AWS_PROFILE = paste0("xxxxx", aws_account),
  AWS_CONFIG = "xxxxx"
)

s3_client = paws::s3() 

BLUEs_dir = "xxxxx/"

objects = s3_client$list_objects_v2(
  Bucket = bucket,
  Prefix = BLUEs_dir
)

file_list = sapply(objects$Contents, function(x) x$Key)


blues_files = file_list[str_detect(file_list, "xxxx\\.csv$")]
blues_files90 = file_list[str_detect(file_list, "xxxx\\.csv$")]
blues_files120 = file_list[str_detect(file_list, "xxxx\\.csv$")]

all_blues_files <- c(blues_files, blues_files90, blues_files120)

blues_data = map_dfr(blues_files, function(key) {
  obj = s3_client$get_object(Bucket = bucket, Key = key)
  read.csv(text = rawToChar(obj$Body), stringsAsFactors = FALSE)
})

write.csv(blues_data, file = file.path(output_path, "xxxx.csv"), row.names = FALSE)

# Phenotype data processing and sub-setting relevant training set data 
Pheno = read.csv("xxxxx.csv", check.names = FALSE) %>%  
  mutate(
    trial = as.factor(trial),
    female = as.factor(Female.name),
    male   = as.factor(Male.name),
    environment = paste0(loc, "_", Year),
    DMC = as.numeric(BLUEs_DMC_RMA),              
    EMY28PLUS = as.numeric(BLUEs_EMY28PLUS_RMA),      
    EMY35PLUS = as.numeric(BLUEs_EMY35PLUS_RMA),       
    EMY40PLUS = as.numeric(BLUEs_EMY40PLUS_RMA),      
    EMY40_60 = as.numeric(BLUEs_EMY40_60_RMA),         
    EMY60PLUS = as.numeric(BLUEs_EMY60PLUS_RMA),       
    MBD = as.numeric(BLUEs_MBD),                       
    TDMC = as.numeric(BLUEs_TDMC_RMA),                 
    TN = as.numeric(BLUEs_TN_PLANT),                   
    TV = as.numeric(BLUEs_TV),  
    RT = as.numeric(BLUEs_RT), 
    STC = as.numeric(BLUEs_STC),
    STY = as.numeric(BLUEs_STY),                       
    SG = as.numeric(BLUEs_SG_RMA),                       
    DOR = as.numeric(BLUEs_DOR_DAYS),                
    TY = as.numeric(BLUEs_TY),                         
    UWW = as.numeric(BLUEs_UWW)  
  )

Pheno_final = Pheno %>%
  dplyr::filter(grepl("xxxxx", trial)) %>%
  dplyr::filter(!grepl("xx", trial)) %>%
  dplyr::filter(loc %in% c(xxxxx
  )) %>%
  dplyr::filter(QC_test == "trustworthy or unknown" | is.na(QC_test)) %>%
  dplyr::mutate(Env = factor(paste0("env", as.integer(factor(environment)), "-", environment)))


# Training data and testing set IDs
TrainIDs = unique(as.factor(c(Pheno_final$male, Pheno_final$female))) %>% droplevels(c())
