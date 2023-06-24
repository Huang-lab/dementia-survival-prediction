#### Preprocessing.R ####
#### Jimmy Zhang @ 1/3/2022 ####

#set working directory
bdir = "~/Box\ Sync/Huang_lab/manuscripts/NACCAgePrediction/analysis/R/jimmy"
setwd(bdir)

#load packages for preprocessing
library(tidyverse)
library(mice)
library(data.table)
library(crunch)
source("utility_functions.R")

#load 2021 NACC data
data_dir = "~/Box\ Sync/Huang_lab/manuscripts/NACCAgePrediction/data/NACC_NeuroClinical_20211019/"
nacc_data_f = paste(data_dir, "huang10192021.txt.gz", sep = "")
nacc_data = fread(nacc_data_f)

#### Part 1: Feature Selection ####
#(1) remove unnecessary header information
COLS_TO_REMOVE_1 = c("PACKET", "FORMVER")

#(2) remove irrelevant Neuropsychological Battery Summary Scores variables 
COLS_TO_REMOVE_2 = c("NPIQINF", "NPIQINFX")

#add all features to final vector
COLS_TO_REMOVE_ALL = c(COLS_TO_REMOVE_1, COLS_TO_REMOVE_2)

#feature selection
nacc_data_cleaned = selectFeatures(nacc_data, COLS_TO_REMOVE_ALL)

#### Part 2: Fix Feature Encoding ####
MMSE_COLS = c("MMSEORDA", "MMSEORLO", "PENTAGON", "NACCMMSE", "LOGIMEM", "MEMUNITS", "MEMTIME", 
              "UDSBENTC", "UDSBENTD", "UDSBENRS", "DIGIF", "DIGIFLEN", "DIGIB", "DIGIBLEN",
              "ANIMALS", "VEG", "TRAILA", "TRAILARR", "TRAILALI", "TRAILB", "TRAILBRR", "TRAILBLI", 
              "WAIS", "BOSTON", "UDSVERFC", "UDSVERFN", "UDSVERNF", "UDSVERLC", "UDSVERLR",
              "UDSVERLN", "UDSVERTN", "UDSVERTE", "UDSVERTI", "COGSTAT", "NACCC1", "MOCACOMP",
              "MOCAREAS", "MOCALOC", "MOCALAN", "MOCAVIS","MOCAHEAR", "MOCATOTS", "MOCATRAI", 
              "MOCACUBE", "MOCACLOC", "MOCACLON", "MOCACLOH", "MOCANAMI", "MOCAREGI", "MOCADIGI",
              "MOCALETT", "MOCASER7", "MOCAREPE", "MOCAFLUE", "MOCAABST", "MOCARECN", "MOCARECC",
              "MOCARECR", "MOCAORDT", "MOCAORMO", "MOCAORYR", "MOCAORDY", "MOCAORPL", "MOCAORCT",
              "CRAFTVRS", "CRAFTURS", "DIGFORCT", "DIGFORSL", "DIGBACCT", "DIGBACLS", "CRAFTDVR",
              "CRAFTDRE", "CRAFTDTI", "CRAFTCUE", "MINTTOTS", "MINTTOTW", "MINTSCNG", "MINTSCNC",
              "MINTPCNG", "MINTPCNC")

nacc_data_cleaned_reencoded = fixVarEncoding(nacc_data_cleaned, MMSE_COLS)

#### Part 3: Filter rows and prepare cohorts for imputation ####
#colMeans(is.na(nacc_data_cleaned_reencoded))

FACTOR_COLS = c("NACCADC", "NACCREAS", "NACCREFR", "SEX","HISPANIC", "HISPOR", 
                "RACE", "RACESEC", "RACETER", "PRIMLANG", "MARISTAT", "NACCLIVS",  
                "INDEPEND", "RESIDENC", "HANDED", "NACCFAM", "TOBAC30", "TOBAC100", "PACKSPER",
                "CVHATT", "CVAFIB", "CVANGIO", "CVBYPASS", "CVPACE", "CVCHF", "CVOTHR", 
                "CBSTROKE", "PD", "PDOTHR", "CBTIA", "SEIZURES", "NACCTBI", "TRAUMBRF", "TRAUMEXT", "TRAUMCHR",
                "NCOTHR", "DIABETES", "HYPERTEN", "HYPERCHO", "B12DEF", "THYROID", "INCONTU", 
                "INCONTF", "ALCOHOL", "ABUSOTHR", "DEP2YRS", "DEPOTHR", "PSYCDIS", "MEMORY", 
                "ORIENT", "JUDGMENT", "COMMUN", "HOMEHOBB", "PERSCARE", "CDRGLOB", "COMPORT",
                "CDRLANG", "DEL", "DELSEV", "HALL", "HALLSEV", "AGIT", "AGITSEV", "DEPD", "DEPDSEV", 
                "ANX", "ANXSEV", "ELAT", "ELATSEV", "APA", "APASEV", "DISN", "DISNSEV", "IRR", 
                "IRRSEV", "MOT", "MOTSEV", "NITE", "NITESEV", "APP", "APPSEV", "NOGDS", "SATIS", 
                "DROPACT", "EMPTY", "BORED", "SPIRITS","AFRAID", "HAPPY", "HELPLESS", "STAYHOME", 
                "MEMPROB", "WONDRFUL", "WRTHLESS", "ENERGY", "HOPELESS", "BETTER", "BILLS", "TAXES", 
                "SHOPPING", "GAMES", "STOVE", "MEALPREP", "EVENTS", "PAYATTN", "REMDATES", "TRAVEL",
                "DECSUB", "DECIN", "DECCLIN", "COGMEM", "COGJUDG", "COGLANG", "COGVIS", "COGATTN", 
                "COGOTHR", "NACCCOGF", "COGMODE", "BEAPATHY", "BEDEP", "BEVHALL", "BEAHALL", "BEDEL",
                "BEDISIN", "BEIRRIT", "BEAGIT", "BEPERCH", "BEOTHR", "NACCBEHF", "BEMODE", "MOGAIT", 
                "MOFALLS", "MOTREM", "MOSLOW", "NACCMOTF", "MOMODE", "COURSE", "FRSTCHG", "MMSELOC", 
                "MMSELAN", "COGSTAT", "NACCC1", "WHODIDDX", "NORMCOG", 
                "DEMENTED", "NACCPPA", "NACCPPME", "NACCBVFT", "NACCLBDS", "NACCTMCI", "NACCMCIL", 
                "NACCMCIA", "NACCMCIE", "NACCMCIV", "NACCMCII", "IMPNOMCI", "NACCALZD", "NACCALZP", 
                "PROBAD", "PROBADIF", "POSSAD", "POSSADIF", "NACCLBDE", "NACCLBDP", "PARK", "PSP", 
                "PSPIF", "CORT", "CORTIF", "FTD", "FTDIF", "PPAPH", "PPAPHIF", "VASC", "VASCIF", 
                "STROKE", "STROKIF", "DOWNS", "DOWNSIF", "HUNT", "HUNTIF", "PRION", "PRIONIF", 
                "BRNINJ", "BRNINJIF", "HYCEPH", "HYCEPHIF", "NEOP", "NEOPIF", "DEP", "DEPIF", 
                "OTHPSY", "OTHPSYIF", "ALCDEM", "ALCDEMIF", "DYSILL", "DYSILLIF", "MEDS", "MEDSIF", 
                "DEMUN", "DEMUNIF", "COGOTH", "COGOTHIF", "COGOTH2", "COGOTH2F", "COGOTH3", "COGOTH3F", 
                "NACCETPR", "NACCADMU", "NACCFTDM", "NACCUDSD", "NACCDIED", "NACCAUTP", "NACCACTV",
                "NACCNOVS", "NACCNURP", "NACCFTD", "NACCMDSS", "NACCPAFF", "NACCMRSA", "NACCAPSA",
                "NACCNE4S", "NACCAPOE", "COGFLUC", "BEREM", "NACCNORM", "NACCIDEM", "NACCNIHR",
                "NACCLBDM", "NACCACSF", "NACCPCSF", "NACCTCSF")

nacc_data_cleaned_reencoded_prepped = prepForImputation(nacc_data_cleaned_reencoded, FACTOR_COLS)

#### Part 4: Imputation with MICE ####
nacc_data_complete = imputeCohort(data = nacc_data_cleaned_reencoded_prepped, 
                                      pred = quickpred(nacc_data_cleaned_reencoded_prepped, 
                                                       minpuc = 0.25, #only use predictors with absolute correlation of over 0.25 with imputed variable
                                                       mincor = 0.7, #only use predictors with > 70% complete cases
                                                       exclude = c("NACCID") #exclude NACCID from imputation
                                      )
)


#### Part 5: Generate OYS, TYS, FYS, and TenYS datasets ####
#first, merge NACCYOD and NACCMOD back into the dataset
nacc_data_complete_final = merge(nacc_data[,c('NACCYOD', 'NACCMOD', 'NACCID', 'NACCVNUM')], 
                                 nacc_data_complete,
                                 by =c('NACCID', 'NACCVNUM')
                                 )

OYS_datasets = splitData(nacc_data_complete_final, 1, "OYS")
TYS_datasets = splitData(nacc_data_complete_final, 3, "TYS")
FYS_datasets = splitData(nacc_data_complete_final, 5, "FYS")
TenYS_datasets = splitData(nacc_data_complete_final, 10, "TenYS")

OYS_2019_data = OYS_datasets[[1]]
OYS_2021_data = OYS_datasets[[2]]
TYS_2019_data = TYS_datasets[[1]]
TYS_2021_data = TYS_datasets[[2]]
FYS_2019_data = FYS_datasets[[1]]
FYS_2021_data = FYS_datasets[[2]]
TenYS_2019_data = TenYS_datasets[[1]]
TenYS_2021_data = TenYS_datasets[[2]]

#### Part 6: Write imputed datasets ####
write.csv.gz(x = OYS_2019_data, file = "out/imputed_data_v2/OYS/NACC_2019_OYS_data_imputed.csv.gz")
write.csv.gz(x = OYS_2021_data, file = "out/imputed_data_v2/OYS/NACC_2021_OYS_data_imputed.csv.gz")
write.csv.gz(x = TYS_2019_data, file = "out/imputed_data_v2/TYS/NACC_2019_TYS_data_imputed.csv.gz")
write.csv.gz(x = TYS_2021_data, file = "out/imputed_data_v2/TYS/NACC_2021_TYS_data_imputed.csv.gz")
write.csv.gz(x = FYS_2019_data, file = "out/imputed_data_v2/FYS/NACC_2019_FYS_data_imputed.csv.gz")
write.csv.gz(x = FYS_2021_data, file = "out/imputed_data_v2/FYS/NACC_2021_FYS_data_imputed.csv.gz")
write.csv.gz(x = TenYS_2019_data, file = "out/imputed_data_v2/TenYS/NACC_2019_TenYS_data_imputed.csv.gz")
write.csv.gz(x = TenYS_2021_data, file = "out/imputed_data_v2/TenYS/NACC_2021_TenYS_data_imputed.csv.gz")