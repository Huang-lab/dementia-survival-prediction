#### utility_functions.R ####
#### Jimmy Zhang @ 1/3/2022 ####

#function for selecting features for each dataset
## @ parameters: (1) data.frame, (2) vector of columns to remove
## @ returns: data frame with selected features
selectFeatures = function(data, remove_cols) {
  #start time
  start.time = Sys.time() 
  print(paste("Feature Selection Start Time:", start.time))
  
  #remove selected columns
  data_cleaned = data %>% select(-all_of(remove_cols))
  
  #remove character/text field variables
  character_vars = names(which(lapply(data_cleaned, class) == "factor" | lapply(data_cleaned, class) == "character"))
  character_vars = character_vars[-1] #keep NACCID
  data_cleaned = data_cleaned %>% select(-all_of(character_vars))
  print(paste("Number of Features Remaining:", length(data_cleaned)))
  
  #print elapsed time
  print(paste("Feature Selection End Time:", Sys.time()))
  print(Sys.time() - start.time)
  
  return(data_cleaned)
}


#function for fixing feature encoding -> convert missing/not assessed/not available codes to NA
## @ parameters: (1) data.frame, cols for MMSE variables (need to be reencoded as well)
## @ returns: data frame with NAs properly encoded
fixVarEncoding = function(data, mmse_cols) {
  start.time = Sys.time() 
  print(paste("Re-Encoding Start Time:", start.time))
  
  #convert all variables except NACCID to numeric
  reencodedData = data.frame(lapply(data[,2:length(data)], as.numeric))
  reencodedData = cbind(data.frame(data$NACCID), reencodedData)
  colnames(reencodedData)[1] = "NACCID"
  
  #convert all -4s and -4.4s (code for not available) to NA
  reencodedData = replace(reencodedData, reencodedData == -4, NA)
  reencodedData = replace(reencodedData, reencodedData == -4.4, NA)
  
  #drop all columns with 100% missing values
  reencodedData = reencodedData %>%
    select(c(names(which(colMeans(is.na(reencodedData)) < 1)))) 
  
  #convert not assessed/not applicable (8/88/888/8888) and unknown (9/99/999/9999) codes to NA
  #if max characters in column = 1, remove 8, 9, if max = 2, remove 88, 99 etc.
  for(i in 6:length(reencodedData)){ #skip header columns
    ##skip NACCETPR
    if(colnames(reencodedData)[i] == "NACCETPR") {
      next
    }
    
    if(nchar(max(reencodedData[i], na.rm = TRUE)) == 1) {
      reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 8 | reencodedData[i] == 9, NA)
    } 
    else if(nchar(max(reencodedData[i], na.rm = TRUE)) == 2) {
      reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 88 | reencodedData[i] == 99, NA)
    }
    else if(nchar(max(reencodedData[i], na.rm = TRUE)) == 3) {
      reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 888 | reencodedData[i] == 999, NA)
    }
    else if(nchar(max(reencodedData[i], na.rm = TRUE)) == 4) {
      reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 8888 | reencodedData[i] == 9999 | reencodedData[i] == 88.8 | reencodedData[i] == 99.9, NA)
    }
    else if(nchar(max(reencodedData[i], na.rm = TRUE)) == 5) {
      reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 888.8 | reencodedData[i] == 999.9, NA)
    }  
  
    #remove extra categories in MMSE variables (95, 96, 97, 98, 995, etc.)
    if(colnames(reencodedData)[i] %in% mmse_cols) {
      if(nchar(max(reencodedData[i], na.rm = TRUE)) == 2) {
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 95, NA)
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 96, NA)
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 97, NA)
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 98, NA)
      }
      else if(nchar(max(reencodedData[i], na.rm = TRUE)) == 3) {
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 995, NA)
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 996, NA)
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 997, NA)
        reencodedData[i] = replace(reencodedData[i], reencodedData[i] == 998, NA)
      }
    }
  }
  
  #fix encoding for race (remove code 50: Other(specify))
  reencodedData$RACE = replace(reencodedData$RACE, reencodedData$RACE == 50, NA)
  
  print(paste("Re-Encoding End Time:", Sys.time()))
  print(Sys.time() - start.time)
  
  return(reencodedData)
}


#function for prepping data for imputation
## @ parameters: (1) data.frame, (2) vector of categorical features
## @ returns: dataset ready for imputation
prepForImputation = function(data, factor_cols) {
  start.time = Sys.time() 
  print(paste("Prepping For Imputation Start Time:", start.time))
  
  
  #remove variables with >40% missing values
  data_filtered = data %>% 
    select(c(names(which(colMeans(is.na(data)) <= 0.4))))

  #convert categorical variables from type numeric to factor
  for(i in 2:length(data_filtered)) { 
    if(names(data_filtered)[i] %in% factor_cols) {
      data_filtered[,i] = as.factor(data_filtered[,i])
    }
  }

  print(paste("Prepping Cohorts End Time:", Sys.time()))
  print(Sys.time() - start.time)
  
  return(data_filtered)
}


#function for imputing each cohort
## @ parameters: (1) data.frame or tbl, (2) quickpred, (3) predictor matrix, (4) int number of iterations
## @ returns: complete, imputed dataset
imputeCohort = function(data = NULL, pred = NULL, predM = NULL, iter = 1) {
  #imputation
  start.time = Sys.time() #track time
  print(paste("Imputation Start Time:", start.time))
  
  if(!is.null(predM)) {   
    tried = try( #try to use default methods for imputation: pmm, logreg, polyreg
        mice(
        data = as.data.frame(data), #convert data to data.frame in case in tbl format
        method = meth,
        predictorMatrix = predM,
        m = 1, 
        seed = 444,
        maxit = iter,
        printFlag = F
      ),
      silent = TRUE
    )
      
    if(inherits(tried, "try-error")) { #if error, use pmm for all vars
      writeLines("Switching to pmm for all.")
      imp_out = mice(
        data = as.data.frame(data),
        defaultMethod = c(rep("pmm", 4)),
        predictorMatrix = predM,
        m = 1, 
        seed = 444,
        maxit = iter,
        printFlag = F
      )
    } else {
      writeLines("Using pmm for numeric, logreg for binary factor, 
                and polyreg for multi-level factor.")
      imp_out = mice(
        data = as.data.frame(data),
        method = meth,
        predictorMatrix = predM,
        m = 1, 
        seed = 444,
        maxit = iter,
        printFlag = F
      )
    }
    
  } else {
    tried = try( #try to use default methods for imputation: pmm, logreg, polyreg
      mice(
        data = as.data.frame(data), #convert data to data.frame in case in tbl format
        method = meth,
        pred = pred,
        m = 1, 
        seed = 444,
        maxit = iter,
        printFlag = F
      ),
      silent = TRUE
    )
    
    if(inherits(tried, "try-error")) { #if error, use pmm for all vars
      writeLines("Switching to pmm for all.")
      imp_out = mice(
        data = as.data.frame(data),
        defaultMethod = c(rep("pmm", 4)),
        pred = pred,
        m = 1, 
        seed = 444,
        maxit = iter,
        printFlag = F
      )
    } else {
      writeLines("Using pmm for numeric, logreg for binary factor, 
                and polyreg for multi-level factor.")
      imp_out = mice(
        data = as.data.frame(data),
        method = meth,
        pred = pred,
        m = 1, 
        seed = 444,
        maxit = iter,
        printFlag = F
      )
    }
  }
    
  #get complete data
  imp_complete = mice::complete(imp_out, include = FALSE) 
  
  #print elapsed time
  print(paste("Imputation End Time:", Sys.time()))
  print(Sys.time() - start.time)
  
  #return imputed dataset
  return(imp_complete)
}


#function for splitting data into 2019 and 2021 datasets
## @ parameters: (1) data.frame, (2) life expectancy (integer), (3) name of new survival variable (string)
## @ returns: a list containing the 2019 and 2021 datasets
splitData = function(data, life_exp, surv_var_name) {
  threshold2019 = 2019.5 - life_exp
  threshold2021 = 2021.89 - life_exp
  
  data$SURV2019 = ifelse(data$VISITYR + (data$VISITMO / 12) > threshold2019,
                               ifelse(data$NACCDIED == 1 & data$NACCYOD + (data$NACCMOD / 12.0) < 2019.5, 1, NA),
                               ifelse(data$NACCYOD + (data$NACCMOD / 12.0) > data$VISITYR + life_exp + data$VISITMO / 12.0, 0, 1))
  
  data$SURV2021 = ifelse(data$VISITYR + (data$VISITMO / 12) > threshold2021,
                         ifelse(data$NACCDIED == 1, 1, NA),
                         ifelse(data$NACCYOD + (data$NACCMOD / 12.0) > data$VISITYR + life_exp + data$VISITMO / 12.0, 0, 1))
  
  #remove NAs
  data_2019 = data[which(!is.na(data$SURV2019)),]
  data_2021 = data[which(is.na(data$SURV2019) & !is.na(data$SURV2021)),]
  
  data_2019 = data_2019 %>%
    select(-SURV2021) %>% 
    rename_with(~surv_var_name, "SURV2019")
  
  data_2021 = data_2021 %>% 
    select(-SURV2019) %>% 
    rename_with(~surv_var_name, "SURV2021")
  
  return(list(data_2019, data_2021))
}


#function for plotting NA threshold vs. number of variables preserved
## @ parameters: (1) data.frame or tbl, (2) vector of thresholds, (3) save (boolean), (4) save file
## @ returns: ggplot
threshVarsPlot = function(data, thresholds, save = FALSE, saveFile = NULL) {
  #keep track of number of variables at each threshold
  varsCount = list()
  i = 1
  
  for(thresh in thresholds) {
    varsCount[i] = length(which(colMeans(is.na(data)) <= thresh))
    i = i + 1
  }
  
  g = ggplot(data = data.frame(thresholds = as.numeric(thresholds), varsCount = as.numeric(varsCount)), aes(x=thresholds, y=varsCount))
  g = g + geom_line(linetype = "dashed", color = "blue") 
  g = g + geom_point()
  g = g + ggtitle("Number of Variables at each Threshold")
  
  if(save) {
    ggsave(saveFile)
  }
  
  return(g)
}