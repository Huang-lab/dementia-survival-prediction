library("survival")
library("survminer")
library("randomcoloR")
#Set working directory
bdir = "~/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/NACCAgePrediction/"
# analysis/Luo/survival_rate

setwd(bdir)
#Read Data
naccdata_f = "data/NACC_NeuroClinical_20211019/huang10192021.txt.gz"
naccdata_2021 = read.csv(naccdata_f)

#Get patient's first visit only entries
remove_vector <- c()
for (row in 1:(nrow(naccdata_2021)-1)){
  if (naccdata_2021[row, "NACCID"] == naccdata_2021[row + 1, "NACCID"]) {
    remove_vector <- append(remove_vector, row + 1)
  }
}
naccdata_1stVisit = naccdata_2021[-remove_vector,]
naccdata_1stVisit$time <- ifelse(naccdata_1stVisit$NACCYOD == "8888",(2021 - naccdata_1stVisit$VISITYR) * 365 + (10 - naccdata_1stVisit$VISITMO) * 30, (naccdata_1stVisit$NACCYOD - naccdata_1stVisit$VISITYR) * 365 + (naccdata_1stVisit$NACCMOD - naccdata_1stVisit$VISITMO) * 30)

table(naccdata_1stVisit$NACCETPR)

#get rid of sub-types with n < 100
naccetpr.group.remove <- c("3","6","9","10","11","14","15","16","17","20","21","22","24","26","27","29","99") #include 88 for comparison
remove_vector = c()
for (row in 1:(nrow(naccdata_1stVisit))){
  for ( i in naccetpr.group.remove) {
    if (naccdata_1stVisit[row,"NACCETPR"] == i) {
      remove_vector <- append(remove_vector, row)
    }
  }
}
naccdata_1stVisit = naccdata_1stVisit[-remove_vector,]


#survival analysis model
fit <- survfit(Surv(time, NACCDIED) ~ NACCETPR, data = naccdata_1stVisit)
print(fit)
res.sum <- surv_summary(fit)
head(res.sum)
ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "Survival Analysis Based On Dementia Subtypes",
           font.main = c(30,"bold","black"),
           font.x = c(20,"bold","black"),
           font.y = c(20,"bold","black"),
           font.tickslab =  c(20,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("Alzheimer’s Disease", "Lewy Body Disease", 
                           "Progressive Supranuclear \nPalsy", 
                           "Corticobasal Degeneration", "FTLD, other", 
                           "Vascular Brain Injury",
                           "Prion Disease", "Traumatic Brain Injury",
                           "Other Neurologic, Genetic, \nor Infectious condition", 
                           "Depression", "Other Psychiatric Disease",
                           "Systemic Disease or Medical Illness", 
                           "Other Specified Reasons", 
                           "Not Applicable, \nNot Cognitively Impaired"),
           legend.title = "Primary Etiologic Diagnosis",
           font.legend = c(20,"bold","black"),
           ggtheme = theme_survminer(base_size = 20) + theme(legend.key.size = unit(15,"mm")),
           palette = distinctColorPalette(14)
) # Change ggplot2 theme
