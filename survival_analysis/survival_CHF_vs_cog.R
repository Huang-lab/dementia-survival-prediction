library("survival")
library("survminer")
library("randomcoloR")
library(data.table)
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442") #### "#0072B2", "#D55E00", "#CC79A7"

#Set working directory
bdir = "/Users/songluo/Library/CloudStorage/Box-Box/Huang_lab/manuscripts/NACCAgePrediction"
setwd(bdir)
#Read Data
naccdata_f = "data/NACC_NeuroClinical_20211019/huang10192021.txt.gz"
naccdata = fread(naccdata_f)

#Get patient's first visit only entries
remove_vector <- c()
for (row in 1:(nrow(naccdata)-1)){
  if (naccdata[row, "NACCID"] == naccdata[row + 1, "NACCID"]) {
    remove_vector <- append(remove_vector, row + 1)
  }
}
naccdata_1stVisit = naccdata[-remove_vector,]
naccdata_1stVisit$time <- ifelse(naccdata_1stVisit$NACCYOD == "8888",(2021 - naccdata_1stVisit$VISITYR) * 365 + (10 - naccdata_1stVisit$VISITMO) * 30, (naccdata_1stVisit$NACCYOD - naccdata_1stVisit$VISITYR) * 365 + (naccdata_1stVisit$NACCMOD - naccdata_1stVisit$VISITMO) * 30)

#divide dataset by CONGHRT status

set.seed(77)

###No CONGHRT
noCHF <- naccdata_1stVisit[naccdata_1stVisit$CONGHRT == "0"]

#survival analysis model
fit.noCHF <- survfit(Surv(time, NACCDIED) ~ CDRGLOB, data = noCHF)
print(fit.noCHF)
noCHF.res.sum <- surv_summary(fit.noCHF)
head(nocaner.res.sum)
ggsurvplot(fit.noCHF,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "Non-CHF Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "0.5", "1", "2", "3"),
           legend.title = "Global CDR Scores",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme

###primary CONGHRT
CHF <-  naccdata_1stVisit[naccdata_1stVisit$CONGHRT == "1"]
#survival analysis model
fit.CHF <- survfit(Surv(time, NACCDIED) ~ CDRGLOB, data = CHF)
print(fit.CHF)
CHF.res.sum <- surv_summary(fit.CHF)
head(CHF.res.sum)
ggsurvplot(fit.CHF,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "CHF Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "0.5", "1", "2", "3"),
           legend.title = "Global CDR Scores",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme


#==========compare with same dementia===============

cdr0 = naccdata_1stVisit[naccdata_1stVisit$CDRGLOB == "0" & 
                           (naccdata_1stVisit$CONGHRT == "0" | naccdata_1stVisit$CONGHRT == "1"  ) ]

#survival analysis model
fit.cdr0 <- survfit(Surv(time, NACCDIED) ~ CONGHRT, data = cdr0)
print(fit.cdr0)
cdr0.res.sum <- surv_summary(fit.cdr0)
head(cdr0.res.sum)
ggsurvplot(fit.cdr0, data = cdr0,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "No Cognitive Impaired Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "1"),
           legend.title = "CONGHRT Status",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme




#CDR=0.5

cdr05 = naccdata_1stVisit[naccdata_1stVisit$CDRGLOB == "0.5" & 
                            (naccdata_1stVisit$CONGHRT == "0" | naccdata_1stVisit$CONGHRT == "1"  ) ]

#survival analysis model
fit.cdr05 <- survfit(Surv(time, NACCDIED) ~ CONGHRT, data = cdr05)
print(fit.cdr05)
cdr05.res.sum <- surv_summary(fit.cdr05)
head(cdr05.res.sum)
ggsurvplot(fit.cdr05, data = cdr05,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "Questionable Impaired Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "1"),
           legend.title = "CONGHRT Status",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme


#CDR=1
cdr1 = naccdata_1stVisit[naccdata_1stVisit$CDRGLOB == "1" & 
                           (naccdata_1stVisit$CONGHRT == "0" | naccdata_1stVisit$CONGHRT == "1"  ) ]

#survival analysis model
fit.cdr1 <- survfit(Surv(time, NACCDIED) ~ CONGHRT, data = cdr1)
print(fit.cdr1)
cdr1.res.sum <- surv_summary(fit.cdr1)
head(cdr1.res.sum)
ggsurvplot(fit.cdr1, data = cdr1,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "Mild Impaired Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "1"),
           legend.title = "CONGHRT Status",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme



#CDR=2
cdr2 = naccdata_1stVisit[naccdata_1stVisit$CDRGLOB == "2" & 
                           (naccdata_1stVisit$CONGHRT == "0" | naccdata_1stVisit$CONGHRT == "1"  ) ]

#survival analysis model
fit.cdr2 <- survfit(Surv(time, NACCDIED) ~ CONGHRT, data = cdr2)
print(fit.cdr2)
cdr2.res.sum <- surv_summary(fit.cdr2)
head(cdr2.res.sum)
ggsurvplot(fit.cdr2, data = cdr2,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "Moderate Impaired Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "1"),
           legend.title = "CONGHRT Status",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme

#CDR3
cdr3 = naccdata_1stVisit[naccdata_1stVisit$CDRGLOB == "3" & 
                           (naccdata_1stVisit$CONGHRT == "0" | naccdata_1stVisit$CONGHRT == "1"  ) ]

#survival analysis model
fit.cdr3 <- survfit(Surv(time, NACCDIED) ~ CONGHRT, data = cdr3)
print(fit.cdr3)
cdr3.res.sum <- surv_summary(fit.cdr3)
head(cdr3.res.sum)
ggsurvplot(fit.cdr3, data = cdr3,
           pval = TRUE, conf.int = TRUE,
           # risk.table = TRUE, # Add risk table
           # risk.table.col = "strata", # Change risk table color by groups
           # linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           break.time.by = 1000,
           title = "Severe Impaired Patients",
           font.main = c(45,"bold","black"),
           font.x = c(35,"bold","black"),
           font.y = c(35,"bold","black"),
           font.tickslab =  c(35,"bold","black"),
           xlab = "Time in Days",
           legend = "right",
           legend.labs = c("0", "1"),
           legend.title = "CONGHRT Status",
           font.legend = c(35,"bold","black"),
           ggtheme = theme_survminer(base_size = 35) + theme(legend.key.size = unit(15,"mm")),
           palette = cbbPalette
) # Change ggplot2 theme


