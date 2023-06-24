library("survival")
library("survminer")
library("randomcoloR")
library(data.table)
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

#survival analysis model
fit <- survfit(Surv(time, NACCDIED) ~ CDRGLOB, data = naccdata_1stVisit)
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
           title = "Survival Analysis Based On Global CDR Scores",
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
           palette = distinctColorPalette(5)
) # Change ggplot2 theme
