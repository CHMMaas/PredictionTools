# clean environment
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) dev.off()
cat('\014')

###
### Validation of Epic's Sepsis Model, prospective evaluation
### Author: C.H.M. Maas
###
library(dplyr)
library(caret)     # sensitivity, specificity, PPV, NPV
library(ggplot2)
library(gridExtra) # risk table
library(cowplot)   # add two plots together
source("C:/Users/carol/OneDrive - Erasmus MC/Projects Tufts/Project Epic Prediciton Models/Project Validation Epic Sepsis Model/Validation 2025/Sepsis Evaluation/Code/HelperFunctions.R")
file.path <- "C:/Users/carol/OneDrive - Erasmus MC/Projects Tufts/Project Epic Prediciton Models/Project Validation Epic Sepsis Model/Validation 2025/Sepsis Evaluation/"
load(paste0(file.path, "Data/patient.data.RData"))

# aggregated full data set
times <- patient.data |> filter(max.time.point==time.point) |> 
  select(admission.date.and.time,
         discharge.date.and.time,
         time.sep3.outcome) |>
  mutate(hosp.time=as.numeric(difftime(discharge.date.and.time,
                                       admission.date.and.time,
                                       units="hours")))
hist(times$hosp.time, 
     xlab="Time in hours", ylab="Frequency", breaks=500, 
     main="Time of hospitalization")

times.sepsis <- patient.data |> filter(time.sep3.outcome==1) |>
  pull(max.time.point)/4/24
hist(times.sepsis, 
     xlab="Time in days", ylab="Frequency", breaks=500, 
     main="Times at which sepsis occurred")
cat("Total number of sepsis:", length(times.sepsis), 
    "(", round(length(times.sepsis)/length(unique(patient.data$LCSN))*100, 1), "%) \n", 
    "Number of cases of sepsis up to 12 hours from admission:", 
    sum(times.sepsis<0.5), 
    "(", round(sum(times.sepsis<0.5)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis between 12-24 hours after admission:", 
    sum(times.sepsis>=0.5 & times.sepsis<1), 
    "(", round(sum(times.sepsis>=0.5 & times.sepsis<1)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis on the second day:", 
    sum(times.sepsis>=1 & times.sepsis<2),
    "(", round(sum(times.sepsis>=1 & times.sepsis<2)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis on the third day:", 
    sum(times.sepsis>=2 & times.sepsis<3), 
    "(", round(sum(times.sepsis>=2 & times.sepsis<3)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis on the fourth day:", 
    sum(times.sepsis>=3 & times.sepsis<4), 
    "(", round(sum(times.sepsis>=3 & times.sepsis<4)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis on the fifth day:", 
    sum(times.sepsis>=4 & times.sepsis<5), 
    "(", round(sum(times.sepsis>=4 & times.sepsis<5)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis on the sixth day:", 
    sum(times.sepsis>=5 & times.sepsis<6), 
    "(", round(sum(times.sepsis>=5 & times.sepsis<6)/length(times.sepsis)*100, 1), "% of positives)", "\n",
    "Number of cases of sepsis after the sixth day:", 
    sum(times.sepsis>=6), 
    "(", round(sum(times.sepsis>=6)/length(times.sepsis)*100, 1), "% of positives)", "\n")

# report descriptives
cat(" Number of rows:", nrow(patient.data), "\n",
    "Number of unique hospitalizations:", length(unique(patient.data$LCSN)), "\n",
    "Number of unique individuals:", length(unique(patient.data$MRN)), "\n",
    "Average time of hospitalization:", round(mean(times$hosp.time), 1), "hours \n",
    "Average time of hospitalization for negative patients:", 
    round(mean(times$hosp.time[times$time.sep3.outcome==0]), 1), "hours \n",
    "Average time of hospitalization for positive patients:", 
    round(mean(times$hosp.time[times$time.sep3.outcome==1])/24, 1), "days \n",
    "Median time of hospitalization:", round(median(times$hosp.time), 1), "hours \n",
    "Median time of hospitalization for negative patients:", 
    round(median(times$hosp.time[times$time.sep3.outcome==0]), 1), "hours \n",
    "Median time of hospitalization for positive patients:", 
    round(median(times$hosp.time[times$time.sep3.outcome==1])/24, 1), "days \n",
    "Maximum time of hospitalization for positive patients:", 
    round(max(times$hosp.time[times$time.sep3.outcome==1])/24, 1), "days \n")

###
### 2-hour look back for predictions 
### 8-hours look forward for outcome
###
final.threshold <- 25
look.back <- 2    # in hours 
look.forward <- 8 # in hours 

# calculate AUC with look-back and look-forward period
Prevalence <- c()
at.risk <- c()
TPs <- c()
FPs <- c()
TNs <- c()
FNs <- c()
AUC <- c()
Sensitivity <- c()
Specificity <- c()
PPV <- c()
NPV <- c()
max.days <- 4 # in days
time.points <- seq(0, max.days*24, by=0.25)  # in hours
for (current.time.point in time.points){
  # DISCUSS: take into account prediction at the time point?
  # extract maximum prediction in look-back period
  predictions <- patient.data |> 
    group_by(LCSN) |> 
    filter(max.time.point >= current.time.point*4 &
             (time.point >= max(0, current.time.point*4 - look.back*4) & 
                time.point <= current.time.point*4)) |> 
    summarize(aggr.pred = max(model.score))
  
  # DISCUSS: take into account outcome at the time point?
  # extract outcome in look-forward period
  outcomes <- patient.data |>
    group_by(LCSN) |>
    filter(max.time.point >= current.time.point*4 &
             (time.point >= current.time.point*4 & 
                time.point <= current.time.point*4 + look.forward*4)) |>
    summarize(aggr.outcome = max(time.sep3.outcome))
  
  # join predictions and outcomes
  joined <- full_join(predictions, outcomes, by="LCSN")
  
  # extract predictions and labels  
  labels <- joined$aggr.outcome
  predictions <- as.numeric(joined$aggr.pred>=final.threshold/100)
  
  # save AUC
  if (sum(labels)==0){
    AUC <- c(AUC, 0.5)
  } else{
    AUC <- c(AUC, as.numeric(pROC::auc(labels, predictions, quiet=TRUE)))
  }
  
  # save other measures
  TP <- sum(predictions==1 & labels==1)
  FP <- sum(predictions==1 & labels==0)
  TN <- sum(predictions==0 & labels==0)
  FN <- sum(predictions==0 & labels==1)
  Sensitivity <- c(Sensitivity, TP/(TP+FN))
  Specificity <- c(Specificity, TN/(TN+FP))
  PPV <- c(PPV, TP/(TP+FP))
  NPV <- c(NPV, TN/(TN+FN))
  at.risk <- c(at.risk, nrow(joined))
  TPs <- c(TPs, TP)
  FPs <- c(FPs, FP)
  TNs <- c(TNs, TN)
  FNs <- c(FNs, FN)
  Prevalence <- c(Prevalence, sum(labels==1)/nrow(joined)*100)
}
###
### check the above
###
# try.ID <- 34371223
# # full record
# View(patient.data |> filter(LCSN==try.ID) |>
#        select(LCSN, time.point, model.score, time.sep3.outcome))
# # look back
# View(patient.data |>
#        filter(LCSN==try.ID &
#                 max.time.point >= current.time.point*4 &
#                 (time.point >= max(0, current.time.point*4 - look.back*4) &
#                    time.point <= current.time.point*4)) |>
#        select(LCSN, time.point, model.score, time.sep3.outcome))
# print(predictions |> filter(LCSN==try.ID))
# # look forward
# View(patient.data |>
#        filter(LCSN==try.ID &
#                 max.time.point >= current.time.point*4 &
#                 (time.point >= current.time.point*4 &
#                    time.point <= current.time.point*4 + look.forward*4)) |>
#        select(LCSN, time.point, model.score, time.sep3.outcome))
# print(outcomes |> filter(LCSN==try.ID))

# save results in data.frame
results.df <- data.frame(Time = time.points/24, # in days
                         Prevalence = Prevalence,
                         Risk = at.risk,
                         TP = TPs,
                         FP = FPs,
                         TN = TNs,
                         FN = FNs,
                         AUC = AUC,
                         Sensitivity = Sensitivity,
                         Specificity = Specificity,
                         PPV = PPV,
                         NPV = NPV)

# risk table
risk_table <- ggplot(results.df |> 
                       filter(Time %in% seq(0, max.days, by=0.5)) |> 
                       select(Time, Prevalence, Risk, TP, FP, TN, FN), 
                     aes(x=Time, y=0)) +
  geom_text(aes(y=5, label=round(Prevalence, 1)), size=3.5) +
  geom_text(aes(y=4, label=Risk), size=3.5) +
  geom_text(aes(y=3, label=TP), size=3.5) +
  geom_text(aes(y=2, label=FP), size=3.5) +
  geom_text(aes(y=1, label=TN), size=3.5) +
  geom_text(aes(y=0, label=FN), size=3.5) +
  scale_x_continuous(breaks = seq(0, max.days, by=0.5)) +
  scale_y_continuous(breaks = 0:5,
                     labels = c("FN", "TN", "FP", "TP", "# at risk", "Prevalence")) +
  labs(title=NULL, 
       y=NULL,
       x="Time (in days)", y=NULL) +
  theme(plot.background = element_rect(fill="white", linewidth=0),
        panel.background = element_blank(),
        panel.grid = element_blank(),
        axis.line.x = element_line(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10, angle=0, vjust=0.6),
        axis.ticks.y = element_blank())

# plots 
plot_list <- list()
for (measure in c("Prevalence", "AUC", 
                  "Sensitivity", "Specificity", 
                  "PPV", "NPV")){
  # extract measure
  results.measure <- results.df[, c("Time", measure)]
  colnames(results.measure) <- c("Time", "Measure")
  
  # make plot for measure
  plot_measure <- ggplot(results.measure, 
                         aes(x = Time, y = Measure)) +
    geom_line(color = "black", linewidth = 0.5) +  
    scale_x_continuous(breaks = seq(0, max.days, 1)) + # set x-axis ticks
    labs(title = paste(measure, "calculated at each hour"), 
         x = NULL, 
         y = measure) +
    theme_minimal() +
    theme(plot.background = element_rect(fill="white", linewidth=0, color="white"),
          axis.text = element_text(size = 10))
  
  # adjust y-axis
  if (measure=="AUC"){
    plot_measure <- plot_measure + 
      geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +  # minimum AUC
      scale_y_continuous(limits=c(0.45, 1),
                         breaks=seq(0.5, 1, 0.1))
  } else{
    plot_measure <- plot_measure +
      ylim(0, 1)
  }
  
  # save plot
  plot_list[[measure]] <- plot_measure
}

# save plots
ggsave(paste0(file.path, "Results/Time dependent measures at a threshold of ", final.threshold, ".png"), 
       plot = ggpubr::ggarrange(plot_list[["Prevalence"]], 
                                plot_list[["AUC"]], 
                                plot_list[["Sensitivity"]], 
                                plot_list[["Specificity"]]+scale_y_continuous(limits=c(0.75, 1)), 
                                plot_list[["PPV"]]+scale_y_continuous(limits=c(0, 0.25)), 
                                plot_list[["NPV"]]+scale_y_continuous(limits=c(0.9, 1)), 
                                risk_table, nrow=7, ncol = 1, align="hv"),
       width = 6, height = 12, dpi = 300)