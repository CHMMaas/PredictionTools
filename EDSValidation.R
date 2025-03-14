# clean environment
rm(list = ls(all.names = TRUE))
if(!is.null(dev.list())) dev.off()
cat('\014')

###
### Validation of Epic's Sepsis Model
### Author: C.H.M. Maas
###
# install PredictionTools: https://github.com/CHMMaas/PredictionTools
# remotes::install_github("CHMMaas/PredictionTools")

### libraries
library(dplyr)
library(ggplot2)
library(PredictionTools)
source("C:/Users/carol/OneDrive - Erasmus MC/Projects Tufts/Project Epic Prediciton Models/Project Validation Epic Sepsis Model/Validation 2025/Sepsis Evaluation/Code/HelperFunctions.R")

### Make sure there are three folders: Code, Data, Results to read from and write to
file.path <- "C:/Users/carol/OneDrive - Erasmus MC/Projects Tufts/Project Epic Prediciton Models/Project Validation Epic Sepsis Model/Validation 2025/Sepsis Evaluation/"

### load data
original.data <- read.csv(paste0(file.path, "Data/SepsisData(X).csv"))
head(original.data)
# predictions (model.score) are generated every 15 minutes
# MRN = patient ID
# CSN = hospitalization ID
# LCSN = Linked CSN, transfers are not considered separate hospitalizations
colnames(original.data) <- c("MRN", "CSN", "LCSN", "sex", "race", "ethnicity", "age", 
                             "admission.date.and.time", 
                             "discharge.date.and.time", "center",
                             "prediction.date.and.time", "model.score",
                             "model.prediction", "sep3.time.zero", "sep3.outcome")

###
### check for duplicates
###
# if the sum is zero there are no duplications so no need to adjust
# sum(duplicated(original.data))

###
### Add columns to data
###
# separate date and time
patient.data <- cbind(original.data, 
                      reshape2::colsplit(original.data$admission.date.and.time, " ", 
                                         names=c("admission.date", "admission.time")),
                      reshape2::colsplit(original.data$discharge.date.and.time, " ", 
                                         names=c("discharge.date", "discharge.time")),
                      reshape2::colsplit(original.data$prediction.date.and.time, " ", 
                                         names=c("prediction.date", "prediction.time")))

###
### Admission inclusion criteria
###
# 16 August is the date at which predictions started generating
cat("Number of patients admitted before 16 August:", sum(patient.data |> 
                                                           filter(admission.date<"2024-08-16") |> 
                                                           summarize(unique_ids=n_distinct(LCSN)) |> 
                                                           pull(unique_ids)), "\n",
    "Number of patients admitted after 31 October:", sum(patient.data |> 
                                                             filter(admission.date>"2024-10-31") |> 
                                                             summarize(unique_ids=n_distinct(LCSN)) |> 
                                                             pull(unique_ids)), "\n")
patient.data <- patient.data |> 
  filter(admission.date>="2024-08-16" & admission.date<"2024-11-01")
cat(" Number of rows:", nrow(patient.data), "\n",
    "Number of unique hospitalizations:", length(unique(patient.data$LCSN)), "\n",
    "Number of unique individuals:", length(unique(patient.data$MRN)), "\n")

# For negative patients predictions are pulled from admission until discharge
# For positive patients predictions are pulled form admission until sepsis onset
# TODO: move this down?
cat('Admission dates lie between', 
    min(patient.data$admission.date), 'and', 
    max(patient.data$admission.date), '\n',
    'Discharge dates lie between', 
    min(patient.data$discharge.date), 'and', 
    max(patient.data$discharge.date), '\n',
    'Prediction dates lie between', 
    min(patient.data$prediction.date), 'and', 
    max(patient.data$prediction.date), '\n')

# Check that full hospitalization was recorded for the patient with discharge at 28 January
View(original.data[which(original.data$MRN=="34145872"),])
# TODO: check the patient with predictions until 21 January

# set time points and maximum time point
patient.data <- patient.data |> 
  group_by(LCSN) |> 
  mutate(time.point=0:(n()-1),
         max.time.point=max(time.point))

# make outcome column time-dependent
patient.data$time.sep3.outcome <- rep(0, nrow(patient.data))
patient.data[which(patient.data$time.point==patient.data$max.time.point & 
                     patient.data$sep3.outcome==1), "time.sep3.outcome"] <- 1

# look at a patient without the outcome
# View(patient.data |> filter(MRN==31225483))
# look at a patient with the outcome
# View(patient.data |> filter(MRN==31259363))
# look at a patient who had two hospitalizations one month apart
# View(patient.data |> filter(MRN==31228576))
# look at patient who had transfer
# View(patient.data |> filter(MRN==31819490))

# Exclude patients that immediately have the outcome at the first prediction
cat("There are", sum(patient.data$max.time.point==0 & patient.data$sep3.outcome==1),
    "(",
    round(sum(patient.data$max.time.point==0 & patient.data$sep3.outcome==1)/sum(patient.data$time.sep3.outcome==1)*100, 1), 
    "% of positives)",
    "patients that experience the outcome immediately at the time of the first prediction.\n")
patient.data <- patient.data |> filter(!(max.time.point==0 & sep3.outcome==1))
cat(" Number of rows:", nrow(patient.data), "\n",
    "Number of unique hospitalizations:", length(unique(patient.data$LCSN)), "\n",
    "Number of unique individuals:", length(unique(patient.data$MRN)), "\n")

# look at patient who had two hospitalizations and the second was immediately sepsis
# View(original.data |> filter(MRN==31235059))
# View(patient.data |> filter(MRN==31235059))

# make predictions between 0 and 100
patient.data$model.score <- patient.data$model.score/100

###
### CENTER
###
patient.data$center.grouped <- ifelse(patient.data$center=="LOWELL GENERAL HOSPITAL MAIN CAMPUS" | patient.data$center=="LOWELL GENERAL HOSPITAL SAINTS CAMPUS", "Lowell",
                                       ifelse(patient.data$center=="MELROSE WAKEFIELD HOSPITAL", "Melrose",
                                              ifelse(patient.data$center=="TUFTS MEDICAL CENTER", "Tufts",
                                                     ifelse(patient.data$center=="HHF HSPC High Pointe House" | patient.data$center=="LAWRENCE MEMORIAL HOSPITAL", "Other", NA))))

###
### SEX
###
patient.data <- patient.data |> mutate(sex = na_if(sex, "Unknown"),
                                         sex = na_if(sex, "X"),
                                         sex = na_if(sex, "Nonbinary"))

###
### RACE
###
patient.data$race.groups <- ifelse(patient.data$race=="White", "White",
                                   ifelse(patient.data$race=="Asian", "Asian",
                                          ifelse(patient.data$race=="Black or African American", "Black",
                                                 ifelse(patient.data$race=="Other" | 
                                                          patient.data$race=="American Indian or Alaska Native" | 
                                                          patient.data$race=="Native Hawaiian / Pacific Islander", 
                                                        "Other", NA))))

###
### AGE
###
# set John Doe's age to NA
patient.data <- patient.data |> mutate(age = na_if(age, 125),
                                       age = na_if(age, 124))

# create age groups
quantiles.age <- quantile(aggregate(patient.data |> select(age), 
                                    list(patient.data$LCSN), 
                                    max)[, "age"], na.rm=TRUE)
young.age <- quantiles.age[2]
old.age <- quantiles.age[4]
patient.data$age.groups <- rep(paste("Between", young.age, "and",
                                      old.age, "years"), 
                               nrow(patient.data))
patient.data$age.groups <- ifelse(patient.data$age<young.age, 
                                  paste("Younger than", young.age, "years"),
                                   ifelse(patient.data$age>old.age, 
                                          paste("Older than", old.age, "years"),
                                          ifelse(is.na(patient.data$age), NA, 
                                                 paste("Between", young.age, "and", old.age, "years"))))

###
### Data checks
###
# for all patients, check admission time and first prediction
first.prediction <- original.data |> 
  group_by(LCSN) |> 
  slice(1) |> 
  filter(admission.date.and.time>"2024-08-16" & admission.date.and.time<"2024-09-30") |>
  select(admission.date.and.time, prediction.date.and.time)
diff.admission.first.prediction <- as.numeric(difftime(first.prediction$prediction.date.and.time,
                                                       first.prediction$admission.date.and.time,
                                                      units="hours"))
hist(diff.admission.first.prediction, xlab="Time in hours", 
     main="Histogram of time between admission and first prediction")

# investigate maximum time between admission and first prediction
max.ID.admission <- as.numeric(first.prediction[which(diff.admission.first.prediction==max(diff.admission.first.prediction)), "LCSN"])
cat("Maximum time between admission and first prediction is", 
    round(max(diff.admission.first.prediction), 1), "hours \n")
# View(original.data |> filter(LCSN==max.ID.admission) |>
#      select(MRN, CSN, LCSN, admission.date.and.time, prediction.date.and.time))

# for negatives, check the time between last prediction and discharge
df.negatives <- patient.data |> filter(sep3.time.zero=="NULL" & time.point==max.time.point)
diff.discharge.last.prediction <- as.numeric(difftime(df.negatives$discharge.date.and.time, 
                                                      df.negatives$prediction.date.and.time,
                                                      units="hours"))
hist(diff.discharge.last.prediction, xlab="Time in hours", 
     main="Histogram of time between last prediction and discharge")

# investigate maximum time between last prediction and discharge
max.ID.discharge <- as.numeric(df.negatives[which(diff.discharge.last.prediction==max(diff.discharge.last.prediction)), "LCSN"])
cat("Maximum time between last prediction and discharge is", 
    round(max(diff.discharge.last.prediction), 1), "hours \n")
# View(original.data |> filter(LCSN==max.ID.discharge) |>
#        select(MRN, CSN, LCSN, discharge.date.and.time, prediction.date.and.time))

# for positives, check the time between last prediction and sepsis outcome
df.positives <- patient.data |> filter(time.sep3.outcome==1)
diff.outcome.last.prediction <- as.numeric(difftime(df.positives$sep3.time.zero, 
                                                    df.positives$prediction.date.and.time,
                                                    units="hours"))
hist(diff.outcome.last.prediction, xlab="Time in hours", 
     main="Histogram of time between last prediction and outcome")

# investigate maximum time between last prediction and sepsis outcome
max.ID.outcome <- as.numeric(df.positives[which(diff.outcome.last.prediction==max(diff.outcome.last.prediction)), "LCSN"])
cat("Maximum time between last prediction and sepsis outcome is", 
    round(max(diff.outcome.last.prediction), 1), "hours \n")
# View(original.data |> filter(LCSN==max.ID.outcome) |>
#        select(MRN, CSN, LCSN, sep3.time.zero, prediction.date.and.time))

# missing predictions
cat(" For", sum(diff.admission.first.prediction>2), "(",
    round(sum(diff.admission.first.prediction>2)/length(unique(patient.data$LCSN))*100, 1), "%)",
    "patients the time between admission and their first prediction is more than 2 hours. \n",
    "For", sum(diff.discharge.last.prediction>1), "(",
    round(sum(diff.discharge.last.prediction>1)/nrow(df.negatives)*100, 1), "%)",
    "negative patients the time between discharge and the last prediction is more than 1 hour. \n",
    "For", sum(diff.outcome.last.prediction>1), "(",  
    round(sum(diff.outcome.last.prediction>1)/nrow(df.positives)*100, 1), "%)",
    "positive patients the time between sepsis and the last prediction is more than 1 hour. \n")

# exclude patients with missing predictions
IDs <- c(c(first.prediction[diff.admission.first.prediction>2,"LCSN"])$LCSN,
         c(df.positives[diff.outcome.last.prediction>1, "LCSN"])$LCSN,
         c(df.negatives[diff.discharge.last.prediction>1, "LCSN"])$LCSN)
patient.data <- patient.data |> filter(!(LCSN %in% IDs))
cat(" Number of rows:", nrow(patient.data), "\n",
    "Number of unique hospitalizations:", length(unique(patient.data$LCSN)), "\n",
    "Number of unique individuals:", length(unique(patient.data$MRN)), "\n")
save(patient.data, file=paste0(file.path, "Data/patient.data.RData"))

###
### SENSITIVITY ANALYSIS: SEPSIS-FREE SURVIVAL 
### 
# patient.data <- patient.data |>
#   group_by(LCSN) |>
#   filter(time.point >= 12*4)

###
### 8-hour time window for patients with the outcome
###
# Negative patients: full hopsitalization
# Positive patients: only those predictions 8 hours before sepsis onset
patient.data.8h <- patient.data |> 
  group_by(LCSN) |> 
  filter(sep3.outcome==0 | (sep3.outcome==1 & time.point >= max.time.point - 8*4))

# look at a patient without the outcome
# View(patient.data.8h |> filter(MRN==31226594))
# look at a patient with the outcome
# View(patient.data.8h |> filter(MRN==31252158))

###
### Aggregate data
###
aggregated.data <- aggregate(patient.data.8h, list(patient.data.8h$LCSN), max)

###
### Descriptive statistics
###
# center
aggregated.data$center.grouped <- as.factor(aggregated.data$center.grouped)
table(aggregated.data$center.grouped)
cat("Number of unknown center:", sum(is.na(aggregated.data$center.grouped)), "\n")

# sex
table(aggregated.data$sex)
round(table(aggregated.data$sex)/nrow(aggregated.data)*100, 0)
cat("Number of unknown sex:", sum(is.na(aggregated.data$sex)), "\n")

# age
cat("Maximum age:", max(aggregated.data$age, na.rm=TRUE), "\n")
hist(aggregated.data$age, xlab="Age", 
     main="Histogram of age (in years)")
quantile(aggregated.data$age, na.rm=TRUE)
cat("Number of unknown age:", sum(is.na(aggregated.data$age.groups)), "\n")

# race
sort(table(aggregated.data$race.groups))
cat("Number of unknown race:", sum(is.na(aggregated.data$race.groups)), "\n")

# report outcomes
cat("Number of events", sum(patient.data.8h$time.sep3.outcome), "(", 
    round(sum(patient.data.8h$time.sep3.outcome)/nrow(aggregated.data)*100, 1), "%) \n")
hist(original.data$model.score/100, xlab="Risk score", 
     main="Histogram of risk scores of full data set")
hist(aggregated.data$model.score, xlab="Risk score", 
     main="Histogram of risk scores of aggregated data set (maximum prediction)")
round(quantile(aggregated.data$model.score)*100, 0)

###
### calculate lead time: time from alarm until sepsis onset
###
# set final threshold as recommended by Epic
final.threshold <- 25

# obtain times at which outcomes occur
# results in dataframe with one row per hospitalization and time at which sepsis occured
outcomes <- patient.data.8h[patient.data.8h$time.sep3.outcome==1, 
                            c("LCSN", "time.point")]

# plot median lead time across range of thresholds
plot.lead.times <- c()
threshold.range <- 0:80
# calculate lead time in hours for different thresholds
for (try.threshold in threshold.range){
  # indicate when alarm went off by new variable 'alarm'
  patient.data.lead.time <- patient.data.8h |> 
    group_by(LCSN) |> 
    mutate(alarm=as.numeric(model.score>try.threshold/100))
  
  # at each threshold, select alarms and the time at which that alarm occured
  alarms <- patient.data.lead.time[patient.data.lead.time$alarm==1, 
                                   c("LCSN", "time.point")]
  
  # at each threshold, select only first alarm and the time at which first alarm occured
  first.alarms <- alarms |> 
    group_by(LCSN) |> 
    slice(1)
  
  # join two tables
  joined <- right_join(first.alarms, outcomes, by="LCSN")
  colnames(joined) <- c("LCSN", "time.first.alarm", "time.outcome")
  
  if (try.threshold==final.threshold){
    cat("At a threshold of", final.threshold, ":", "\n",
        "Number of patients without the outcome:", nrow(aggregated.data)-nrow(outcomes), "\n",
        "Number of patients with the outcome:", nrow(outcomes), "\n",
        "Number of patients with alarm before outcome:", sum(is.na(joined$time.first.alarm)), "\n",
        "Number of patients without alarm before outcome:", sum(!is.na(joined$time.first.alarm)), "\n")
  }
  
  # only keep those patients that had an alarm before the outcome 
  joined <- joined |> filter(!is.na(time.first.alarm))
  
  # calculate lead times in hours
  lead.times <- (joined$time.outcome-joined$time.first.alarm)/4
  
  # plot median lead times
  plot.lead.times <- c(plot.lead.times, median(lead.times)) 
}
table.lead.times <- data.frame(threshold=threshold.range, lead.time=plot.lead.times)
print(table.lead.times)

# save plot
ggsave(filename = paste0(file.path, "Results/Lead times.png"), 
       plot = ggplot(table.lead.times, 
                     aes(x = threshold.range, y = lead.time)) +
         geom_line() + 
         geom_hline(yintercept = 1, color = "red") + 
         ylim(0, ceiling(max(table.lead.times$lead.time))) +
         labs(x = "Threshold", y = "Median lead time") +
         theme_bw() +
         theme(text=element_text(size=20)), 
       width = 5, height = 5, units = "in", dpi = 300)

# report median lead time at a threshold of 25
cat("Median lead time at threshold of 25 is", 
    as.numeric(table.lead.times |> filter(threshold==25) |> select(lead.time)), 
    "hours. \n")
###
### Decision-curve analysis
###
grDevices::png(file=paste0(file.path, "Results/DCA.png"),
               width=500, height=500, units="px")
dcurves::dca(sep3.outcome ~ model.score, 
             data=aggregated.data) |>
  plot(smooth=TRUE)
grDevices::dev.off()

###
### Overall alert fatigue
###
cat("Number of alarms with threshold of", final.threshold, 
    "and higher:", sum(aggregated.data$model.score>final.threshold/100),
    "(", sum(aggregated.data$model.score>final.threshold/100)/nrow(aggregated.data)*100, "% ) \n",
    "False alarms:", sum(aggregated.data$model.score>final.threshold/100&aggregated.data$sep3.outcome==0),
    "(", sum(aggregated.data$model.score>final.threshold/100&aggregated.data$sep3.outcome==0)/sum(aggregated.data$model.score>final.threshold/100)*100, "% ) \n")

###
### Average number of alarms in September
###
# add alarms to patient data
patient.data.8h <- patient.data.8h |> 
  group_by(LCSN) |> 
  mutate(alarm=as.numeric(model.score>final.threshold/100))

# select Tufts patients in September
Tufts.sep.patients <- patient.data.8h |>
  filter(center=="TUFTS MEDICAL CENTER") |>
  filter(prediction.date>="2024-09-01" & prediction.date<="2024-09-30")

# report
cat(" There were", sum(Tufts.sep.patients$alarm), "alarms in September at Tufts medical center \n",
    "Of which", sum(Tufts.sep.patients$alarm & Tufts.sep.patients$sep3.outcome==1), "were correct alarms \n",
    "Of which", sum(Tufts.sep.patients$alarm & Tufts.sep.patients$sep3.outcome==0), "were false alarms \n",
    "On average, there were", round(sum(Tufts.sep.patients$alarm)/30, 1), "alarms per day \n",
    "Of which", round(sum(Tufts.sep.patients$alarm & Tufts.sep.patients$sep3.outcome==1)/30, 1), "were correct alarms per day \n",
    "Of which", round(sum(Tufts.sep.patients$alarm & Tufts.sep.patients$sep3.outcome==0)/30, 1), "were false alarms per day \n")

###
### Plot measures for different thresholds
###
nboot <- 10 # FINAL: nboot = 1000
threshold <- final.threshold/100 # TODO: make prettier
measures.table <- data.frame(names=c("n", "p", "sensitivity", "specificity",
                                     "PPV", "NPV", "Cindex", "Eavg", "E90"))
subgroups <- c("All",  "Male", "Female", "Tufts", "Melrose", "Lowell",
               paste("Younger than", young.age, "years"), 
               paste("Between", young.age, "and", old.age, "years"), 
               paste("Older than", old.age, "years"),
               "White", "Asian", "Black", "Other")
for (subgroup in subgroups){
  cat("Now calculating for", subgroup, "\n")
  
  # select the relevant hospitalizations
  if (subgroup=="All"){
    # Note: we're using the maximum prediction for each hospitalization
    data.subgroup <- aggregated.data
  } else if (subgroup=="Male" | subgroup=="Female"){
    data.subgroup <- aggregated.data |> filter(sex==subgroup)
  } else if (subgroup=="Tufts" | subgroup=="Melrose" | subgroup=="Lowell"){
    data.subgroup <- aggregated.data |> filter(center.grouped==subgroup)
  } else if (subgroup==paste("Younger than", young.age, "years") |
             subgroup==paste("Between", young.age, "and", old.age, "years") |
             subgroup==paste("Older than", old.age, "years")){
    data.subgroup <- aggregated.data |> filter(age.groups==subgroup)
  } else if (subgroup=="White" | subgroup=="Asian" | subgroup=="Black" |
             subgroup=="Other" | subgroup=="Unknown"){
    data.subgroup <- aggregated.data |> filter(race.groups==subgroup)
  }

  # calculate measures
  measures.CI <- measures.with.CI(predictions=data.subgroup$model.score,
                                  labels=data.subgroup$sep3.outcome,
                                  nboot=nboot,
                                  min.pred=0)
  # make a figure of measures
  create.figure(measures=measures.CI, 
                threshold=threshold,
                predictions=data.subgroup$model.score,
                min.pred=0,
                title=subgroup,
                file.path=paste0(file.path, "Results/"))

  # calibration plot
  calibrated.predictions <- glm(sep3.outcome ~ model.score,
                                  data=data.subgroup,
                                  family="binomial")
  lp <- predict(calibrated.predictions)
  grDevices::png(file=paste0(file.path, "Results/Calibration ", subgroup, ".png"),
                 width=500, height=500, units="px")
  val.prob <- PredictionTools::val.prob.mi(lp.mi=lp, 
                                           lim=c(0, 0.2),
                                           y=data.subgroup$sep3.outcome,
                                           g=ifelse(subgroup=="All", 10, 5), 
                                           main=subgroup, dist=TRUE)
  dev.off()

  # add to table
  measures.table <- cbind(measures.table,
                          c(print.measures(threshold=threshold,
                                         predictions=data.subgroup$model.score,
                                         labels=data.subgroup$sep3.outcome),
                            sprintf("%.3f", val.prob$cindex),
                            sprintf("%.3f", val.prob$e.avg),
                            sprintf("%.3f", val.prob$e.90)))
}
colnames(measures.table) <- c("", subgroups)
openxlsx::write.xlsx(measures.table, file=paste0(file.path, "Results/Measures Table.xlsx"))
