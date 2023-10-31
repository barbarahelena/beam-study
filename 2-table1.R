# Table 1, compliance and adverse events
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Library
library(tidyverse)
library(tableone)
library(kableExtra)

# Data
df <- readRDS("data/demographics_BEAM.RDS")
officebp <- readRDS("data/officebp_summary.RDS") %>% filter(visit == "V2")
bmi <- readRDS("data/bia_data.RDS") %>% filter(visit == "V2") %>% select(ID, BMI)
diet <- readRDS("data/diet_summary.RDS") %>% filter(visit == "V2")
lab <- readRDS("data/lab_results.RDS") %>% filter(visit == "V2") %>% select(ID, GFR, TC, HDL, LDL, TG)
abpm <- readRDS("data/abpm_total.RDS") %>% filter(visit == "V2")
df_total <- right_join(officebp, df, by = "ID") %>% 
                right_join(bmi, ., by = "ID") %>% 
                right_join(lab, ., by = "ID") %>% 
                right_join(diet, ., by = "ID") %>% 
                right_join(abpm, ., by = "ID") %>% 
                filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
                ungroup()

# Table 1
table1 <- df_total %>%
    select(Age, Sex, BMI, Smoking, 
           BPlowMed, 
           Vasomed_Type,
           Systolic, Diastolic, Pulse, 
           GFR, TC, HDL, LDL, TG, 
           Energy, Fibers, Alcohol,
           Treatment_group) %>% 
    CreateTableOne(data=.,
                   strata = 'Treatment_group', 
                   test = TRUE, 
                   addOverall = TRUE) %>% 
    print(nonnormal = "Alcohol")
write.csv2(table1, "results/table1.csv")

## Compliance
df_med <- readRDS("data/compliance_incl_pharmacy.RDS")
sum <- df_med %>% group_by(Group) %>% summarise(mean = mean(Perc_pills_taken),
                                                sd = sd(Perc_pills_taken),
                                                median = median(Perc_pills_taken),
                                                lowestq = quantile(Perc_pills_taken, 1/4),
                                                highestq = quantile(Perc_pills_taken, 3/4))
sum

## Adverse events
df_ae <- readRDS("data/adverse_events.RDS")
df_ae
