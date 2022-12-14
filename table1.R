# Table 1
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Library
library(tidyverse)
library(tableone)
library(kableExtra)

# Data
df <- readRDS("data/demographics_BEAM.RDS")
officebp <- readRDS("data/officebp_summary.RDS")
officebp_v2 <- officebp %>% filter(visit == "V2") %>% select(-visit)
bia <- readRDS("data/bia_data.RDS")
diet <- readRDS("data/diet_summary.RDS")
bmi <- bia %>% filter(visit == "V2") %>% select(ID, BMI)
lab <- readRDS("data/lab_results.RDS")
lab_v2 <- lab %>% filter(visit == "V2") %>% select(ID, GFR, TC, HDL, LDL, TG)
nexfin <- readRDS("data/nexfin_data.RDS")
nexfin_v2 <- nexfin %>% filter(visit == "V2") %>% select(ID, meanBRS, SDNN)
df_total <- right_join(df, officebp_v2, by = "ID")
df_total <- right_join(df_total, bmi, by = "ID")
df_total <- right_join(df_total, lab_v2, by = "ID")
df_total <- right_join(df_total, nexfin_v2, by = "ID")
df_total <- right_join(df_total, diet, by = "ID")
head(df_total)
names(df_total)
df_total <- df_total %>% filter(! ID %in% c("BEAM_664", "BEAM_713", "BEAM_299"))

# Table 1
table1 <- df_total %>% 
    select(Age, Sex, BMI, Smoking, 
           # VG? CVD?
           BPlowMed, 
           # BPlowMed specs
           Systolic, Diastolic, Pulse, # this is baseline office BP
           # V1_Systolic, V1_Diastolic, V1_Pulse, # this is screening BP
           # Nexfin
           GFR, TC, HDL, LDL, TG, 
           Energy, Fibers,
           Treatment_group) %>% 
    CreateTableOne(data=., strata = 'Treatment_group', test = TRUE, addOverall = TRUE) %>% 
    print()
write.csv2(table1, "results/table1.csv")

