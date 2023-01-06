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
df_total <- right_join(officebp_v2, df, by = "ID")
df_total <- right_join(bmi, df_total, by = "ID")
df_total <- right_join(lab_v2, df_total, by = "ID")
df_total <- right_join(nexfin_v2, df_total, by = "ID")
df_total <- right_join(diet, df_total, by = "ID")
head(df_total)
names(df_total)
df_total <- df_total %>% filter(! ID %in% c("BEAM_664")) %>% ungroup()

# Table 1
table1 <- df_total %>%
    filter(! visit %in% c("V4", "V5")) %>% 
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

tab <- df %>% select(ID, Treatment_group)
write.csv2(tab, "data/randomisation.csv")
