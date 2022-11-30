## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)

# Functions
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.line.y = element_line(colour="black"),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 

# Data
df <- rio::import("data/BEAM_export_20221124.csv")
groups <- rio::import("data/treatment_groups.xlsx")
## repeated measurements to be inserted here (AE and medication)

# Change participant Id into ID
df <- df %>% select(ID = `Participant Id`, everything(.), -`V829`)
names(df)
names(groups)

# Merge with treatment allocation df and set NA values
df <- right_join(df, groups, by = "ID") %>% 
    naniar::replace_with_na_all(condition = ~(.x == -99 | .x == -96)) %>% 
    mutate(across(where(is.character), ~ na_if(.,"")))
Amelia::missmap(df)
names(df)

# Separate into different dataframes per datasource
demographics <- df %>% select(ID, Sex, Age, AgeStrata, Smoking, PackYears,
                              AlcoholUse, AlcoholUnits, History_List,
                              eGFR, Date_eGFR, BPlowMed = Medication_BPlowering,
                              HT_years = History_Hypertension_Years,
                              V1_DateTime, V2_DateTime, V3_DateTime, V4_DateTime, V5_datetime) # check
homebp <- df %>% select(ID, contains("HomeBP")) # check
diet <- df %>% select(ID, contains("Diet")) # check
bp_measurement <- df %>% select(ID, contains("BP_Measurement"))
bia <- df %>% select(ID, contains("BIA"))
nexfin <- df %>% select(ID, contains("Nexfin"))
abpm <- df %>% select(ID, contains("ABPM")) # check
pbmc <- df %>% select(ID, contains("PBMC"))
lab <- df %>% select(ID, contains("Lab_"),
                            contains("Leukodiff"),
                            contains("UrineLab"))
samplestorage <- df %>% select(ID, contains("Cryovials"),
                               contains("Proc"),
                               contains("Feces"))

#### Demographics clean and check ####
names(demographics)
str(demographics)

demographics <- demographics %>% 
    mutate(BPlowMed = case_when(
        BPlowMed == 1 ~ paste0("Yes"),
        BPlowMed == 0 ~ paste0("No")
    ),
    AlcoholUse = case_when(
        AlcoholUse == 1 ~ paste0("Yes"),
        AlcoholUse == 0 ~ paste0("No"))) %>% 
    mutate(across(c("Sex", "Smoking", "BPlowMed", "AlcoholUse"), as.factor))
    

str(demographics)
plot(demographics)
graphics.off()
hist(demographics$Age)
plot(demographics$Sex)
plot(demographics$HT_years)
any(demographics$AlcoholUse == "Yes" & demographics$AlcoholUnits == 0)
any(demographics$AlcoholUse == "No" & demographics$AlcoholUnits != 0)
any(demographics$BPlowMed == "Yes" & is.na(demographics$HT_years))

#### Home BP clean and check ####
names(homebp)
Amelia::missmap(homebp)
stripnames_homebp <- function(varname) {
    varname <- str_remove(varname, "_min_")
    varname <- str_remove(varname, "_mmHg_")
    varname <- str_remove(varname, "HomeBP_")
    varname <- str_remove(varname, "Measurement_")
    varname <- str_replace(varname, "_BP", "BP")
    return(varname)
}

# Strip/change names for parts defined above; get rid of vars that we don't need
homebp <- homebp %>% rename_with(., stripnames_homebp) %>% 
    select(!(contains("Remarks") | contains("_1_"))) %>% 
    mutate(across(contains("Pulse") | contains("BP"), as.numeric)) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
names(homebp)
str(homebp)

# min needs to pivot separately (has no timing); dates need to pivot separately
homebp_dates <- homebp %>% select(ID, contains("Date"))
homebp_nodates <- homebp %>% select(ID, !contains("Date"), -contains("min")) 
homebp_min <- homebp %>% select(ID, contains("min")) %>% select(-contains("Date"))

homebp_long_min <- pivot_longer(homebp_min, cols = 2:ncol(homebp_min),
                            names_to = c("visit", "week", "measurement", "item"),
                            names_sep = "_",
                            values_to = "value")
homebp_long <- pivot_longer(homebp_nodates, cols = 2:ncol(homebp_nodates),
                            names_to = c("visit", "week", "measurement", "timing", "item"),
                            names_sep = "_",
                            values_to = "value")
names(homebp_long)
head(homebp_long)
head(homebp_long_min)
homebp_long <- bind_rows(homebp_long, homebp_long_min)
head(homebp_long)
homebp_wider <- pivot_wider(homebp_long, id_cols = c(1:5), values_from = value,
                            names_from = item)
head(homebp_wider)

head(homebp_dates)
homebp_long_dates <- homebp_dates %>% 
    rename_with(., ~ str_remove_all(.x, "_Date")) %>% 
    pivot_longer(homebp_dates, cols = 2:ncol(homebp_dates),
                 names_to = c("visit", "week", "measurement", "timing"),
                 names_sep = "_",
                 values_to = "date") %>% 
    select(-measurement, -visit) %>% 
    distinct()
head(homebp_long_dates)

homebp_summary_pertime <- homebp_wider %>% group_by(ID, week, timing) %>% 
    summarise(mean_sbp = mean(SystolicBP, na.rm = TRUE), 
              mean_dbp = mean(DiastolicBP, na.rm = TRUE),
              mean_pulse = mean(Pulse, na.rm = TRUE))
head(homebp_summary_pertime)
homebp_summary_perday <- homebp_summary_pertime %>% group_by(ID, week) %>% 
    summarise(mean_sbp = mean(mean_sbp, na.rm = TRUE), 
              mean_dbp = mean(mean_dbp, na.rm = TRUE),
              mean_pulse = mean(mean_pulse, na.rm = TRUE)) %>% 
    ungroup()
head(homebp_summary_perday)

homebp_long_dates_morning <- homebp_long_dates %>% filter(timing == "morning" | is.na(timing)) %>% 
    select(-timing)
head(homebp_long_dates_morning)

homebp_total_pertime <- right_join(homebp_summary_pertime, homebp_long_dates, 
                           by = c("ID", "week", "timing")) %>% 
                        mutate(date = lubridate::dmy(date),
                               timing = as.factor(timing),
                               week = case_when(
                                   str_detect(week, "min") ~ str_replace(week, "min", "-"),
                                   str_detect(week, "week") ~ str_replace(week, "week", "")
                               ),
                               week = as.numeric(week)) %>% 
                        arrange(ID, date)
head(homebp_total_pertime)

homebp_total_perweek <- right_join(homebp_summary_perday, 
                                   homebp_long_dates_morning, 
                                   by = c("ID", "week")) %>% 
                        mutate(date = lubridate::dmy(date),
                               week = case_when(
                                   str_detect(week, "min") ~ str_replace(week, "min", "-"),
                                   str_detect(week, "week") ~ str_replace(week, "week", "")
                               ),
                               week = as.numeric(week)) %>% 
                        arrange(ID, date) 
head(homebp_total_perweek)

#### ABPM clean and check ####
names(abpm)
stripnames_abpm <- function(varname) {
    varname <- str_remove(varname, "_min_")
    varname <- str_remove(varname, "_mmHg_")
    varname <- str_remove(varname, "ABPM_Measurements_")
    return(varname)
}
abpm <- abpm %>% select(ID, contains("TotalNo"), contains("SuccessNo"),
                        contains("Arm"), contains("Bedtime"), contains("WakeUpTime"),
                        contains("Measurements")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
str(abpm)
abpm <- abpm %>% rename_with(., stripnames_abpm, contains("Measurements"))
names(abpm)

# no of successful measurement should not exceed total
any(abpm$V2_ABPM_TotalNo < abpm$V2_ABPM_SuccessNo) 
any(abpm$V4_ABPM_TotalNo < abpm$V4_ABPM_SuccessNo)
any(abpm$V5_ABPM_TotalNo < abpm$V5_ABPM_SuccessNo)

# calculation perc of successful measurements
abpm <- abpm %>% mutate(
    V2_Success_perc = (V2_ABPM_SuccessNo / V2_ABPM_TotalNo) *100,
    V4_Success_perc = (V4_ABPM_SuccessNo / V4_ABPM_TotalNo) *100,
    V5_Success_perc = (V5_ABPM_SuccessNo / V5_ABPM_TotalNo) *100
)

str(abpm)

# pivot longer in two parts - short and long variables
abpm_longvars <- abpm %>% select(ID, contains("SD"), contains("Mean"))
abpm_shortvars <- abpm %>% select(ID, contains("ABPM")) %>% mutate(across(everything(.), as.character))

abpm_longvars <- pivot_longer(abpm_longvars, cols = 2:ncol(abpm_longvars),
                                   names_to = c("visit", "timing", "statistic", "outcome"),
                                   names_sep = "_",
                                   values_to = "value") %>% 
                        pivot_wider(., id_cols = 1:2,
                                    values_from = value, names_from = c(timing, outcome, statistic),
                                    names_sep = "_")
head(abpm_longvars)
abpm_shortvars <- pivot_longer(abpm_shortvars, cols = 2:ncol(abpm_shortvars),
                                   names_to = c("visit", "abpm", "variable"),
                                   names_sep = "_",
                                   values_to = "value") %>% 
                        select(-abpm) %>% 
                        pivot_wider(., id_cols = 1:2, values_from = value,
                                    names_from = variable)

abpm_total <- right_join(abpm_longvars, abpm_shortvars, by = c("ID", "visit")) %>% 
                mutate(across(contains("SD") | contains("Mean") | c("TotalNo", "SuccessNo"),
                                         as.numeric),
                   across(c("Arm", "visit"), as.factor),
                   across(c("Bedtime", "WakeUpTime"), ~lubridate::hm(.x))) %>% 
                droplevels(.)
str(abpm_total)

#### Dietary data ####
str(diet)
diet <- diet %>% select(-(contains("File") | contains("Remarks"))) %>% 
    rename_with(., ~str_replace_all(., "Saturated_fat", "SaturatedFat"))
str(diet)
diet_long <- diet %>% 
    pivot_longer(., cols = 2:ncol(diet), names_sep = "_",
                 names_to = c("visit", "drop1", "variable", "drop2", "day"),
                 values_to = "value") %>% 
    select(-drop1, -drop2)
head(diet_long)



#### Nexfin ####


#### Office BP ####