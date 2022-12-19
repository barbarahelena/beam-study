## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(lubridate)

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
fecalscfa <- rio::import("data/221201_Fecal_SCFA_tidy.xlsx") %>%
    select(ID = Studienummer, visit = Visit, contains("DW"), contains("WW"))
reninaldo <- readxl::read_xlsx("data/15122022_ReninAldo_ErasmusMC.xlsx", skip = 3)
## repeated measurements to be inserted here (AE and medication)

# Change participant Id into ID
df <- df %>% select(ID = `Participant Id`, everything(.), -`V829`)
names(df)
names(groups)

# Merge with treatment allocation df and set NA values
groups <- groups %>%  mutate(Treatment_group = na_if(Treatment_group, "-"))
df <- right_join(df, groups, by = "ID") %>% 
    naniar::replace_with_na_all(condition = ~(.x == -99 | .x == -96 | .x == -95 |
                                                  .x == -98)) %>% 
    mutate(across(where(is.character), ~ na_if(.,"")))
Amelia::missmap(df)
names(df)

# Separate into different dataframes per datasource
demographics <- df %>% select(ID, Sex, Age, AgeStrata, Smoking, PackYears,
                              AlcoholUse, AlcoholUnits, History_List,
                              eGFR, Date_eGFR, BPlowMed = Medication_BPlowering,
                              HT_years = History_Hypertension_Years, V1_Weight, V1_Height,
                              V1_DateTime, V2_DateTime, V3_DateTime, V4_DateTime, V5_datetime,
                              Treatment_group, (contains("V1") & contains("BP"))) 
homebp <- df %>% select(ID, contains("HomeBP")) 
diet <- df %>% select(ID, contains("Diet")) 
bp_measurement <- df %>% select(ID, contains("BP_Measurement")) 
bia <- df %>% select(ID, contains("BIA")) 
nexfin <- df %>% select(ID, contains("Nexfin")) 
abpm <- df %>% select(ID, contains("ABPM")) 
pbmc <- df %>% select(ID, contains("PBMC")) 
lab <- df %>% select(ID, contains("Lab_"),
                            contains("Leukodiff"),
                            contains("UrineLab")) 
samplestorage <- df %>% select(ID, contains("Cryovials"),
                               contains("Proc"),
                               contains("Feces"),
                               contains("Urine_"),
                               contains("Blood")) 
intervention <- df %>% select(ID, contains("Capsules"))

#### Demographics clean and check ####
names(demographics)
str(demographics)

bp_screening <- demographics %>% select(ID, contains("BP_"), 
                            -(contains("Arm") | contains("Remarks") | contains("_1_"))) %>% 
    rename_with(., ~ str_remove(.x, "BP_Pressure_Measurement_")) %>% 
    rename_with(., ~ str_remove(.x, "_mmHg_")) %>% 
    rename_with(., ~ str_remove(.x, "_min_")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "measurement", "variable"),
                 values_to = "value", names_sep = "_") %>% 
    pivot_wider(., id_cols = 1:3, names_from = "variable", values_from = "value") %>% 
    group_by(ID, visit) %>% 
    summarise(across(c("Systolic", "Diastolic", "Pulse"), ~ mean(.x)))

demographics <- demographics %>% 
    select(-contains("BP_")) %>% 
    mutate(BPlowMed = case_when(
        BPlowMed == 1 ~ paste0("Yes"),
        BPlowMed == 0 ~ paste0("No")
    ),
    AlcoholUse = case_when(
        AlcoholUse == 1 ~ paste0("Yes"),
        AlcoholUse == 0 ~ paste0("No"))) %>% 
    mutate(across(c("Sex", "Smoking", "BPlowMed", "AlcoholUse", "Treatment_group"), 
                  as.factor)) %>% 
    right_join(., bp_screening %>% select(-visit), by = "ID") %>%
    rename(V1_Systolic = Systolic, V1_Diastolic = Diastolic, V1_Pulse = Pulse)
    

str(demographics)
plot(demographics)
graphics.off()
hist(demographics$Age)
plot(demographics$Sex)
plot(demographics$HT_years)
any(demographics$AlcoholUse == "Yes" & demographics$AlcoholUnits == 0)
any(demographics$AlcoholUse == "No" & demographics$AlcoholUnits != 0)
any(demographics$BPlowMed == "Yes" & is.na(demographics$HT_years))

saveRDS(demographics, "data/demographics_BEAM.RDS")
write.csv2(demographics, "data/demographics_BEAM.csv")

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

# Found typo in data; remove this upon next data export
homebp$V3_week1_2_morning_DiastolicBP
homebp$V3_week1_2_morning_DiastolicBP[which(homebp$V3_week1_2_morning_DiastolicBP == 960)] <- 96
homebp$V3_week1_2_morning_DiastolicBP

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

homebp_long <- bind_rows(homebp_long, homebp_long_min)
homebp_wider <- pivot_wider(homebp_long, id_cols = c(1:5), values_from = value,
                            names_from = item)
head(homebp_wider)

names(homebp_dates)
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

saveRDS(homebp_total_pertime, "data/homebp_pertime_BEAM.RDS")
write.csv2(homebp_total_pertime, "data/homebp_pertime_BEAM.csv")

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

saveRDS(homebp_total_perweek, "data/homebp_perweek_BEAM.RDS")
write.csv2(homebp_total_perweek, "data/homebp_perweek_BEAM.csv")

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
plot(abpm_total %>% select(contains("systolic")))
plot(abpm_total %>% select(contains("diastolic")))
plot(abpm_total %>% select(contains("HR")))

saveRDS(abpm_total, "data/abpm_total.RDS")
write.csv2(abpm_total, "data/abpm_total.csv")

#### Dietary data ####
str(diet)
diet <- diet %>% select(-(contains("File") | contains("Remarks"))) %>% 
    rename_with(., ~str_replace_all(., "Saturated_fat", "SaturatedFat")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
str(diet)
diet_long <- diet %>% 
    pivot_longer(., cols = 2:ncol(diet), names_sep = "_",
                 names_to = c("visit", "drop1", "variable", "drop2", "day"),
                 values_to = "value") %>% 
    select(-drop1, -drop2) %>% 
    pivot_wider(., id_cols = c("ID", "visit", "day"), names_from = "variable",
                values_from = "value") %>% 
    select(-Date) %>% 
    mutate(visit = as.factor(visit))
head(diet_long)

# found typo in data; corrected in Castor, so remove upon next export
diet_long$Salt[which(diet_long$Salt == 1104)] <- 11.04

print(diet_long[which(as.numeric(diet_long$Energy) < 1000),], n = 100)

diet_summary <- diet_long %>% mutate(across(.cols = c(3:11), as.numeric)) %>% 
    group_by(ID, visit) %>% 
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
    select(-day)

head(diet_summary)
plot(diet_summary[,3:ncol(diet_summary)])

saveRDS(diet_summary, "data/diet_summary.RDS")
write.csv2(diet_summary, "data/diet_summary.csv")

#### Nexfin ####
str(nexfin)
# correct wrong value (changed in castor, remove upon next export)
nexfin$V2_Nexfin_MAP[which(nexfin$V2_Nexfin_MAP == "93.9262093.92620137137")] <-
    "93.92620137"

# select vars
nexfin <- nexfin %>% select(-contains("Remark"), -contains("yesno"), -contains("File"),
                            -contains('hand'), -contains("reason")) %>% 
    mutate(V2_Nexfin_MAP = as.numeric(V2_Nexfin_MAP)) %>% 
    rename_with(., ~ str_remove(.x, "_Nexfin")) %>% 
    rename_with(., ~ str_replace(.x, "N_physiocals", "physiocal"))
names(nexfin)

nexfin_long <- nexfin %>% pivot_longer(., cols = c(2:ncol(nexfin)), 
                                       names_to = c("visit", "variable"),
                                       names_sep = "_",
                                       values_to = "value") %>% 
                    pivot_wider(., id_cols = c(1,2), names_from = "variable",
                                values_from = "value")

# two outliers below have been corrected in Castor; remove upon next export
nexfin_long$nLargestStableBeats[which(nexfin_long$nLargestStableBeats == 503)] <- 70
nexfin_long$SVR[which(nexfin_long$SVR == 13)] <- 1359.5487804878

nexfin_long$ID[which(nexfin_long$nLargestStableBeats < 30)]
head(nexfin_long)
dim(nexfin_long)
plot(nexfin_long[,3:15])
plot(nexfin_long[,16:30])

saveRDS(nexfin_long, "data/nexfin_data.RDS")
write.csv2(nexfin_long, "data/nexfin_data.csv")

#### Office BP ####
names(bp_measurement)
stripnames_officebp <- function(varname) {
    varname <- str_remove(varname, "_min_")
    varname <- str_remove(varname, "_mmHg_")
    varname <- str_remove(varname, "OfficeBP_")
    varname <- str_remove(varname, "Measurement_")
    varname <- str_remove(varname, "BP_")
    return(varname)
}

# Strip/change names for parts defined above; get rid of vars that we don't need
officebp <- bp_measurement %>% rename_with(., stripnames_officebp) %>% 
    select(-contains("_1_")) %>% 
    mutate(across(contains("Pulse") | contains("Systolic") | contains("Diastolic"), as.numeric)) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
names(officebp)
str(officebp)

officebp_long <- pivot_longer(officebp, cols = 2:ncol(officebp), 
                              names_to = c("visit", "measurement", "variable"),
                              names_sep = "_",
                              values_to = "value") %>% 
    pivot_wider(., id_cols = 1:3, names_from = "variable", values_from = "value")
head(officebp_long)

officebp_summary <- officebp_long %>% group_by(ID, visit) %>% 
    summarise(across(c("Systolic", "Diastolic", "Pulse"), mean))
plot(officebp_summary[,3:ncol(officebp_summary)])

saveRDS(officebp_summary, "data/officebp_summary.RDS")
write.csv2(officebp_summary, "data/officebp_summary.csv")

#### PBMC ####
str(pbmc)
pbmc %>% select(contains("Reason"))
pbmc_sel <- pbmc %>% select(ID, contains("Cellcount"), contains("Cryovials")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
pbmc_sel

pbmc_long <- pivot_longer(pbmc_sel, cols = 2:ncol(pbmc_sel), names_sep = "_",
                          names_to = c("visit", "drop", "variable"), values_to = "value") %>% 
            select(-drop) %>% 
            pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value") %>% 
            select(-Cryovials) %>% 
            mutate(Residual_count = Cellcount -3)
head(pbmc_long)

nrow(pbmc_long[which(pbmc_long$Residual_count < 15),]) / nrow(pbmc_long)

gghistogram(pbmc_long$Residual_count, bins = 15, fill = pal_jama()(1)) + theme_Publication() +
    xlab("Cell count vial #4") +
    ggtitle("PBMC cell count vial 4") +
    geom_vline(xintercept = 15, color = pal_jama()(2)[2], size = 1.0) +
    scale_x_continuous(n.breaks = 6)
ggsave("results/pbmc/cellcountvial4.pdf", width = 5, height = 5)

saveRDS(pbmc_long, "data/pbmc_storage.RDS")
write.csv2(pbmc_long, "data/pbmc_storage.csv")

#### Lab ####
str(lab)
print(lab %>% select(contains("Remarks")), n=24)
lab_sel <- lab %>% select(ID, !contains("Remarks")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    rename_with(., ~ str_remove(.x, "Lab_")) %>% 
    rename_with(., ~ str_remove(.x, "Leukodiff_"))

names(lab_sel)

lab_long <- lab_sel %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "variable"),
                 names_sep = "_", values_to = "value") %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
head(lab_long)    
lab_compl <- right_join(lab_long, fecalscfa, by = c("ID", "visit"))
plot(lab_long[,3:11])
plot(lab_long[,12:15])
plot(lab_long[,16:18])
plot(lab_long[,19:24])

saveRDS(lab_compl, "data/lab_results.RDS")
write.csv2(lab_compl, "data/lab_results.csv")

#### BIA ####
str(bia)
bia <- bia %>% select(ID, !contains("Remarks")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    rename_with(., ~ str_remove(.x, "_BIA"))
str(bia)

bia_long <- bia %>% pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "variables"),
                                 names_sep = "_", values_to = "value") %>% 
                    pivot_wider(., id_cols = 1:2, names_from = "variables", values_from = "value") %>% 
                    mutate(visit = as.factor(visit))
str(bia_long)
plot(bia_long[,3:13])

saveRDS(bia_long, "data/bia_data.RDS")
write.csv2(bia_long, "data/bia_data.csv")

#### Sample storage ####
str(samplestorage)
names(samplestorage)
sampleremarks <- samplestorage %>% select(ID, contains("Remarks"))
samplestorage <- samplestorage %>% select(!(contains("PBMC") | contains("Check"))) %>% 
    rename_with(., ~ str_remove(.x, "Cryovials_")) %>% 
    rename_with(., ~ str_remove(.x, "_Number_of_cryovials")) %>% 
    rename_with(., ~ str_replace(.x, "Lithium_heparin", "heparin")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))

# blood samples
print(samplestorage %>% select(ID, contains("Blood") & contains("SOP")), n = 21)
blood <- samplestorage %>% select(ID, contains("Blood") & !contains("SOP")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
str(blood)

write.csv2(blood, "data/bloodsamples.csv")

# fecal samples; date pivoting separately
feces_date <- samplestorage %>% select(ID, contains("Feces_DateTime")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
feces <- samplestorage %>% select(ID, contains("Feces"), -(contains("FecesCollection")), 
                                                         -(contains("Feces_DateTime"))) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value") %>% 
    rename(Samples_stored = Storage)
feces <- right_join(feces, feces_date, by = c("ID", "visit"))
head(feces)

write.csv2(feces, "data/fecessamples.csv")

# urine samples; date pivoting separately
print(samplestorage %>% select(ID, contains("Urine") & contains("SOP")), n = 21)
urine_dates <- samplestorage %>% select(ID, contains("Urine") & contains("DateTime")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value") %>% 
    mutate(across(contains("Date"), ~dmy_hm(.x))) %>% 
    mutate(CollectionTime = EndDateTime - StartDateTime)
urine <- samplestorage %>% select(ID, contains("Urine") & !contains("DateTime") &
                                      !contains("SOP")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
urine <- right_join(urine, urine_dates, by = c("ID", "visit"))
print(urine, n = 63)

saveRDS(urine, "data/urinesamples.RDS")
write.csv2(urine, "data/urinesamples.csv")

#### Renin aldo ####
reninaldo <- reninaldo %>% select(Sample_ID, ID = Subject_ID, Renin = Kolom1, 
                                  Aldosterone = Kolom2) %>% 
    mutate(visit = str_extract(Sample_ID, "V[0-9]")) %>% 
    filter(str_detect(ID, "BEAM"))
saveRDS(reninaldo, "data/reninaldo.RDS")
write.csv2(reninaldo, "data/reninaldo.csv")

#### Intervention ####
str(intervention)
print(intervention %>% select(ID, contains("Remarks")))
intervention <- intervention %>% select(ID, !(contains("Checks") | contains("Remarks"))) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Total_NonCompliantDays = V3_Capsules_NonCompliantDays + V4_Capsules_NonCompliantDays)
intervention$Total_NonCompliantDays
names(intervention)

saveRDS(intervention, "data/intervention.RDS")
write.csv2(intervention, "data/intervention.csv")

