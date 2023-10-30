## Data cleaning
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(lubridate)

# Data
df <- rio::import("data/BEAM_export_20221219.csv")
fecalscfa <- rio::import("data/221201_Fecal_SCFA_tidy.xlsx") %>%
    dplyr::select(ID = Studienummer, visit = Visit, contains("DW"), contains("WW"))
plasmascfa <- rio::import("data/221128_EDTASamples_BEAM_Plasma_SCFA report_Xinmin Li-02-08-2023.xlsx") %>% select(ID = Subject_ID, visit = Sample_ID, 8:13)
reninaldo <- readxl::read_xlsx("data/15122022_ReninAldo_ErasmusMC.xlsx", skip = 3)
alcoholuse <- rio::import("data/230215_AlcoholUse.xlsx")
serotonin <- rio::import('data/230421_BEAM_serotonin_EDTA_tidy.xlsx')
calprotectin <- rio::import("data/230324_BEAM_calprotectine.xlsx")
samplelijst_calprotectin <- rio::import("data/230320_Samplelijst_FecesvoorCalprotectin.xlsx")
inflelisa <- rio::import("data/infl_elisa_tidy_2.xlsx")
vasomed <- read_csv2("data/BEAM_csv_export_20221219104054/BEAM_Vasoactive_medication_export_20221219.csv") %>% 
    select(ID = `Participant Id`, everything()) %>% filter(`Participant Status` == "Completed") 
adverse_events <- read_csv2("data/BEAM_Adverse_event_export_20221219.csv") %>% 
    select(ID = `Participant Id`, AE_Description)

# Change participant Id into ID
df <- df %>% dplyr::select(ID = `Participant Id`, everything(.), -`V829`)
names(df)

# Set NA values, plot missing data with Amelia missmap
df <- df %>% naniar::replace_with_na_all(condition = ~(.x == -99 | .x == -96 | .x == -95 |
                                                  .x == -98)) %>% 
    mutate(across(where(is.character), ~ na_if(.,"")))
Amelia::missmap(df)
names(df)

# Separate into different dataframes per data source
demographics <- df %>% dplyr::select(ID, Sex, Age, AgeStrata, Smoking, PackYears,
                              AlcoholUse, AlcoholUnits, History_List,
                              eGFR, Date_eGFR, BPlowMed = Medication_BPlowering,
                              HT_years = History_Hypertension_Years, V1_Weight, V1_Height,
                              V1_DateTime, V2_DateTime, V3_DateTime, V4_DateTime, V5_datetime,
                              Treatment_group = `Randomization Group`, (contains("V1") & contains("BP"))) 
homebp <- df %>% dplyr::select(ID, contains("HomeBP")) 
diet <- df %>% dplyr::select(ID, contains("Diet")) 
alcohol <- alcoholuse %>% dplyr::select(ID = Subject, 2:4)
bp_measurement <- df %>% dplyr::select(ID, contains("BP_Measurement")) 
bia <- df %>% dplyr::select(ID, contains("BIA"), contains("Weight")) 
nexfin <- df %>% dplyr::select(ID, contains("Nexfin")) 
abpm <- df %>% dplyr::select(ID, contains("ABPM")) 
pbmc <- df %>% dplyr::select(ID, contains("PBMC")) 
lab <- df %>% dplyr::select(ID, contains("Lab_"),
                            contains("Leukodiff"),
                            contains("UrineLab")) 
samplestorage <- df %>% dplyr::select(ID, contains("Cryovials"),
                               contains("Proc"),
                               contains("Feces"),
                               contains("Urine_"),
                               contains("Blood")) 
intervention <- df %>% dplyr::select(ID, contains("Capsules"))

#### Demographics clean and check ####
names(demographics)
str(demographics)

bp_screening <- demographics %>% dplyr::select(ID, contains("BP_"), 
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
    dplyr::select(-contains("BP_")) %>% 
    mutate(BPlowMed = case_when(
        BPlowMed == 1 ~ paste0("Yes"),
        BPlowMed == 0 ~ paste0("No")
    ),
    AlcoholUse = case_when(
        AlcoholUse == 1 ~ paste0("Yes"),
        AlcoholUse == 0 ~ paste0("No"))) %>% 
    mutate(across(c("Sex", "Smoking", "BPlowMed", "AlcoholUse", "Treatment_group"), 
                  as.factor),
           Treatment_group = fct_relevel(Treatment_group, "Butyrate", after = 1L)) %>% 
    right_join(., bp_screening %>% dplyr::select(-visit), by = "ID") %>%
    full_join(., vasomed, by = "ID") %>% 
    rename(V1_Systolic = Systolic, V1_Diastolic = Diastolic, V1_Pulse = Pulse) %>% 
    mutate(V2_time = hms::as_hms(lubridate::dmy_hm(V2_DateTime)),
            V4_time = hms::as_hms(lubridate::dmy_hm(V4_DateTime)),
            V5_time = hms::as_hms(lubridate::dmy_hm(V5_datetime)),
            V4_time_bin = case_when(
                V4_time > hms::as_hms("09:00:00") ~ paste("late"),
                V4_time <= hms::as_hms("09:00:00") ~ paste("early")
            ),
            V4_time_bin = as.factor(V4_time_bin),
            V4_hourdiff = (V4_time - hms::as_hms("07:30:00"))/3600,
           Vasomed_Type = fct_recode(Vasomed_Type, "Diuretics" = "DiU",
                                      "ARB" = "ARB",
                                      "CaAnt" = "CAD")
           )

# checks
any(demographics$AlcoholUse == "Yes" & demographics$AlcoholUnits == 0)
any(demographics$AlcoholUse == "No" & demographics$AlcoholUnits != 0)
any(demographics$BPlowMed == "Yes" & is.na(demographics$HT_years))

# save dataset
saveRDS(demographics, "data/demographics_BEAM.RDS")
write.csv2(demographics, "data/demographics_BEAM.csv")


#### ABPM clean and check ####
names(abpm)
stripnames_abpm <- function(varname) {
    varname <- str_remove(varname, "_min_")
    varname <- str_remove(varname, "_mmHg_")
    varname <- str_remove(varname, "ABPM_Measurements_")
    return(varname)
}
abpm <- abpm %>% dplyr::select(ID, contains("TotalNo"), contains("SuccessNo"),
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
abpm_longvars <- abpm %>% dplyr::select(ID, contains("SD"), contains("Mean"))
abpm_shortvars <- abpm %>% dplyr::select(ID, contains("ABPM")) %>% mutate(across(everything(.), as.character))

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
                        dplyr::select(-abpm) %>% 
                        pivot_wider(., id_cols = 1:2, values_from = value,
                                    names_from = variable)

abpm_total <- right_join(abpm_longvars, abpm_shortvars, by = c("ID", "visit")) %>% 
                mutate(across(contains("SD") | contains("Mean") | c("TotalNo", "SuccessNo"),
                                         as.numeric),
                   across(c("Arm", "visit"), as.factor),
                   across(c("Bedtime", "WakeUpTime"), ~lubridate::hm(.x))) %>% 
                droplevels(.) %>% 
                mutate(Awake_MAP_Mean = (Awake_systolic_Mean+(2*Awake_diastolic_Mean))/3,
                       Awake_PP_Mean = (Awake_systolic_Mean-Awake_diastolic_Mean),
                       Asleep_MAP_Mean = (Asleep_systolic_Mean+(2*Asleep_diastolic_Mean))/3,
                       Asleep_PP_Mean = (Asleep_systolic_Mean-Asleep_diastolic_Mean),
                       Total_MAP_Mean = (Total_systolic_Mean+(2*Total_diastolic_Mean))/3,
                       Total_PP_Mean = (Total_systolic_Mean-Total_diastolic_Mean),
                       )
str(abpm_total)
plot(abpm_total %>% dplyr::select(contains("systolic")))
plot(abpm_total %>% dplyr::select(contains("diastolic")))
plot(abpm_total %>% dplyr::select(contains("HR")))

saveRDS(abpm_total, "data/abpm_total.RDS")
write.csv2(abpm_total, "data/abpm_total.csv")

#### Dietary data ####
str(diet)
diet <- diet %>% dplyr::select(-(contains("File") | contains("Remarks"))) %>% 
    rename_with(., ~str_replace_all(., "Saturated_fat", "SaturatedFat")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
str(diet)
diet_long <- diet %>% 
    pivot_longer(., cols = 2:ncol(diet), names_sep = "_",
                 names_to = c("visit", "drop1", "variable", "drop2", "day"),
                 values_to = "value") %>% 
    dplyr::select(-drop1, -drop2) %>% 
    pivot_wider(., id_cols = c("ID", "visit", "day"), names_from = "variable",
                values_from = "value") %>% 
    dplyr::select(-Date) %>% 
    mutate(visit = as.factor(visit))
head(diet_long)

print(diet_long[which(as.numeric(diet_long$Energy) < 1000),], n = 100)

diet_summary <- diet_long %>% mutate(across(.cols = c(3:11), as.numeric)) %>% 
    group_by(ID, visit) %>% 
    summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>% 
    dplyr::select(-day)

head(diet_summary)
plot(diet_summary[,3:ncol(diet_summary)])

alcohol_long <- alcohol %>% 
    pivot_longer(., cols = 2:4, names_sep = "_",
                 names_to = c("visit", "drop1"),
                 values_to = "Alcohol") %>% 
    dplyr::select(-drop1) 

diet_summary <- right_join(diet_summary, alcohol_long, by = c("ID", "visit"))
head(diet_summary)

saveRDS(diet_summary, "data/diet_summary.RDS")
write.csv2(diet_summary, "data/diet_summary.csv")

#### Nexfin ####
str(nexfin)

# dplyr::select vars
nexfin <- nexfin %>% dplyr::select(-contains("Remark"), -contains("yesno"), -contains("File"),
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
    dplyr::select(-contains("_1_")) %>% 
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

#### Lab ####
str(lab)
print(lab %>% dplyr::select(contains("Remarks")), n=24)
lab_sel <- lab %>% dplyr::select(ID, !contains("Remarks")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    rename_with(., ~ str_remove(.x, "Lab_")) %>% 
    rename_with(., ~ str_remove(.x, "Leukodiff_"))

names(lab_sel)

lab_long <- lab_sel %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "variable"),
                 names_sep = "_", values_to = "value") %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
head(lab_long)    

plasmascfa <- plasmascfa %>% mutate(visit = str_extract(visit, "V[0-9]"))
lab_compl <- right_join(lab_long, fecalscfa, by = c("ID", "visit")) %>% 
    right_join(., plasmascfa, by = c("ID", "visit"))
plot(lab_long[,3:11])
plot(lab_long[,12:15])
plot(lab_long[,16:18])
plot(lab_long[,19:24])

saveRDS(lab_compl, "data/lab_results.RDS")
write.csv2(lab_compl, "data/lab_results.csv")

serotonin <- serotonin %>% mutate(visit = str_c("V", Visit), Visit = NULL)
colnames(serotonin)[3] <- "serotonin_uM"
saveRDS(serotonin, "data/serotonin.RDS")
write.csv2(serotonin, "data/serotonin.csv")

samplelijst_calprotectin <- samplelijst_calprotectin %>% select(ID = Studienummer, visit = Visit, BARCODE = ID)
calprotectin <- left_join(calprotectin, samplelijst_calprotectin, by = "BARCODE") 
calprotectin <- calprotectin %>% select(ID, visit, calprotectin_ug_g = UITSLAG)
saveRDS(calprotectin, "data/calprotectin.RDS")
write.csv2(calprotectin, "data/calprotectin.csv")

reninaldo <- reninaldo %>% dplyr::select(Sample_ID, ID = Subject_ID, Renin = Kolom1, 
                                         Aldosterone = Kolom2) %>% 
    mutate(visit = str_extract(Sample_ID, "V[0-9]")) %>% 
    filter(str_detect(ID, "BEAM"))
saveRDS(reninaldo, "data/reninaldo.RDS")
write.csv2(reninaldo, "data/reninaldo.csv")

inflelisa <- inflelisa %>% select(SampleID = ID, everything(.)) %>% 
    mutate(visit = str_extract(SampleID, "V[0-9]"),
           ID = str_extract(SampleID, "BEAM_[0-9]*"))
saveRDS(inflelisa, "data/inflelisa.RDS")
write.csv2(inflelisa, "data/inflelisa.csv")

#### BIA ####
str(bia)
bia <- bia %>% dplyr::select(ID, !contains("Remarks")) %>% 
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
sampleremarks <- samplestorage %>% dplyr::select(ID, contains("Remarks"))
samplestorage <- samplestorage %>% dplyr::select(!(contains("PBMC") | contains("Check"))) %>% 
    rename_with(., ~ str_remove(.x, "Cryovials_")) %>% 
    rename_with(., ~ str_remove(.x, "_Number_of_cryovials")) %>% 
    rename_with(., ~ str_replace(.x, "Lithium_heparin", "heparin")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))

# blood samples
print(samplestorage %>% dplyr::select(ID, contains("Blood") & contains("SOP")), n = 21)
blood <- samplestorage %>% dplyr::select(ID, contains("Blood") & !contains("SOP")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    dplyr::select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
str(blood)

write.csv2(blood, "data/bloodsamples.csv")

# fecal samples; date pivoting separately
feces_date <- samplestorage %>% dplyr::select(ID, contains("Feces_DateTime")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    dplyr::select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
feces <- samplestorage %>% dplyr::select(ID, contains("Feces"), -(contains("FecesCollection")), 
                                                         -(contains("Feces_DateTime"))) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    dplyr::select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value") %>% 
    rename(Samples_stored = Storage)
feces <- right_join(feces, feces_date, by = c("ID", "visit"))
head(feces)

write.csv2(feces, "data/fecessamples.csv")

# urine samples; date pivoting separately
print(samplestorage %>% dplyr::select(ID, contains("Urine") & contains("SOP")), n = 21)
urine_dates <- samplestorage %>% dplyr::select(ID, contains("Urine") & contains("DateTime")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    dplyr::select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value") %>% 
    mutate(across(contains("Date"), ~dmy_hm(.x))) %>% 
    mutate(CollectionTime = EndDateTime - StartDateTime)
urine <- samplestorage %>% dplyr::select(ID, contains("Urine") & !contains("DateTime") &
                                      !contains("SOP")) %>% 
    pivot_longer(., cols = 2:ncol(.), names_to = c("visit", "drop", "variable"),
                 names_sep = "_", values_to = c("value")) %>% 
    dplyr::select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value")
urine <- right_join(urine, urine_dates, by = c("ID", "visit"))
print(urine, n = 63)

saveRDS(urine, "data/urinesamples.RDS")
write.csv2(urine, "data/urinesamples.csv")

#### PBMC ####
str(pbmc)
pbmc %>% dplyr::select(contains("Reason"))
pbmc_sel <- pbmc %>% dplyr::select(ID, contains("Cellcount"), contains("Cryovials")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713"))
pbmc_sel

pbmc_long <- pivot_longer(pbmc_sel, cols = 2:ncol(pbmc_sel), names_sep = "_",
                          names_to = c("visit", "drop", "variable"), values_to = "value") %>% 
    dplyr::select(-drop) %>% 
    pivot_wider(., id_cols = 1:2, names_from = "variable", values_from = "value") %>% 
    dplyr::select(-Cryovials) %>% 
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

#### Capsules ####
str(intervention)
print(intervention %>% dplyr::select(ID, contains("Remarks")))
intervention <- intervention %>% dplyr::select(ID, !(contains("Checks") | contains("Remarks"))) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Total_NonCompliantDays = V3_Capsules_NonCompliantDays + V4_Capsules_NonCompliantDays)
intervention$Total_NonCompliantDays
names(intervention)

saveRDS(intervention, "data/intervention.RDS")
write.csv2(intervention, "data/intervention.csv")

#### Adverse events ####
str(adverse_events)
adverse_events <- adverse_events %>% left_join(., demographics %>% select(ID, Treatment_group), by = "ID")
saveRDS(intervention, "data/adverse_events.RDS")
write.csv2(intervention, "data/adverse_events.csv")
