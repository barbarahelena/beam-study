# BP analyses
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(lme4)
library(afex)

# Data
df <- readRDS("data/demographics_BEAM.RDS")
bia <- readRDS("data/bia_data.RDS") %>% select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS") %>% 
    select(ID, visit, GFR, LDL, CRP, UrineSodium)
urinedata <- readRDS("data/urinesamples.RDS") %>% select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium)
officebp <- readRDS("data/officebp_summary.RDS")
homebp <- readRDS("data/homebp_perweek_BEAM.RDS")
abpm <- readRDS("data/abpm_total.RDS")

# ABPM
covariates <- right_join(bia, 
                         right_join(lab, 
                                    right_join(dietarydata, urinedata, by = c("ID", "visit")), 
                                    by = c("ID", "visit")),
                            by = c("ID", "visit"))
df_total <- right_join(abpm, right_join(covariates, df, by = "ID"), by= c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks)) %>% 
    droplevels(.) 

df_means <- df_total %>% 
    select(ID, Total_systolic_Mean, Total_diastolic_Mean,
           Awake_systolic_Mean, Awake_diastolic_Mean,
           Asleep_systolic_Mean, Asleep_diastolic_Mean, weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(contains("Mean"), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

model1 <- lmer(Total_systolic_Mean ~ Treatment_group*weeks + (1|ID), 
               data = df_total %>% filter(weeks %in% c(0, 4)))
res <- summary(model1)
pval_v4 <- summary(model1)[3,5]

(plota <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 160), 
              fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_means, aes(x = weeks, y = Total_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1) +
    geom_line(data = df_total, aes(x = weeks, y = Total_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2) +
    geom_errorbar(data = df_means,
                  aes(ymin = Total_systolic_Mean_mean - (Total_systolic_Mean_sd/sqrt(21)),
                      ymax = Total_systolic_Mean_mean + (Total_systolic_Mean_sd/sqrt(21)),
                      x = weeks,
                      color = Treatment_group), width=0.1) +
    scale_color_jama(guide = "none") + 
    scale_y_continuous(limits = c(100,160)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic blood pressure (mmHg)", title = "Average systolic BP"))


statres <- c()

group1 <- 0
group2 <- 4
model1_v4 <- lmer(Awake_systolic_Mean ~ Treatment_group*weeks + (1|ID), 
               data = df_total %>% filter(weeks %in% c(0, 4)))
res_v4 <- summary(model1_v4)
pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
pval <- as.numeric(pval)
statres_line1 <- cbind(group1, group2, pval)

group1 <- 4
group2 <- 5
model1_v5 <- lmer(Awake_systolic_Mean ~ Treatment_group*weeks + (1|ID), 
                  data = df_total %>% filter(weeks %in% c(4, 5)))
res_v5 <- summary(model1_v5)
pval <- format(round(res_v5$coefficients[4,5], 3), nsmall = 3)
pval <- as.numeric(pval)
statres_line2 <- cbind(group1, group2, pval)

statres <- rbind(statres_line1, statres_line2)
statres <- tibble::as_tibble(statres)
statres$p_signif <- case_when(
    statres$pval < 0.05 ~paste0("*"),
    statres$pval < 0.01 ~paste0("**"),
    statres$pval < 0.001 ~paste0("***"),
    statres$pval > 0.05 ~paste0("")
)

(plotb <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 160), 
                  fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_means, aes(x = weeks, y = Awake_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1) +
    geom_line(data = df_total, aes(x = weeks, y = Awake_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2) +
    geom_errorbar(data = df_means,
                  aes(ymin = Awake_systolic_Mean_mean - (Awake_systolic_Mean_sd/sqrt(21)),
                      ymax = Awake_systolic_Mean_mean + (Awake_systolic_Mean_sd/sqrt(21)),
                      x = weeks,
                      color = Treatment_group), width=0.1) +
    stat_pvalue_manual(statres, y.position = 150, label = "p_signif", remove.bracket = TRUE) +
    scale_color_jama(guide = "none") + 
    scale_y_continuous(limits = c(100,160)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic blood pressure (mmHg)", title = "Daytime systolic BP"))

(plotc <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 160), 
              fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_means, aes(x = weeks, y = Asleep_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1) +
    geom_line(data = df_total, aes(x = weeks, y = Asleep_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2) +
    geom_errorbar(data = df_means,
                  aes(ymin = Asleep_systolic_Mean_mean - (Asleep_systolic_Mean_sd/sqrt(21)),
                      ymax = Asleep_systolic_Mean_mean + (Asleep_systolic_Mean_sd/sqrt(21)),
                      x = weeks,
                      color = Treatment_group), width=0.1) +
    scale_color_jama(guide = "none") + 
    scale_y_continuous(limits = c(100,160)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic blood pressure (mmHg)", title = "Nighttime systolic BP"))


(plotd <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 160), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Total_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Total_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Total_diastolic_Mean_mean - (Total_diastolic_Mean_sd/sqrt(21)),
                          ymax = Total_diastolic_Mean_mean + (Total_diastolic_Mean_sd/sqrt(21)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(45,100)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic blood pressure (mmHg)", title = "Average diastolic BP"))

(plote <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 160), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Awake_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Awake_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Awake_diastolic_Mean_mean - (Awake_diastolic_Mean_sd/sqrt(21)),
                          ymax = Awake_diastolic_Mean_mean + (Awake_diastolic_Mean_sd/sqrt(21)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(45, 100)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic blood pressure (mmHg)", title = "Daytime diastolic BP"))

(plotf <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 160), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Asleep_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Asleep_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Asleep_diastolic_Mean_mean - (Asleep_diastolic_Mean_sd/sqrt(21)),
                          ymax = Asleep_diastolic_Mean_mean + (Asleep_diastolic_Mean_sd/sqrt(21)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(45,100)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic blood pressure (mmHg)", title = "Nighttime diastolic BP"))

