# ABPM adjusted analyses
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(lme4)
library(afex)
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

linearmixed_abpm_cov <- function(data, var){
    data1 <- data %>% filter(weeks %in% c(0,4)) %>% 
        mutate(var = {{ var }})
    ## Model 0
    model0 <- lmer(var ~ (1|ID) + Treatment_group*visit, 
                      data = data1)
    res0 <- summary(model0)
    confint_model0 <- confint(model0)
    print(res0)
    estimate <- as.numeric(format(round(res0$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model0[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model0[6,2], 3), nsmall = 3))
    stderr <- as.numeric(format(round(res0$coefficients[4,2], 3), nsmall = 3))
    pval <- as.numeric(format(round(res0$coefficients[4,5], 3), nsmall = 3))
    statres_line1 <- cbind(model = "model0", estimate, conflow, confhigh, stderr, pval)
    ## Model 1
    model1 <- lmer(var ~ (1|ID) + Age + Sex + BMI + Treatment_group*visit, 
                   data = data1)
    res1 <- summary(model1)
    confint_model1 <- confint(model1)
    print(res1)
    estimate <- as.numeric(format(round(res1$coefficients[7,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model1[9,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model1[9,2], 3), nsmall = 3))
    stderr <- as.numeric(format(round(res1$coefficients[7,2], 3), nsmall = 3))
    pval <- as.numeric(format(round(res1$coefficients[7,5], 3), nsmall = 3))
    statres_line2 <- cbind(model = "model1", estimate, conflow, confhigh, stderr, pval)
    ## Model 2
    model2 <- lmer(var ~ Treatment_group*visit + (1|ID) + Age + Sex + BMI + 
                       Capsules_left + Sodium, 
                   data = data1)
    res2 <- summary(model2)
    print(res2)
    estimate <- as.numeric(format(round(res2$coefficients[9,1], 3), nsmall = 3))
    confint_model2 <- confint(model2)
    conflow <- as.numeric(format(round(confint_model2[11,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model2[11,2], 3), nsmall = 3))
    stderr <- as.numeric(format(round(res2$coefficients[9,2], 3), nsmall = 3))
    pval <- as.numeric(format(round(res2$coefficients[9,5], 3), nsmall = 3))
    statres_line3 <- cbind(model = "model2", estimate, conflow, confhigh, stderr, pval)

    statres <- rbind(statres_line1, statres_line2, statres_line3)
    statres <- as.data.frame(statres) %>% 
        mutate(p_signif = case_when(
                    pval < 0.001 ~paste0("***"),
                    pval < 0.01 ~paste0("**"),
                    pval < 0.05 ~paste0("*"),
                    pval > 0.05 ~paste0("")),
            adjust_for = case_when(
                    model == "model0" ~paste0("unadjusted"),
                    model == "model1" ~paste0("Age, sex, BMI"),
                    model == "model2" ~paste0("adjusted") # age, sex, BMI, sodium, compliance
            ),
            adjust_for = fct_relevel(adjust_for, "unadjusted", after = 0L),
        across(c(2:6), as.numeric))
    return(statres)
}

plot_lmm <- function(results, dfname){
    ylab <- "estimate butyrate*time"
    colors <- pal_lancet()(4)[c(1,3)]
    results <- results %>% filter(model %in% c("model0", "model2"))
    pl1 <- ggplot(results, aes(x=outcome, y=estimate, color=adjust_for, shape = as.factor(p_signif))) +
        geom_hline(yintercept = 0, color = "grey40") +
        geom_point(position=position_dodge(-0.5)) +
        geom_errorbar(aes(ymin=conflow,ymax=confhigh,width=.3), 
                      position=position_dodge(-0.5)) +
        coord_flip()+
        theme_Publication()+
        labs(x = "", y = ylab, 
             title = str_c(dfname)) +
        scale_y_continuous(breaks = seq(from = -10, to = 16, by = 2)) +
        scale_color_manual(values = colors) + 
        scale_shape_manual(values = c(21,19), guide = "none") +
        theme(legend.title = element_blank(),
              legend.position = "bottom")
    return(pl1)
}

save_function_bp <- function(plot, group, name, width = 5, height = 4){
    ggsave(plot = plot, 
           filename = str_c("results/", group, "/", name, ".pdf"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/", group, "/", name, ".svg"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/", group, "/", name, ".png"), width = width, height = height)
}

# Data
df <- readRDS("data/demographics_BEAM.RDS") %>% 
    mutate(V2_time = hms::as_hms(lubridate::dmy_hm(V2_DateTime)),
           V4_time = hms::as_hms(lubridate::dmy_hm(V4_DateTime)),
           V5_time = hms::as_hms(lubridate::dmy_hm(V5_datetime)),
           V4_time_bin = case_when(
               V4_time > hms::as_hms("09:00:00") ~ paste("late"),
               V4_time <= hms::as_hms("09:00:00") ~ paste("early")
           ),
           V4_time_bin = as.factor(V4_time_bin),
           V4_hourdiff = (V4_time - hms::as_hms("07:30:00"))/3600)

df2 <- df %>% pivot_longer(., cols = 25:27,
                           names_to = c("visit", "rest"),
                           names_sep = "_",
                           values_to = "time") %>% 
    dplyr::select(-rest)
bia <- readRDS("data/bia_data.RDS") %>% dplyr::select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS") %>% 
    dplyr::select(ID, visit, GFR, LDL, CRP, UrineSodium)
urinedata <- readRDS("data/urinesamples.RDS") %>% dplyr::select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS") %>% dplyr::select(ID, visit, Sodium, Fibers, Alcohol)
officebp <- readRDS("data/officebp_summary.RDS")
homebp <- readRDS("data/homebp_perweek_BEAM.RDS")
compliance <- readRDS("data/compliance_incl_pharmacy.RDS")
abpm <- readRDS("data/abpm_total.RDS")

#### ABPM ####
covariates <- right_join(bia, 
                         right_join(lab, 
                                    right_join(dietarydata, 
                                               right_join(urinedata,
                                                          right_join(df, compliance))))) %>% 
    mutate(Alc_log = log10(Alcohol+0.1))
df_total <- right_join(abpm, covariates, by= c("ID", "visit")) %>% 
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

#### LMM ####
totsys_lm <- df_total %>% linearmixed_abpm_cov(Total_systolic_Mean) %>% mutate(outcome = "Total systolic BP")
totdia_lm <- df_total %>% linearmixed_abpm_cov(Total_diastolic_Mean) %>% mutate(outcome = "Total diastolic BP")
totpulse_lm <- df_total %>% linearmixed_abpm_cov(Total_HR_Mean) %>% mutate(outcome = "Total heart rate")
daysys_lm <- df_total %>% linearmixed_abpm_cov(Awake_systolic_Mean) %>% mutate(outcome = "Daytime systolic BP")
daydia_lm <- df_total %>% linearmixed_abpm_cov(Awake_diastolic_Mean) %>% mutate(outcome = "Daytime diastolic BP")
daypulse_lm <- df_total %>% linearmixed_abpm_cov(Awake_HR_Mean) %>% mutate(outcome = "Daytime heart rate")
nightsys_lm <- df_total %>% linearmixed_abpm_cov(Asleep_systolic_Mean) %>% mutate(outcome = "Nighttime systolic BP")
nightdia_lm <- df_total %>% linearmixed_abpm_cov(Asleep_diastolic_Mean) %>% mutate(outcome = "Nighttime diastolic BP")
nightpulse_lm <- df_total %>% linearmixed_abpm_cov(Asleep_HR_Mean) %>% mutate(outcome = "Nighttime heart rate")

total_res <- bind_rows(totsys_lm, totdia_lm, totpulse_lm,
                       daysys_lm, daydia_lm, daypulse_lm,
                       nightsys_lm, nightdia_lm, nightpulse_lm) %>% 
    mutate(outcome = fct_inorder(outcome),
           outcome = fct_rev(outcome))
(plot_abpm <- plot_lmm(total_res, dfname = "ABPM linear mixed models"))
save_function_bp(plot_abpm, group = "abpm", width = 5, height = 5, name = "abpm_lmm_2")

#### Office BP ####
df_office <- right_join(officebp, right_join(covariates, df, by = "ID"), by= c("ID", "visit")) %>% 
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

sysoffice_lm <- df_office %>% linearmixed_abpm_cov(Systolic) %>% mutate(outcome = "Office systolic BP")
diaoffice_lm <- df_office %>% linearmixed_abpm_cov(Diastolic) %>% mutate(outcome = "Office diastolic BP")

total_office_res <- bind_rows(sysoffice_lm, diaoffice_lm)
(plot_office <- plot_lmm(total_office_res, dfname = "Office BP linear mixed models"))
save_function_bp(plot_abpm, group = "officebp", width = 6, height = 4, name = "officebp_lmm")

total_bp_res <- bind_rows(totsys_lm, totdia_lm, daysys_lm, daydia_lm, nightsys_lm, nightdia_lm,
                         sysoffice_lm, diaoffice_lm)
(plot_all <- plot_lmm(total_bp_res, dfname = "Linear mixed models"))
save_function_bp(plot_all, group = "abpm", width = 6, height = 5, name = "allbp_lmm")

