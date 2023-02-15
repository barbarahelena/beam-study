# BP plots: ABPM, home BP and office BP
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(lme4)
library(afex)
library(ggpubr)
library(ggsci)
library(patchwork)

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
                legend.position = "bottom",
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

linearmixed <- function(data, var){
    data1 <- data %>% filter(weeks %in% c(0,4)) %>% 
        mutate(var = {{ var }})
    model1_v4 <- lmer(var ~ Treatment_group*weeks + (1|ID), 
                      data = data1)
    res_v4 <- summary(model1_v4)
    print(res_v4)
    pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line1 <- cbind(group1 = 0, group2 = 4, pval)
    
    data2 <- data %>%
        filter(weeks %in% c(4,5)) %>% 
        mutate(var = {{ var }})
    model1_v5 <- lmer(var ~ Treatment_group*visit + (1|ID),
                      data = data2)
    res_v5 <- summary(model1_v5)
    print(res_v5)
    pval <- format(round(res_v5$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line2 <- cbind(group1 = 4, group2 = 5, pval)
    
    statres <- rbind(statres_line1, statres_line2)
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    return(statres)
}

linearmixed_office <- function(data, var){
    data0 <- data %>% filter(weeks %in% c(0,2)) %>% 
        mutate(var = {{ var }})
    model1_v2 <- lmer(var ~ Treatment_group*weeks + (1|ID), 
                      data = data0)
    res_v2 <- summary(model1_v2)
    print(res_v2)
    pval <- format(round(res_v2$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line0 <- cbind(group1 = 0, group2 = 2, pval)
    
    data1 <- data %>% filter(weeks %in% c(0,4)) %>% 
        mutate(var = {{ var }})
    model1_v4 <- lmer(var ~ Treatment_group*weeks + (1|ID), 
                      data = data1)
    res_v4 <- summary(model1_v4)
    print(res_v4)
    pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line1 <- cbind(group1 = 2, group2 = 4, pval)
    
    data2 <- data %>%
        filter(weeks %in% c(4,5)) %>% 
        mutate(var = {{ var }})
    model1_v5 <- lmer(var ~ Treatment_group*visit + (1|ID),
                      data = data2)
    res_v5 <- summary(model1_v5)
    print(res_v5)
    pval <- format(round(res_v5$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line2 <- cbind(group1 = 4, group2 = 5, pval)
    
    statres <- rbind(statres_line1, statres_line2)
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    return(statres)
}
linearmixed_homebp <- function(data, var){
    statres <- c()
    for(a in c(2:5)){
        data1 <- data %>% filter(week %in% c(a-1,a)) %>% 
            mutate(var = {{ var }})
        model <- lmer(var ~ Treatment_group*week + (1|ID), 
                      data = data1)
        res <- summary(model)
        pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
        pval <- as.numeric(pval)
        statres_line <- cbind(group1 = a-1, group2 = a, pval)
        statres <- rbind(statres, statres_line)
        print(statres)
        print(res)
    }
    
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    return(statres)
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
df <- readRDS("data/demographics_BEAM.RDS")
bia <- readRDS("data/bia_data.RDS") %>% dplyr::select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS") %>% 
    dplyr::select(ID, visit, GFR, LDL, CRP, UrineSodium)
urinedata <- readRDS("data/urinesamples.RDS") %>% dplyr::select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS") %>% dplyr::select(ID, visit, Sodium)
officebp <- readRDS("data/officebp_summary.RDS")
homebp <- readRDS("data/homebp_perweek_BEAM.RDS")
abpm <- readRDS("data/abpm_total.RDS")

#### ABPM ####
covariates <- right_join(bia, 
                         right_join(lab, 
                                    right_join(dietarydata, urinedata, by = c("ID", "visit")), 
                                    by = c("ID", "visit")),
                            by = c("ID", "visit"))
df_total <- right_join(abpm, right_join(covariates, df, by = "ID"), by= c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks)) %>% 
    droplevels(.) 

df_means <- df_total %>% 
    dplyr::select(ID, Total_systolic_Mean, Total_diastolic_Mean, Total_HR_Mean,
           Awake_systolic_Mean, Awake_diastolic_Mean, Awake_HR_Mean,
           Asleep_systolic_Mean, Asleep_diastolic_Mean, Asleep_HR_Mean,
           weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(contains("Mean"), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

totsys_lm <- c()
totsys_lm <- df_total %>% linearmixed(Total_systolic_Mean)

(plota <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 90, ymax = 165), 
              fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_means, aes(x = weeks, y = Total_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1) +
    geom_line(data = df_total, aes(x = weeks, y = Total_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2) +
    geom_errorbar(data = df_means,
                  aes(ymin = Total_systolic_Mean_mean - (Total_systolic_Mean_sd/sqrt(11)),
                      ymax = Total_systolic_Mean_mean + (Total_systolic_Mean_sd/sqrt(11)),
                      x = weeks,
                      color = Treatment_group), width=0.1) +
    stat_pvalue_manual(totsys_lm, y.position = 150, label = "p_signif", 
                       remove.bracket = TRUE, bracket.size = 0) +
    scale_color_jama() + 
    scale_y_continuous(limits = c(90,165), breaks = seq(from = 90, to = 160, by = 10)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Average systolic BP", 
         color = "Treatment"))

awsys_lm <- c()
awsys_lm <- df_total %>% linearmixed(Awake_systolic_Mean)

(plotb <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 105, ymax = 165), 
                  fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_means, aes(x = weeks, y = Awake_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1) +
    geom_line(data = df_total, aes(x = weeks, y = Awake_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2) +
    geom_errorbar(data = df_means,
                  aes(ymin = Awake_systolic_Mean_mean - (Awake_systolic_Mean_sd/sqrt(11)),
                      ymax = Awake_systolic_Mean_mean + (Awake_systolic_Mean_sd/sqrt(11)),
                      x = weeks,
                      color = Treatment_group), width=0.1) +
    stat_pvalue_manual(awsys_lm, y.position = 150, label = "p_signif", 
                       remove.bracket = TRUE, bracket.size = 0) +
    scale_color_jama() + 
    scale_y_continuous(limits = c(105,165), breaks = seq(from = 105, to = 170, by = 10)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Daytime systolic BP",
         color = "Treatment"))

aslsys_lm <- c()
aslsys_lm <- df_total %>% linearmixed(Asleep_systolic_Mean)

(plotc <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 90, ymax = 165), 
                  fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_means, aes(x = weeks, y = Asleep_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1) +
    geom_line(data = df_total, aes(x = weeks, y = Asleep_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2) +
    geom_errorbar(data = df_means,
                  aes(ymin = Asleep_systolic_Mean_mean - (Asleep_systolic_Mean_sd/sqrt(11)),
                      ymax = Asleep_systolic_Mean_mean + (Asleep_systolic_Mean_sd/sqrt(11)),
                      x = weeks,
                      color = Treatment_group), width=0.1) +
    stat_pvalue_manual(aslsys_lm, y.position = 150, label = "p_signif", 
                       remove.bracket = TRUE, bracket.size = 0) +
    scale_color_jama() + 
    scale_y_continuous(limits = c(90,165), breaks = seq(from = 90, to = 160, by = 10)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Nighttime systolic BP", 
         color = "Treatment"))

totdia_lm <- c()
(totdia_lm <- df_total %>% linearmixed(Total_diastolic_Mean))

(plotd <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 105), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Total_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Total_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Total_diastolic_Mean_mean - (Total_diastolic_Mean_sd/sqrt(11)),
                          ymax = Total_diastolic_Mean_mean + (Total_diastolic_Mean_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(totdia_lm, y.position = 100, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,105), breaks = seq(from = 45, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Average diastolic BP", 
             color = "Treatment"))

awdia_lm <- c()
(awdia_lm <- df_total %>% linearmixed(Awake_diastolic_Mean))

(plote <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 105), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Awake_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Awake_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Awake_diastolic_Mean_mean - (Awake_diastolic_Mean_sd/sqrt(11)),
                          ymax = Awake_diastolic_Mean_mean + (Awake_diastolic_Mean_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(awdia_lm, y.position = 100, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,105), breaks = seq(from = 45, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Daytime diastolic BP", 
             color = "Treatment"))

asldia_lm <- c()
(asldia_lm <- df_total %>% linearmixed(Asleep_diastolic_Mean))

(plotf <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 105), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_means, aes(x = weeks, y = Asleep_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Asleep_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Asleep_diastolic_Mean_mean - (Asleep_diastolic_Mean_sd/sqrt(11)),
                          ymax = Asleep_diastolic_Mean_mean + (Asleep_diastolic_Mean_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(asldia_lm, y.position = 100, label = "p_signif", 
                           hide.ns = TRUE, remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,105), breaks = seq(from = 45, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Nighttime diastolic BP", 
             color = "Treatment"))


totp_lm <- c()
(totp_lm <- df_total %>% linearmixed(Total_HR_Mean))

(plotg <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 85), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Total_HR_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Total_HR_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Total_HR_Mean_mean - (Total_HR_Mean_sd/sqrt(11)),
                          ymax = Total_HR_Mean_mean + (Total_HR_Mean_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(totp_lm, y.position = 82, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,85), breaks = seq(from = 45, to = 85, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Heart rate (/min)", title = "Total heart rate", 
             color = "Treatment"))

awp_lm <- c()
(awp_lm <- df_total %>% linearmixed(Awake_HR_Mean))

(ploth <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 85), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Awake_HR_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Awake_HR_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Awake_HR_Mean_mean - (Awake_HR_Mean_sd/sqrt(11)),
                          ymax = Awake_HR_Mean_mean + (Awake_HR_Mean_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(awp_lm, y.position = 82, label = "p_signif", 
                           remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,85), breaks = seq(from = 45, to = 85, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Heart (/min)", title = "Daytime heart rate", 
             color = "Treatment"))

aslp_lm <- c()
(aslp_lm <- df_total %>% linearmixed(Asleep_HR_Mean))

(ploti <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 85), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_means, aes(x = weeks, y = Asleep_HR_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_total, aes(x = weeks, y = Asleep_HR_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Asleep_HR_Mean_mean - (Asleep_HR_Mean_sd/sqrt(11)),
                          ymax = Asleep_HR_Mean_mean + (Asleep_HR_Mean_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(aslp_lm, y.position = 82, label = "p_signif", 
                           hide.ns = TRUE, remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,85), breaks = seq(from = 45, to = 85, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Heart rate (/min)", title = "Nighttime heart rate", 
             color = "Treatment"))


(plot_total <- ggarrange(plota, plotb, plotc, plotd, plote, plotf, 
                         plotg, ploth, ploti,
                         nrow = 3, ncol = 3, 
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                         common.legend = TRUE, legend = "bottom"))
save_function_bp(plot_total, group = "abpm", name = "abpm_lineplots_lmm", width = 11, height = 12)

plot_pulse <- ggarrange(plotg, ploth, ploti, labels = c("A", "B", "C"), 
                        common.legend = TRUE,
                        nrow = 1,
                        legend = "bottom")
save_function_bp(plot_pulse, group = "abpm", name = "abpm_pulse", width = 9, height = 4)

# plot_emphasis <- ((plotb / plote ) | 
#     ((plota | plotc) / (plotd | plotf)) & theme(legend.position = "bottom")) + 
#     plot_layout(guides = "collect")
# save_function_bp(plot_emphasis, group = "abpm", name = "abpm_lineplots_emphasis", 
#                  width = 10, height = 8)

save_function_bp(plota, group = "abpm", name = "total_sbp_lmm")
save_function_bp(plotb, group = "abpm", name = "day_sbp_lmm")
save_function_bp(plotc, group = "abpm", name = "night_sbp_lmm")
save_function_bp(plotd, group = "abpm", name = "total_dbp_lmm")
save_function_bp(plote, group = "abpm", name = "day_dbp_lmm")
save_function_bp(plotf, group = "abpm", name = "night_dbp_lmm")

#### Home BP ####
str(homebp)
df_homebp <- full_join(homebp, df, by= c("ID")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group)) %>% 
    mutate(week = as.numeric(week),
        week = case_when(
            week < 0 ~ paste0(week + 1),
            week > 0 ~ paste0(week),
    ), week = as.numeric(week)) %>% 
    filter(week %in% c(1:5)) %>%
    droplevels(.) 

df_homebp_filter <- df_homebp %>% select(ID, contains("sbp")) %>% 
    group_by(ID) %>% 
    summarise(., nas = sum(is.na(mean_sbp))) %>%
    filter(., nas > 0)

df_homebp <- df_homebp %>% filter(!ID %in% df_homebp_filter$ID)

df_mean_homebp <- df_homebp %>% 
    select(ID, mean_sbp, mean_dbp, mean_pulse, week, Treatment_group) %>% 
    group_by(Treatment_group, week) %>% 
    summarise(across(c("mean_sbp", "mean_dbp", "mean_pulse"), list(mean = mean, sd = sd), na.rm = TRUE,
                     .names = "{.col}_{.fn}"))

lm_home_sbp <- linearmixed_homebp(df_homebp, mean_sbp)
lm_home_dbp <- linearmixed_homebp(df_homebp, mean_dbp)
lm_home_pulse <- linearmixed_homebp(df_homebp, mean_pulse)

(plot_homebp_1 <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 115, ymax = 180), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_mean_homebp, aes(x = week, y = mean_sbp_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_homebp, aes(x = week, y = mean_sbp,
                                       color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_mean_homebp,
                      aes(ymin = mean_sbp_mean - (mean_sbp_sd/sqrt(11)),
                          ymax = mean_sbp_mean + (mean_sbp_sd/sqrt(11)),
                          x = week,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(lm_home_sbp, y.position = 160, label = "p_signif", 
                           hide.ns = TRUE, remove.bracket = TRUE) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(115,180), breaks = seq(from = 115, to = 180, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Home systolic BP"))

(plot_homebp_2 <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 70, ymax = 107), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_mean_homebp, aes(x = week, y = mean_dbp_mean, 
                            color = Treatment_group, group = Treatment_group), 
                            alpha = 1) +
        geom_line(data = df_homebp, aes(x = week, y = mean_dbp,
                                        color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_mean_homebp,
                      aes(ymin = mean_dbp_mean - (mean_dbp_sd/sqrt(11)),
                          ymax = mean_dbp_mean + (mean_dbp_sd/sqrt(11)),
                          x = week,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(lm_home_dbp, y.position = 100, label = "p_signif", 
                           remove.bracket = TRUE, hide.ns = TRUE) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(70, 107), breaks = seq(from = 70, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Home diastolic BP"))

(plot_homebp_3 <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 90), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_mean_homebp, aes(x = week, y = mean_pulse_mean, 
                                             color = Treatment_group, group = Treatment_group), 
                  alpha = 1) +
        geom_line(data = df_homebp, aes(x = week, y = mean_pulse,
                                        color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_mean_homebp,
                      aes(ymin = mean_pulse_mean - (mean_pulse_sd/sqrt(11)),
                          ymax = mean_pulse_mean + (mean_pulse_sd/sqrt(11)),
                          x = week,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(lm_home_dbp, y.position = 85, label = "p_signif", 
                           remove.bracket = TRUE, hide.ns = TRUE) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(45, 90), breaks = seq(from = 45, to = 90, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Pulse / min", title = "Home pulse"))

ggarrange(plot_homebp_1, plot_homebp_2, plot_homebp_3, nrow = 1, ncol = 3)
ggsave(filename = "results/homebp/homebp_lineplots_with_lmm.svg", width = 12, height = 4)
ggsave(filename = "results/homebp/homebp_lineplots_with_lmm.pdf", width = 12, height = 4)

save_function_bp(plot_homebp_1, "homebp", "home_sbp")
save_function_bp(plot_homebp_1, "homebp", "home_dbp")
save_function_bp(plot_homebp_1, "homebp", "home_pulse")

#### Office BP ####
df_office <- officebp %>% left_join(., df, by = "ID") %>% 
    #right_join(right_join(covariates, df, by = "ID"), officebp, by= c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V3" ~ paste0(2),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks)) %>% 
    droplevels(.) 

df_officemean <- df_office %>% 
    select(ID, Systolic, Diastolic, Pulse, weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("Systolic", "Diastolic", "Pulse"), list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))

officesys_lm <- df_office %>% linearmixed_office(Systolic)

(plot_officesbp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 105, ymax = 180), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_officemean, aes(x = weeks, y = Systolic_mean, 
                                             color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_office, aes(x = weeks, y = Systolic,
                                        color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_officemean,
                      aes(ymin = Systolic_mean - (Systolic_sd/sqrt(11)),
                          ymax = Systolic_mean + (Systolic_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(officesys_lm, y.position = 160, label = "p_signif", 
                           hide.ns = TRUE, remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(105,180), breaks = seq(from = 105, to = 180, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Office systolic BP", color = ""))

officedia_lm <- df_office %>% linearmixed_office(Diastolic)

(plot_officedbp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 70, ymax = 115), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_officemean, aes(x = weeks, y = Diastolic_mean, 
                                            color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_office, aes(x = weeks, y = Diastolic,
                                        color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_officemean,
                      aes(ymin = Diastolic_mean - (Diastolic_sd/sqrt(11)),
                          ymax = Diastolic_mean + (Diastolic_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(officedia_lm, y.position = 105, label = "p_signif",
                           hide.ns = TRUE, remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(70,115), breaks = seq(from = 70, to = 115, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Office diastolic BP", color = ""))

officepulse_lm <- df_office %>% linearmixed_office(Pulse)

(plot_officepulse <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 90), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_officemean, aes(x = weeks, y = Pulse_mean, 
                                            color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_office, aes(x = weeks, y = Pulse,
                                        color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_officemean,
                      aes(ymin = Pulse_mean - (Pulse_sd/sqrt(11)),
                          ymax = Pulse_mean + (Pulse_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(officepulse_lm, y.position = 85, label = "p_signif", 
                           hide.ns = TRUE, remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,90), breaks = seq(from = 45, to = 90, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Pulse / min", title = "Office pulse", color = ""))

ggarrange(plot_officesbp, plot_officedbp, plot_officepulse, nrow = 1, ncol = 3,
          labels = c("A", "B", "C"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/officebp/officebp_lineplots_with_lmm.svg", width = 12, height = 4)
ggsave(filename = "results/officebp/officebp_lineplots_with_lmm.pdf", width = 12, height = 4)
ggsave(filename = "results/officebp/officebp_lineplots_with_lmm.png", width = 12, height = 4)

officebp_boxplots <- df_office %>% 
    mutate(before_after = case_when(
        visit == "V2" ~ paste0("Before"),
        visit == "V4" ~ paste0("Treatment"),
        visit == "V5" ~ paste0("After")
    ),
    before_after = fct_relevel(before_after, "Before", after = 0L),
    before_after = fct_relevel(before_after, "After", after = 2L))

(boxplot_officesbp <- ggplot(officebp_boxplots, aes(x = before_after, y = Systolic)) +
        geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        stat_compare_means(method = "t.test", label = "p.value",
                           comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
        facet_wrap(~Treatment_group) + 
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(x = "", title = "Office systolic BP", y = "Systolic BP (mmHg)") +
        theme_Publication() )
save_function_bp(boxplot_officesbp, "officebp", "boxplot_officesbp", height = 5)

(boxplot_officedbp <- ggplot(officebp_boxplots, aes(x = before_after, y = Diastolic)) +
        geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        stat_compare_means(method = "t.test", label = "p.value",
                           comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
        facet_wrap(~Treatment_group) + 
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(x = "", title = "Office diastolic BP", y = "Diastolic BP (mmHg)") +
        theme_Publication() )
save_function_bp(boxplot_officesbp, "officebp", "boxplot_officedbp", height = 5)

