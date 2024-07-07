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
    statres$p.signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
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
    
    statres <- rbind(statres_line0, statres_line1, statres_line2)
    statres <- tibble::as_tibble(statres)
    statres$p.signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
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
    statres$p.signif <- case_when(
        statres$pval < 0.001 ~paste0("***"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.05 ~paste0("*"),
        statres$pval >= 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
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
           Total_MAP_Mean, Awake_MAP_Mean, Asleep_MAP_Mean,
           Total_PP_Mean, Awake_PP_Mean, Asleep_PP_Mean,
           weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(contains("Mean"), list(mean = mean, sd = sd, n = length), .names = "{.col}_{.fn}"))

## LMMs ABPM
(totsys_lm <- df_total %>% linearmixed(Total_systolic_Mean))
(awsys_lm <- df_total %>% linearmixed(Awake_systolic_Mean))
(aslsys_lm <- df_total %>% linearmixed(Asleep_systolic_Mean))
(totdia_lm <- df_total %>% linearmixed(Total_diastolic_Mean))
(awdia_lm <- df_total %>% linearmixed(Awake_diastolic_Mean))
(asldia_lm <- df_total %>% linearmixed(Asleep_diastolic_Mean))
(totp_lm <- df_total %>% linearmixed(Total_HR_Mean))
(awp_lm <- df_total %>% linearmixed(Awake_HR_Mean))
(aslp_lm <- df_total %>% linearmixed(Asleep_HR_Mean))
(totpp_lm <- df_total %>% linearmixed(Total_PP_Mean))
(awpp_lm <- df_total %>% linearmixed(Awake_PP_Mean))
(aslpp_lm <- df_total %>% linearmixed(Asleep_PP_Mean))
(totmap_lm <- df_total %>% linearmixed(Total_MAP_Mean))
(awmap_lm <- df_total %>% linearmixed(Awake_MAP_Mean))
(aslmap_lm <- df_total %>% linearmixed(Asleep_MAP_Mean))

(plota <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 90, ymax = 170), 
              fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_total, aes(x = weeks, y = Total_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
    geom_point(data = df_total, aes(x = weeks, y = Total_systolic_Mean,
                                    color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
    geom_line(data = df_means, aes(x = weeks, y = Total_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
    geom_point(data = df_means, aes(x = weeks, y = Total_systolic_Mean_mean, 
                                    color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
    geom_errorbar(data = df_means,
                  aes(ymin = Total_systolic_Mean_mean - (Total_systolic_Mean_sd/sqrt(Total_systolic_Mean_n)),
                      ymax = Total_systolic_Mean_mean + (Total_systolic_Mean_sd/sqrt(Total_systolic_Mean_n)),
                      x = weeks,
                      color = Treatment_group), width=0.1, linewidth = 0.8) +
    stat_pvalue_manual(totsys_lm, y.position = 160, label = "p.signif", 
                       tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
    scale_color_jama() + 
    scale_y_continuous(limits = c(90,170), breaks = seq(from = 90, to = 170, by = 10)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Average systolic BP", 
         color = ""))

(plotb <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 90, ymax = 170), 
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Awake_systolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Awake_systolic_Mean,
                                        color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Awake_systolic_Mean_mean, 
                                     color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Awake_systolic_Mean_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Awake_systolic_Mean_mean - (Awake_systolic_Mean_sd/sqrt(Awake_systolic_Mean_n)),
                          ymax = Awake_systolic_Mean_mean + (Awake_systolic_Mean_sd/sqrt(Awake_systolic_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(awsys_lm, y.position = 165, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(90,170), breaks = seq(from = 90, to = 170, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Daytime systolic BP",
             color = ""))

(plotc <- ggplot() +
    geom_rect(aes(xmin = 0, xmax = 4, ymin = 90, ymax = 170), 
              fill = "#CDCDCD", alpha = 0.3) +
    geom_line(data = df_total, aes(x = weeks, y = Asleep_systolic_Mean,
                                   color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
    geom_point(data = df_total, aes(x = weeks, y = Asleep_systolic_Mean,
                                    color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
    geom_line(data = df_means, aes(x = weeks, y = Asleep_systolic_Mean_mean, 
                                   color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
    geom_point(data = df_means, aes(x = weeks, y = Asleep_systolic_Mean_mean, 
                                    color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
    geom_errorbar(data = df_means,
                  aes(ymin = Asleep_systolic_Mean_mean - (Asleep_systolic_Mean_sd/sqrt(Asleep_systolic_Mean_n)),
                      ymax = Asleep_systolic_Mean_mean + (Asleep_systolic_Mean_sd/sqrt(Asleep_systolic_Mean_n)),
                      x = weeks,
                      color = Treatment_group), width=0.1, linewidth = 0.8) +
    stat_pvalue_manual(aslsys_lm, y.position = 165, label = "p.signif", 
                       tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
    scale_color_jama() + 
    scale_y_continuous(limits = c(90,170), breaks = seq(from = 90, to = 170, by = 10)) +
    theme_Publication() +
    labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Nighttime systolic BP", 
         color = ""))


(plotd <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 105), fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Total_diastolic_Mean,
                                     color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Total_diastolic_Mean,
                                    color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Total_diastolic_Mean_mean, 
                                     color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Total_diastolic_Mean_mean,
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Total_diastolic_Mean_mean - (Total_diastolic_Mean_sd/sqrt(Total_diastolic_Mean_n)),
                          ymax = Total_diastolic_Mean_mean + (Total_diastolic_Mean_sd/sqrt(Total_diastolic_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(totdia_lm, y.position = 103, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,105), breaks = seq(from = 45, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Average diastolic BP", 
             color = ""))

(plote <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 105), fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Awake_diastolic_Mean,
                                    color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Awake_diastolic_Mean,
                                    color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Awake_diastolic_Mean_mean, 
                                     color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Awake_diastolic_Mean_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Awake_diastolic_Mean_mean - (Awake_diastolic_Mean_sd/sqrt(Awake_diastolic_Mean_n)),
                          ymax = Awake_diastolic_Mean_mean + (Awake_diastolic_Mean_sd/sqrt(Awake_diastolic_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(awdia_lm, y.position = 103, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,105), breaks = seq(from = 45, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Daytime diastolic BP", 
             color = "Treatment"))
# save_function_bp(plote, group = "abpm", name = "dbp_daytime_abpm", width = 5, height = 4)

(plotf <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 105), fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Asleep_diastolic_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Asleep_diastolic_Mean,
                                        color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Asleep_diastolic_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Asleep_diastolic_Mean_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Asleep_diastolic_Mean_mean - (Asleep_diastolic_Mean_sd/sqrt(Asleep_diastolic_Mean_n)),
                          ymax = Asleep_diastolic_Mean_mean + (Asleep_diastolic_Mean_sd/sqrt(Asleep_diastolic_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(asldia_lm, y.position = 103, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,105), breaks = seq(from = 45, to = 105, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Nighttime diastolic BP", 
             color = "Treatment"))


(plotg <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 85), fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Total_HR_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Total_HR_Mean,
                                        color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Total_HR_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Total_HR_Mean_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Total_HR_Mean_mean - (Total_HR_Mean_sd/sqrt(Total_HR_Mean_n)),
                          ymax = Total_HR_Mean_mean + (Total_HR_Mean_sd/sqrt(Total_HR_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(totp_lm, y.position = 82, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,85), breaks = seq(from = 45, to = 85, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Heart rate (/min)", title = "Average heart rate", 
             color = "Treatment"))

(ploth <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 85), fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Awake_HR_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Awake_HR_Mean,
                                        color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Awake_HR_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Awake_HR_Mean_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Awake_HR_Mean_mean - (Awake_HR_Mean_sd/sqrt(Awake_HR_Mean_n)),
                          ymax = Awake_HR_Mean_mean + (Awake_HR_Mean_sd/sqrt(Awake_HR_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(awp_lm, y.position = 82, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,85), breaks = seq(from = 45, to = 85, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Heart (/min)", title = "Daytime heart rate", 
             color = "Treatment"))

(ploti <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 85), fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_total, aes(x = weeks, y = Asleep_HR_Mean,
                                       color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_total, aes(x = weeks, y = Asleep_HR_Mean,
                                        color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Asleep_HR_Mean_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Asleep_HR_Mean_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Asleep_HR_Mean_mean - (Asleep_HR_Mean_sd/sqrt(Asleep_HR_Mean_n)),
                          ymax = Asleep_HR_Mean_mean + (Asleep_HR_Mean_sd/sqrt(Asleep_HR_Mean_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(aslp_lm, y.position = 82, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,85), breaks = seq(from = 45, to = 85, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Heart rate (/min)", title = "Nighttime heart rate", 
             color = "Treatment"))


(plot_total <- ggarrange(plota, plotb, plotc, plotd, plote, plotf, 
                         plotg, ploth, ploti,
                         nrow = 3, ncol = 3, 
                         labels = c("A", "", "", "B", "", "", "C", "", ""),
                         common.legend = TRUE, legend = "bottom"))
save_function_bp(plot_total, group = "abpm", name = "abpm_lineplots", width = 12, height = 12)

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
save_function_bp(plotg, group = "abpm", name = "total_pulse_lmm")
save_function_bp(ploth, group = "abpm", name = "day_pulse_lmm")
save_function_bp(ploti, group = "abpm", name = "night_pulse_lmm")

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
    summarise(across(c("Systolic", "Diastolic", "Pulse"), list(mean = mean, sd = sd, n = length), 
                     .names = "{.col}_{.fn}"))

officesys_lm <- df_office %>% linearmixed_office(Systolic)

(plot_officesbp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 105, ymax = 185), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_office, aes(x = weeks, y = Systolic,
                                       color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_office, aes(x = weeks, y = Systolic,
                                        color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_officemean, aes(x = weeks, y = Systolic_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_officemean, aes(x = weeks, y = Systolic_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_officemean,
                      aes(ymin = Systolic_mean - (Systolic_sd/sqrt(Systolic_n)),
                          ymax = Systolic_mean + (Systolic_sd/sqrt(Systolic_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(officesys_lm, y.position = 160, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(105,185), breaks = seq(from = 105, to = 185, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Office systolic BP", color = ""))

officedia_lm <- df_office %>% linearmixed_office(Diastolic)

(plot_officedbp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 70, ymax = 115), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_office, aes(x = weeks, y = Diastolic,
                                        color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_office, aes(x = weeks, y = Diastolic,
                                         color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_officemean, aes(x = weeks, y = Diastolic_mean, 
                                            color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_officemean, aes(x = weeks, y = Diastolic_mean, 
                                             color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_officemean,
                      aes(ymin = Diastolic_mean - (Diastolic_sd/sqrt(Diastolic_n)),
                          ymax = Diastolic_mean + (Diastolic_sd/sqrt(Diastolic_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(officedia_lm, y.position = 110, label = "p.signif",
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(70,115), breaks = seq(from = 70, to = 115, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Office diastolic BP", color = ""))

officepulse_lm <- df_office %>% linearmixed_office(Pulse)

(plot_officepulse <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 90), 
                  fill = "#CDCDCD", alpha = 0.4) +
        geom_line(data = df_office, aes(x = weeks, y = Pulse,
                                        color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_office, aes(x = weeks, y = Pulse,
                                         color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_officemean, aes(x = weeks, y = Pulse_mean, 
                                            color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_officemean, aes(x = weeks, y = Pulse_mean, 
                                             color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_officemean,
                      aes(ymin = Pulse_mean - (Pulse_sd/sqrt(Pulse_n)),
                          ymax = Pulse_mean + (Pulse_sd/sqrt(Pulse_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.8) +
        stat_pvalue_manual(officepulse_lm, y.position = 85, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,90), breaks = seq(from = 45, to = 90, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Pulse / min", title = "Office pulse", color = ""))

ggarrange(plot_officesbp, plot_officedbp, plot_officepulse, nrow = 3, ncol = 1,
          labels = c("A", "B", "C"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/officebp/officebp_lineplots_with_lmm.svg", width = 5, height = 12)
ggsave(filename = "results/officebp/officebp_lineplots_with_lmm.pdf", width = 5, height = 12)
ggsave(filename = "results/officebp/officebp_lineplots_with_lmm.png", width = 5, height = 12)
