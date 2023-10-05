# Serotonin analysis
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

linearmixed_cal <- function(data, var){
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
    model1_v5 <- lmer(var ~ Treatment_group*weeks + (1|ID),
                      data = data2)
    res_v5 <- summary(model1_v5)
    print(res_v5)
    pval <- format(round(res_v5$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line2 <- cbind(group1 = 4, group2 = 5, pval)
    
    statres <- rbind(statres_line1, statres_line2)
    statres <- tibble::as_tibble(statres)
    statres$p_signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        statres$pval > 0.05 ~ paste0("")
    )
    statres$p_signif <- na_if(statres$p_signif, "")
    return(statres)
}


save_function <- function(plot, name){
    ggsave(plot = plot, 
           filename = str_c("results/calprotectin/", name, ".pdf"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/calprotectin/", name, ".svg"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/calprotectin/", name, ".png"), width = 5, height = 4)
}


#### Data ####
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
    
calprotectin <- readRDS("data/calprotectin.RDS")

#### Dataset for serotonin analysis ####
df_calprotectin <- right_join(df, calprotectin, by = c("ID")) %>% 
    mutate(calprotectin_ug_g = na_if(calprotectin_ug_g, "TWM"),
           calprotectin_ug_g = as.numeric(calprotectin_ug_g)) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           log_cal = log10(calprotectin_ug_g)) %>% 
    droplevels(.) 

df_means <- df_calprotectin %>% 
    select(ID, calprotectin_ug_g, log_cal, weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("calprotectin_ug_g","log_cal"),
                 list(mean = ~mean(.x, na.rm = TRUE),
                      sd = ~sd(.x, na.rm = TRUE),
                      n = ~length(.x)
                 ),
                 .names = "{.col}_{.fn}"))

#### Serotonin plot with LMM ####
cal_lm <- c()
cal_lm <- df_calprotectin %>% linearmixed_cal(calprotectin_ug_g)

(plot_calprotectin <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 200),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = calprotectin_ug_g_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_calprotectin, aes(x = weeks, y = calprotectin_ug_g,
                                      color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = calprotectin_ug_g_mean - (calprotectin_ug_g_sd/sqrt(calprotectin_ug_g_n)),
                          ymax = calprotectin_ug_g_mean + (calprotectin_ug_g_sd/sqrt(calprotectin_ug_g_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(cal_lm, y.position = 100, label = "p_signif",
                           remove.bracket = TRUE, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,200), breaks = seq(from = 0, to = 200, by = 20)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Calprotectin (ug/g)", title = "Fecal calprotectin levels", color = ""))

save_function(plot_calprotectin, "fecal_calprotectin")

logcal_lm <- c()
logcal_lm <- df_calprotectin %>% linearmixed_cal(log_cal)

(plot_log_cal <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 2.5),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = log_cal_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_calprotectin, aes(x = weeks, y = log_cal,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = log_cal_mean - (log_cal_sd/sqrt(log_cal_n)),
                          ymax = log_cal_mean + (log_cal_sd/sqrt(log_cal_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(logcal_lm, y.position = 0.75, label = "p_signif", remove.bracket = TRUE, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,2.5), breaks = seq(from = 0, to = 2.5, by = 0.5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Log10(calprotectin (ug/g))", title = "Fecal calprotectin levels", color = ""))

save_function(plot_log_cal, "plasma_logcal")

