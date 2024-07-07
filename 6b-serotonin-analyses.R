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

linearmixed_stn <- function(data, var){
    data1 <- data %>% filter(weeks %in% c(0,4)) %>% 
        mutate(var = {{ var }})
    model1_v4 <- lmer(var ~ Treatment_group*weeks + (1|ID) + (1|Batch), 
                      data = data1)
    res_v4 <- summary(model1_v4)
    print(res_v4)
    pval <- format(round(res_v4$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line1 <- cbind(group1 = 0, group2 = 4, pval)
    
    data2 <- data %>%
        filter(weeks %in% c(4,5)) %>% 
        mutate(var = {{ var }})
    model1_v5 <- lmer(var ~ Treatment_group*weeks + (1|ID) + (1|Batch),
                      data = data2)
    res_v5 <- summary(model1_v5)
    print(res_v5)
    pval <- format(round(res_v5$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line2 <- cbind(group1 = 4, group2 = 5, pval)
    
    statres <- rbind(statres_line1, statres_line2)
    statres <- tibble::as_tibble(statres)
    statres$p.signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        statres$pval > 0.05 ~ paste0("")
    )
    statres$p_signif <- na_if(statres$p.signif, "")
    statres <- statres %>% filter(p.signif != "")
    return(statres)
}


save_function <- function(plot, name){
    ggsave(plot = plot, 
           filename = str_c("results/serotonin/", name, ".pdf"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/serotonin/", name, ".svg"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/serotonin/", name, ".png"), width = 5, height = 4)
}


#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
serotonin <- readRDS("data/serotonin.RDS")

#### Dataset for serotonin analysis ####
df_serotonin <- right_join(df, serotonin, by = c("ID")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           log_stn = log10(serotonin_uM)) %>% 
    droplevels(.) 

df_means <- df_serotonin %>% 
    select(ID, serotonin_uM, log_stn, weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("serotonin_uM","log_stn"),
                 list(mean = ~mean(.x, na.rm = TRUE),
                      sd = ~sd(.x, na.rm = TRUE),
                      n = ~length(.x)
                 ),
                 .names = "{.col}_{.fn}"))

#### Serotonin plot with LMM ####
serotonin_lm <- df_serotonin %>% linearmixed_stn(serotonin_uM) %>% filter(p_signif != "")

(plot_serotonin <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 16),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_serotonin, aes(x = weeks, y = serotonin_uM,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_serotonin, aes(x = weeks, y = serotonin_uM,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = serotonin_uM_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = serotonin_uM_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = serotonin_uM_mean - (serotonin_uM_sd/sqrt(serotonin_uM_n)),
                          ymax = serotonin_uM_mean + (serotonin_uM_sd/sqrt(serotonin_uM_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(serotonin_lm, y.position = 15, label = "p.signif",
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,16), breaks = seq(from = 0, to = 16, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Serotonin (uM)", title = "Serotonin plasma levels", color = ""))

save_function(plot_serotonin, "plasma_serotonin")

logstn_lm <- df_serotonin %>% linearmixed_stn(log_stn)

(plot_log_serotonin <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = -1, ymax = 1.5),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_serotonin, aes(x = weeks, y = log_stn,
                                           color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_serotonin, aes(x = weeks, y = log_stn,
                                            color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = log_stn_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = log_stn_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = log_stn_mean - (log_stn_sd/sqrt(log_stn_n)),
                          ymax = log_stn_mean + (log_stn_sd/sqrt(log_stn_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(logstn_lm, y.position = 0.75, label = "p_signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(-1,1.5), breaks = seq(from = -1, to = 1.5, by = 0.25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Log(Serotonin (uM))", title = "Serotonin plasma levels", color = ""))

save_function(plot_log_serotonin, "plasma_logstn")

