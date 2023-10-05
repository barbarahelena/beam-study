# Dietary analysis
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

linearmixed_diet <- function(data, var){
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
        statres$pval > 0.05 ~paste0("")
    )
    return(statres)
}

save_function_diet <- function(plot, name, a = 5, b = 4){
    ggsave(plot = plot, 
           filename = str_c("results/diet/", name, ".pdf"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/diet/", name, ".svg"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/diet/", name, ".png"), width = a, height = b)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
# bia <- readRDS("data/bia_data.RDS")
# lab <- readRDS("data/lab_results.RDS") %>% 
#     select(ID, visit, GFR, LDL, CRP, UrineSodium, contains("DW"), contains("WW"))
# urinedata <- readRDS("data/urinesamples.RDS") %>% select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS")

#### Dataset for SCFA analysis ####
df_diet <- right_join(df, dietarydata, by = c("ID")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           Alc_log = log10(Alcohol+1)) %>% 
    droplevels(.) 

df_means <- df_diet %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("Energy", "Fat", "SaturatedFat", "Proteins", "Carbohydrates",
                       "Fibers", "Salt", "Sodium", "Alcohol", "Alc_log"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

#### Dietary data plots with LMMs ####
energy_lm <- df_diet %>% linearmixed_diet(Energy)
fat_lm <- df_diet %>% linearmixed_diet(Fat)
satfat_lm <- df_diet %>% linearmixed_diet(SaturatedFat)
protein_lm <- df_diet %>% linearmixed_diet(Proteins)
carbs_lm <- df_diet %>% linearmixed_diet(Carbohydrates)
fibers_lm <- df_diet %>% linearmixed_diet(Fibers)
salt_lm <- df_diet %>% linearmixed_diet(Salt)
sodium_lm <- df_diet %>% linearmixed_diet(Sodium)
alc_lm <- df_diet %>% linearmixed_diet(Alcohol)
alc_log_lm <- df_diet %>% linearmixed_diet(Alc_log)

(plot_energy <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 500, ymax = 3000),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Energy_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = Energy,
        #                              color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Energy_mean - (Energy_sd/sqrt(Energy_n)),
                          ymax = Energy_mean + (Energy_sd/sqrt(Energy_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(energy_lm, y.position = 2800, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(500,3000), breaks = seq(from = 500, to = 3000, by = 500)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Energy (kcal)", title = "Energy intake",
             color = ""))

(plot_fat <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 125),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Fat_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        # geom_line(data = df_diet, aes(x = weeks, y = Fat,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Fat_mean - (Fat_sd/sqrt(Fat_n)),
                          ymax = Fat_mean + (Fat_sd/sqrt(Fat_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(fat_lm, y.position = 100, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,125), breaks = seq(from = 0, to = 125, by = 25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Fat (g)", title = "Fat",
             color = ""))

(plot_satfat <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 75),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = SaturatedFat_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = SaturatedFat,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = SaturatedFat_mean - (SaturatedFat_sd/sqrt(SaturatedFat_n)),
                          ymax = SaturatedFat_mean + (SaturatedFat_sd/sqrt(SaturatedFat_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(satfat_lm, y.position = 60, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,75), breaks = seq(from = 0, to = 75, by = 25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Saturated fat (g)", title = "Saturated fat",
             color = ""))

(plot_protein <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 150),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Proteins_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = Proteins,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Proteins_mean - (Proteins_sd/sqrt(Proteins_n)),
                          ymax = Proteins_mean + (Proteins_sd/sqrt(Proteins_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(protein_lm, y.position = 140, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,150), breaks = seq(from = 0, to = 150, by = 25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Proteins (g)", title = "Proteins",
             color = ""))

(plot_fibers <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 60),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Fibers_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = Fibers,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Fibers_mean - (Fibers_sd/sqrt(Fibers_n)),
                          ymax = Fibers_mean + (Fibers_sd/sqrt(Fibers_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(fibers_lm, y.position = 50, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,60), breaks = seq(from = 0, to = 60, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Fibers (g)", title = "Fibers",
             color = ""))

(plot_sodium <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 6500),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Sodium_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = Sodium,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Sodium_mean - (Sodium_sd/sqrt(Sodium_n)),
                          ymax = Sodium_mean + (Sodium_sd/sqrt(Sodium_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(sodium_lm, y.position = 4000, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,6750), breaks = seq(from = 0, to = 6000, by = 1000)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Sodium (mg)", title = "Sodium",
             color = ""))

(plot_alcohol <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 65),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Alcohol_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = Alcohol,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Alcohol_mean - (Alcohol_sd/sqrt(Alcohol_n)),
                          ymax = Alcohol_mean + (Alcohol_sd/sqrt(Alcohol_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(sodium_lm, y.position = 20, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        #scale_y_continuous(limits = c(0,6750), breaks = seq(from = 0, to = 6000, by = 1000)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Alcohol (g)", title = "Alcohol",
             color = ""))

(plot_alcohol_log <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 2),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Alc_log_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_diet, aes(x = weeks, y = Alc_log,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Alc_log_mean - (Alc_log_sd/sqrt(Alc_log_n)),
                          ymax = Alc_log_mean + (Alc_log_sd/sqrt(Alc_log_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(sodium_lm, y.position = 2, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        #scale_y_continuous(limits = c(0,6750), breaks = seq(from = 0, to = 6000, by = 1000)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Alcohol (g)", title = "Alcohol",
             color = ""))


(pl_diet <- ggarrange(plot_energy, plot_fat, plot_satfat, plot_protein, 
                      plot_fibers, plot_sodium, plot_alcohol,
                     labels = c("A", "B", "C", "D", "E", "F", "G"),
                     nrow = 3, ncol = 3,
                     common.legend = TRUE,
                     legend = "bottom"))
save_function_diet(pl_diet, "diet_plots", a = 9, b = 9)

save_function_diet(plot_energy, "energy_diet")
save_function_diet(plot_fat, "fat_diet")
save_function_diet(plot_satfat, "satfat_diet")
save_function_diet(plot_protein, "protein_diet")
save_function_diet(plot_fibers, "fibers_diet")
save_function_diet(plot_sodium, "sodium_diet")
save_function_diet(plot_alcohol, "alcohol_diet")
save_function_diet(plot_alcohol_log, "alc_log_diet")
