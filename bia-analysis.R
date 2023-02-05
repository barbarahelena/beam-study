# BIA analysis
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

linearmixed_bia <- function(data, var){
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

save_function_bia <- function(plot, name, a = 5, b = 4){
    ggsave(plot = plot, 
           filename = str_c("results/bia/", name, ".pdf"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/bia/", name, ".svg"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/bia/", name, ".png"), width = a, height = b)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
bia <- readRDS("data/bia_data.RDS")
# lab <- readRDS("data/lab_results.RDS") %>% 
#     select(ID, visit, GFR, LDL, CRP, UrineSodium, contains("DW"), contains("WW"))
urinedata <- readRDS("data/urinesamples.RDS") %>% dplyr::select(ID, visit, Volume)
# dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium, Fibers)
covariates <- right_join(urinedata, df, by = "ID")

#### Dataset for SCFA analysis ####
df_bia <- right_join(bia, covariates, by = c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           ECW = case_when(
               Sex == "Male" ~ (4.8) + (0.225* (V1_Height/100)^2)/(Imp5),
                Sex == "Female" ~ (1.7) + (0.2*(V1_Height/100)^2)/(Imp5)+(0.057*Weight)
           ),
           ICW = TBW - ECW,
           ECWperc = ECW / TBW * 100,
           ICWperc = ICW / TBW * 100) %>% 
    droplevels(.) 

df_means <- df_bia %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("Fatmass", "Fatperc", "Leanmass", "Leanperc", "DryLean",
                       "TBW", "TBWperc", "BMR", "BMI", "BFMI", "FFMI", "ECW", "ICW",
                       "ECWperc", "ICWperc", "Weight"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

#### BIA plots with LMMs ####
weight_lm <- df_bia %>% linearmixed_bia(Weight)
bmi_lm <- df_bia %>% linearmixed_bia(BMI)
ffmi_lm <- df_bia %>% linearmixed_bia(FFMI)
fmi_lm <- df_bia %>% linearmixed_bia(BFMI)
fm_lm <- df_bia %>% linearmixed_bia(Fatmass)
fatp_lm <- df_bia %>% linearmixed_bia(Fatperc)
leanp_lm <- df_bia %>% linearmixed_bia(Leanperc)
drylean_lm <- df_bia %>% linearmixed_bia(DryLean)
tbw_lm <- df_bia %>% linearmixed_bia(TBW)
tbwperc_lm <- df_bia %>% linearmixed_bia(TBWperc)
ecw_lm <- df_bia %>% linearmixed_bia(ECW)
icw_lm <- df_bia %>% linearmixed_bia(ICW)
ecwp_lm <- df_bia %>% linearmixed_bia(ECWperc)
icwp_lm <- df_bia %>% linearmixed_bia(ICWperc)

(plot_weight <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 18, ymax = 30),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Weight_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = Weight,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Weight_mean - (Weight_sd/sqrt(Weight_n)),
                          ymax = Weight_mean + (Weight_sd/sqrt(Weight_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(bmi_lm, y.position = 28, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(50, 105), breaks = seq(from = 50, to = 120, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Weight (kg)", title = "Weight", color = ""))

(plot_bmi <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 18, ymax = 30),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = BMI_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = BMI,
                                      color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = BMI_mean - (BMI_sd/sqrt(BMI_n)),
                          ymax = BMI_mean + (BMI_sd/sqrt(BMI_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(bmi_lm, y.position = 28, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(18,30), breaks = seq(from = 18, to = 30, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "BMI (kg/m2)", title = "BMI", color = ""))

(plot_ffmi <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 12, ymax = 22),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = FFMI_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = FFMI,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = FFMI_mean - (FFMI_sd/sqrt(BMI_n)),
                          ymax = FFMI_mean + (FFMI_sd/sqrt(BMI_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(ffmi_lm, y.position = 20, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(12,22), breaks = seq(from = 12, to = 22, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "FFMI (kg/m2)", title = "FFMI", color = ""))

(plot_fmi <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 12),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = BFMI_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = BFMI,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = BFMI_mean - (BFMI_sd/sqrt(BFMI_n)),
                          ymax = BFMI_mean + (BFMI_sd/sqrt(BFMI_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(fmi_lm, y.position = 10, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,12), breaks = seq(from = 0, to = 12, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "FMI (kg/m2)", title = "FMI", color = ""))

(plot_fatp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 15, ymax = 40),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Fatperc_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = Fatperc,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Fatperc_mean - (Fatperc_sd/sqrt(Fatperc_n)),
                          ymax = Fatperc_mean + (Fatperc_sd/sqrt(Fatperc_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(fatp_lm, y.position = 35, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(15,40), breaks = seq(from = 15, to = 40, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Fat percentage (%)", title = "Fat percentage", color = ""))

(plot_leanp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 55, ymax = 90),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Leanperc_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = Leanperc,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Leanperc_mean - (Leanperc_sd/sqrt(Leanperc_n)),
                          ymax = Leanperc_mean + (Leanperc_sd/sqrt(Leanperc_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(leanp_lm, y.position = 85, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(55,90), breaks = seq(from = 55, to = 90, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Lean weight (%)", title = "Lean weight", color = ""))

(plot_tbwp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 40, ymax = 70),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = TBWperc_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = TBWperc,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = TBWperc_mean - (TBWperc_sd/sqrt(TBWperc_n)),
                          ymax = TBWperc_mean + (TBWperc_sd/sqrt(TBWperc_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(tbwperc_lm, y.position = 65, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama(guide = "none") + 
        scale_y_continuous(limits = c(40,70), breaks = seq(from = 45, to = 70, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Total body water (%)", 
             title = "Total body water", color = ""))


(plot_ebwp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 5, ymax = 20),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = ECWperc_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = ECWperc,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = ECWperc_mean - (ECWperc_sd/sqrt(TBWperc_n)),
                          ymax = ECWperc_mean + (ECWperc_sd/sqrt(TBWperc_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(ecw_lm, y.position = 17.5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(5,20), breaks = seq(from = 5, to = 20, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Extracellular body water (%)", 
             title = "Extracellular body water", color = ""))

(plot_ibwp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 80, ymax = 95),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = ICWperc_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_bia, aes(x = weeks, y = ICWperc,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = ICWperc_mean - (ICWperc_sd/sqrt(ICWperc_n)),
                          ymax = ICWperc_mean + (ICWperc_sd/sqrt(ICWperc_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(icw_lm, y.position = 17.5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(80,95), breaks = seq(from = 80, to = 95, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Intracellular body water (%)", 
             title = "Intracellular body water", color = ""))

(pl_bia <- ggarrange(plot_bmi, plot_ffmi, plot_fmi, plot_fatp, plot_tbwp, plot_ebwp,
                        labels = c("A", "B", "C", "D", "E", "F"),
                        nrow = 2, ncol = 3, common.legend = TRUE,
                     legend = "bottom")
    )
save_function_bia(pl_bia, "bia_plots", a = 9, b = 6)

save_function_bia(plot_bmi, "bmi_bia")
save_function_bia(plot_ffmi, "ffmi_bia")
save_function_bia(plot_fmi, "fmi_bia")
save_function_bia(plot_fatp, "fatp_bia")
save_function_bia(plot_tbwp, "tbwp_bia")
save_function_bia(plot_ebwp, "ecw_bia")
save_function_bia(plot_ibwp, "icw_bia")
save_function_bia(plot_drylean, "drylean_bia")
