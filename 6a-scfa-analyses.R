# SCFA analysis
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

linearmixed_scfa <- function(data, var){
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

linearmixed_scfa_cov <- function(data, var){
    data1 <- data %>% filter(weeks %in% c(0,4)) %>% 
        mutate(var = {{ var }})
    model1_v4 <- lmer(var ~ Treatment_group*weeks + (1|ID) + V4_hourdiff + Capsules_left,
                      data = data1)
    res_v4 <- summary(model1_v4)
    print(res_v4)
    pval <- format(round(res_v4$coefficients[6,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line1 <- cbind(group1 = 0, group2 = 4, pval)
    
    data2 <- data %>%
        filter(weeks %in% c(4,5)) %>% 
        mutate(var = {{ var }})
    model1_v5 <- lmer(var ~ Treatment_group*weeks + (1|ID) + V4_hourdiff + Capsules_left,
                      data = data2)
    res_v5 <- summary(model1_v5)
    print(res_v5)
    pval <- format(round(res_v5$coefficients[6,5], 3), nsmall = 3)
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

save_function <- function(plot, name){
    ggsave(plot = plot, 
           filename = str_c("results/fecalscfa/", name, ".pdf"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/fecalscfa/", name, ".svg"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/fecalscfa/", name, ".png"), width = 5, height = 4)
}

save_function_plasma <- function(plot, name){
    ggsave(plot = plot, 
           filename = str_c("results/plasmascfa/", name, ".pdf"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/plasmascfa/", name, ".svg"), width = 5, height = 4)
    ggsave(plot = plot, 
           filename = str_c("results/plasmascfa/", name, ".png"), width = 5, height = 4)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")

df2 <- df %>% pivot_longer(., cols = 25:27,
                 names_to = c("visit", "rest"),
                 names_sep = "_",
                 values_to = "time") %>% 
    dplyr::select(-rest)
    
bia <- readRDS("data/bia_data.RDS") %>% select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS") %>% 
    select(ID, visit, GFR, LDL, CRP, UrineSodium, contains("DW"), contains("WW"),
           plasma_acetate = `Acetic acid(µM)`,
           plasma_butyrate = `Butyric acid(µM)`,
           plasma_propionate = `Propionic acid(µM)`,
           plasma_isoval = `Isovaleric acid (µM)`,
           plasma_lactic = `Lactic acid(µM)`,
           plasma_succ = `Succinic acid(µM)`)
urinedata <- readRDS("data/urinesamples.RDS") %>% select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium, Fibers, Energy)
compliance <- readRDS("data/compliance_incl_pharmacy.RDS") %>% select(1:7, 31, 34)
covariates <- right_join(bia, dietarydata, by = c("ID", "visit")) %>% 
    right_join(., compliance, by = c("ID"))

#### Dataset for SCFA analysis ####
df_scfa <- right_join(df, lab, by = c("ID")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           Total_SCFA_DW = DW_AA_umolg + DW_BA_umolg + DW_PA_umolg + DW_FA_umolg) %>% 
    droplevels(.) 

df_scfa$plasma_acetate[which(df_scfa$plasma_acetate > 400)] <- NA

df_means <- df_scfa %>% 
    select(ID, contains("DW"), contains("WW"), contains("plasma"), weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(contains("DW") | contains("WW") | contains("plasma"),
                 list(mean = ~mean(.x, na.rm = TRUE),
                      sd = ~sd(.x, na.rm = TRUE),
                      n = ~length(.x)
                 ),
                 .names = "{.col}_{.fn}"))

#### SCFA plots with LMMs ####
dw_totalscfa_lm <- c()
dw_totalscfa_lm <- df_scfa %>% linearmixed(Total_SCFA_DW)

(plot_total <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 750),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = Total_SCFA_DW_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = Total_SCFA_DW,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = Total_SCFA_DW_mean - (Total_SCFA_DW_sd/sqrt(Total_SCFA_DW_n)),
                          ymax = Total_SCFA_DW_mean + (Total_SCFA_DW_sd/sqrt(Total_SCFA_DW_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        # stat_pvalue_manual(dw_aa, y.position = 300, label = "p_signif", 
        #                    remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,750), breaks = seq(from = 0, to = 750, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "SCFA (umol/g)", title = "Fecal acetate (DW)", color = ""))


dw_acetate_lm <- c()
dw_acetate_lm <- df_scfa %>% linearmixed(DW_AA_umolg)

(plot_acetate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 800),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = DW_AA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = DW_AA_umolg,
        #                                color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_AA_umolg_mean - (DW_AA_umolg_sd/sqrt(DW_AA_umolg_n)),
                          ymax = DW_AA_umolg_mean + (DW_AA_umolg_sd/sqrt(DW_AA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_aa, y.position = 300, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,800), breaks = seq(from = 0, to = 800, by = 100)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Acetic acid (µmol/g)", title = "Fecal acetic acid", color = ""))

dw_butyrate_lm <- c()
dw_butyrate_lm <- df_scfa %>% linearmixed(DW_BA_umolg)

(plot_butyrate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 220),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = DW_BA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = DW_BA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_BA_umolg_mean - (DW_BA_umolg_sd/sqrt(DW_BA_umolg_n)),
                          ymax = DW_BA_umolg_mean + (DW_BA_umolg_sd/sqrt(DW_BA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_ba, y.position = 150, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,220), breaks = seq(from = 0, to = 220, by = 20)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Butyric acid (µmol/g)", title = "Fecal butyric acid",
             color = ""))

dw_propionate_lm <- c()
dw_propionate_lm <- df_scfa %>% linearmixed(DW_PA_umolg)

(plot_propionate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 400),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = DW_PA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = DW_PA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_PA_umolg_mean - (DW_PA_umolg_sd/sqrt(DW_PA_umolg_n)),
                          ymax = DW_PA_umolg_mean + (DW_PA_umolg_sd/sqrt(DW_PA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_pa, y.position = 150, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,400), breaks = seq(from = 0, to = 400, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Proprionic acid (µmol/g)", title = "Fecal propionic acid", 
                color = ""))

dw_fa_lm <- c()
dw_fa_lm <- df_scfa %>% linearmixed(DW_FA_umolg)

(plot_fa <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 30),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = DW_FA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = DW_FA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_FA_umolg_mean - (DW_FA_umolg_sd/sqrt(DW_FA_umolg_n)),
                          ymax = DW_FA_umolg_mean + (DW_FA_umolg_sd/sqrt(DW_FA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_fa_lm, y.position = 25, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,30), breaks = seq(from = 0, to = 30, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Formic acid (µmol/g)", title = "Fecal formic acid (DW)",
    color = ""))



ww_acetate_lm <- c()
ww_acetate_lm <- df_scfa %>% linearmixed(WW_AA_umolg)

(plot_acetate_ww <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 125),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = WW_AA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = WW_AA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = WW_AA_umolg_mean - (WW_AA_umolg_sd/sqrt(WW_AA_umolg_n)),
                          ymax = WW_AA_umolg_mean + (WW_AA_umolg_sd/sqrt(WW_AA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(ww_acetate_lm, y.position = 90, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,125), breaks = seq(from = 0, to = 120, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Acetic acid (µmol/g)", title = "Fecal acetic acid (WW)",
    color = ""))

ww_butyrate_lm <- c()
ww_butyrate_lm <- df_scfa %>% linearmixed(WW_BA_umolg)

(plot_butyrate_ww <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 50),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = WW_BA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = WW_BA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = WW_BA_umolg_mean - (WW_BA_umolg_sd/sqrt(WW_BA_umolg_n)),
                          ymax = WW_BA_umolg_mean + (WW_BA_umolg_sd/sqrt(WW_BA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(ww_butyrate_lm, y.position = 35, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,50), breaks = seq(from = 0, to = 50, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Butyric acid (umol/g)", title = "Fecal butyric acid (WW)",
             color = ""))

ww_propionate_lm <- c()
ww_propionate_lm <- df_scfa %>% linearmixed(WW_PA_umolg)

(plot_propionate_ww <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 40),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = WW_PA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = WW_PA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = WW_PA_umolg_mean - (WW_PA_umolg_sd/sqrt(WW_PA_umolg_n)),
                          ymax = WW_PA_umolg_mean + (WW_PA_umolg_sd/sqrt(WW_PA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(ww_propionate_lm, y.position = 35, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,40), breaks = seq(from = 0, to = 40, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Propionic acid (µmol/g)", title = "Fecal propionic acid (WW)",
             color = ""))

ww_fa_lm <- c()
ww_fa_lm <- df_scfa %>% linearmixed(WW_FA_umolg)

(plot_fa_ww <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 7),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = WW_FA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = WW_FA_umolg,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = WW_FA_umolg_mean - (WW_FA_umolg_sd/sqrt(WW_FA_umolg_n)),
                          ymax = WW_FA_umolg_mean + (WW_FA_umolg_sd/sqrt(WW_FA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(ww_fa_lm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,7), breaks = seq(from = 0, to = 7, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Formic acid (µmol/g)", title = "Fecal formic acid (WW)",
             color = ""))

ggarrange(plot_acetate, plot_butyrate, plot_propionate, plot_fa, 
          nrow = 2, ncol = 2,
          labels = c("A", "B", "C", "D"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/fecalscfa/fecalscfa_dryweightcorr_all.pdf", width = 8, height = 7)
ggsave(filename = "results/fecalscfa/fecalscfa_dryweightcorr_all.svg", width = 8, height = 7)
ggsave(filename = "results/fecalscfa/fecalscfa_dryweightcorr_all.png", width = 8, height = 7)

ggarrange(plot_acetate_ww, plot_butyrate_ww, plot_propionate_ww, plot_fa_ww, 
          nrow = 2, ncol = 2,
          labels = c("A", "B", "C", "D"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/fecalscfa/fecalscfa_wetweightcorr_all.pdf", width = 8, height = 7)
ggsave(filename = "results/fecalscfa/fecalscfa_wetweightcorr_all.svg", width = 8, height = 7)
ggsave(filename = "results/fecalscfa/fecalscfa_wetweightcorr_all.png", width = 8, height = 7)

save_function(plot_acetate, "acetate_dw")
save_function(plot_butyrate, "butyrate_dw")
save_function(plot_propionate, "propionate_dw")
save_function(plot_acetate, "acetate_ww")
save_function(plot_butyrate, "butyrate_ww")
save_function(plot_propionate, "propionate_ww")

#### LMMs with covariates ####
df_scfa_compl <- right_join(covariates, df_scfa, by = c("ID", "visit"))
names(df_scfa_compl)

dw_aa <- linearmixed_scfa_cov(df_scfa_compl, DW_AA_umolg)
dw_ba <- linearmixed_scfa_cov(df_scfa_compl, DW_BA_umolg)
dw_pa <- linearmixed_scfa_cov(df_scfa_compl, DW_PA_umolg)
# linearmixed_scfa_cov(df_scfa_compl, DW_FA_umolg)

# linearmixed_scfa_cov(df_scfa_compl, WW_AA_umolg)
# linearmixed_scfa_cov(df_scfa_compl, WW_BA_umolg)
# linearmixed_scfa_cov(df_scfa_compl, WW_PA_umolg)
# linearmixed_scfa_cov(df_scfa_compl, WW_FA_umolg)

p_acetate_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_acetate)
p_butyrate_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_butyrate)
p_propionate_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_propionate)
p_isoval_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_isoval)
p_lactic_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_lactic)
p_succ_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_succ)

#### Plasma SCFA ####
(plot_p_acetate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 200),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_acetate_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = plasma_acetate,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_acetate_mean - (plasma_acetate_sd/sqrt(plasma_acetate_n)),
                          ymax = plasma_acetate_mean + (plasma_acetate_sd/sqrt(plasma_acetate_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_acetate_lm, y.position = 150, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,200), breaks = seq(from = 0, to = 200, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Acetic acid (µM)", title = "Plasma acetic acid", color = ""))

(plot_p_butyrate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 8),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_butyrate_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = plasma_butyrate,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_butyrate_mean - (plasma_butyrate_sd/sqrt(plasma_butyrate_n)),
                          ymax = plasma_butyrate_mean + (plasma_butyrate_sd/sqrt(plasma_butyrate_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_butyrate_lm, y.position = 6, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,8), breaks = seq(from = 0, to = 8, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Butyric acid (µM)", title = "Plasma butyric acid", color = ""))

(plot_p_propionate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 10),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_propionate_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = plasma_propionate,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_propionate_mean - (plasma_propionate_sd/sqrt(plasma_propionate_n)),
                          ymax = plasma_propionate_mean + (plasma_propionate_sd/sqrt(plasma_propionate_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_propionate_lm, y.position = 6, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,10), breaks = seq(from = 0, to = 10, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Propionic acid (µM)", title = "Plasma propionic acid", color = ""))

(plot_p_isoval <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 4),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_isoval_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = plasma_isoval,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_isoval_mean - (plasma_isoval_sd/sqrt(plasma_isoval_n)),
                          ymax = plasma_isoval_mean + (plasma_isoval_sd/sqrt(plasma_isoval_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(p_isoval_lm, y.position = 3.5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,4), breaks = seq(from = 0, to = 4, by = 0.5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Isovaleric acid (µM)", title = "Plasma isovaleric acid", color = ""))

(plot_p_lactic <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 500, ymax = 1750),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_lactic_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = plasma_lactic,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_lactic_mean - (plasma_lactic_sd/sqrt(plasma_lactic_n)),
                          ymax = plasma_lactic_mean + (plasma_lactic_sd/sqrt(plasma_lactic_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_lactic_lm, y.position = 1500, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(500,1750), breaks = seq(from = 500, to = 1750, by = 250)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Lactic acid (µM)", title = "Plasma lactic acid", color = ""))

(plot_p_succ <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 16),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_succ_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.75) +
        # geom_line(data = df_scfa, aes(x = weeks, y = plasma_succ,
        #                               color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_succ_mean - (plasma_succ_sd/sqrt(plasma_succ_n)),
                          ymax = plasma_succ_mean + (plasma_succ_sd/sqrt(plasma_succ_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_succ_lm, y.position = 3, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,16), breaks = seq(from = 0, to = 16, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Succinic acid (µM)", title = "Plasma succinnic acid", color = ""))

ggarrange(plot_p_acetate, plot_p_butyrate, plot_p_propionate, 
          plot_p_isoval, plot_p_lactic, plot_p_succ, 
          nrow = 2, ncol = 3,
          labels = c("A", "B", "C", "D", "E", "F"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/plasmascfa/plasmascfa_all.pdf", width = 10, height = 7)
ggsave(filename = "results/plasmascfa/plasmascfa_all.svg", width = 10, height = 7)
ggsave(filename = "results/plasmascfa/plasmascfa_all.png", width = 10, height = 7)

ggarrange(plot_p_acetate, plot_p_butyrate, plot_p_propionate,
          plot_acetate, plot_butyrate, plot_propionate,
          nrow = 2, ncol = 3,
          labels = c("A", "B", "C", "D", "E", "F"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/plasmascfa/plasma_fecalscfa.pdf", width = 10, height = 7)
ggsave(filename = "results/plasmascfa/plasma_fecalscfa.svg", width = 10, height = 7)
ggsave(filename = "results/plasmascfa/plasma_fecalscfa.png", width = 10, height = 7)

save_function_plasma(plot_p_acetate, "plasma_acetic")
save_function_plasma(plot_p_butyrate, "plasma_butyric")
save_function_plasma(plot_p_propionate, "plasma_propionic")
save_function_plasma(plot_p_lactic, "plasma_lactic")
save_function_plasma(plot_p_succ, "plasma_succ")
save_function_plasma(plot_p_isoval, "plasma_isoval")

#### LMMs with covariates ####
linearmixed_scfa_cov(df_scfa_compl, plasma_acetate)
linearmixed_scfa_cov(df_scfa_compl, plasma_butyrate)
linearmixed_scfa_cov(df_scfa_compl, plasma_propionate)


