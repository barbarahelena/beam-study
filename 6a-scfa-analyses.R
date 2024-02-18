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
    statres$p.signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        statres$pval > 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
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

df_scfa$plasma_acetate[which(df_scfa$plasma_acetate > 400)] <- NA # remove one extreme outlier

df_scfa_compl <- right_join(covariates, df_scfa, by = c("ID", "visit"))

df_means <- df_scfa_compl %>% 
    select(ID, contains("DW"), contains("WW"), contains("plasma"), weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(contains("DW") | contains("WW") | contains("plasma"),
                 list(mean = ~mean(.x, na.rm = TRUE),
                      sd = ~sd(.x, na.rm = TRUE),
                      n = ~length(.x)
                 ),
                 .names = "{.col}_{.fn}"))

#### SCFA plots with LMMs ####
dw_aa <- linearmixed_scfa_cov(df_scfa_compl, DW_AA_umolg)
dw_ba <- linearmixed_scfa_cov(df_scfa_compl, DW_BA_umolg)
dw_pa <- linearmixed_scfa_cov(df_scfa_compl, DW_PA_umolg)

(plot_acetate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 800),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_scfa, aes(x = weeks, y = DW_AA_umolg,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_scfa, aes(x = weeks, y = DW_AA_umolg,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = DW_AA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = DW_AA_umolg_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_AA_umolg_mean - (DW_AA_umolg_sd/sqrt(DW_AA_umolg_n)),
                          ymax = DW_AA_umolg_mean + (DW_AA_umolg_sd/sqrt(DW_AA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_aa, y.position = 300, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,800), breaks = seq(from = 0, to = 800, by = 100)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Acetic acid (µmol/g)", title = "Fecal acetic acid", color = ""))

(plot_butyrate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 220),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_scfa, aes(x = weeks, y = DW_BA_umolg,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_scfa, aes(x = weeks, y = DW_BA_umolg,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = DW_BA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = DW_BA_umolg_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_BA_umolg_mean - (DW_BA_umolg_sd/sqrt(DW_BA_umolg_n)),
                          ymax = DW_BA_umolg_mean + (DW_BA_umolg_sd/sqrt(DW_BA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_ba, y.position = 150, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,220), breaks = seq(from = 0, to = 220, by = 20)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Butyric acid (µmol/g)", title = "Fecal butyric acid",
             color = ""))

(plot_propionate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 400),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_scfa, aes(x = weeks, y = DW_PA_umolg,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_scfa, aes(x = weeks, y = DW_PA_umolg,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = DW_PA_umolg_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = DW_PA_umolg_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = DW_PA_umolg_mean - (DW_PA_umolg_sd/sqrt(DW_PA_umolg_n)),
                          ymax = DW_PA_umolg_mean + (DW_PA_umolg_sd/sqrt(DW_PA_umolg_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(dw_pa, y.position = 150, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,400), breaks = seq(from = 0, to = 400, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Proprionic acid (µmol/g)", title = "Fecal propionic acid", 
                color = ""))

save_function(plot_acetate, "acetate_dw")
save_function(plot_butyrate, "butyrate_dw")
save_function(plot_propionate, "propionate_dw")


#### Plasma SCFA ####
p_acetate_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_acetate)
p_butyrate_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_butyrate)
p_propionate_lm <- df_scfa_compl %>% linearmixed_scfa_cov(plasma_propionate)

(plot_p_acetate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 200),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_scfa, aes(x = weeks, y = plasma_acetate,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_scfa, aes(x = weeks, y = plasma_acetate,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_acetate_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = plasma_acetate_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_acetate_mean - (plasma_acetate_sd/sqrt(plasma_acetate_n)),
                          ymax = plasma_acetate_mean + (plasma_acetate_sd/sqrt(plasma_acetate_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_acetate_lm, y.position = 150, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,200), breaks = seq(from = 0, to = 200, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Acetic acid (µM)", title = "Plasma acetic acid", color = ""))

(plot_p_butyrate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 8),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_scfa, aes(x = weeks, y = plasma_butyrate,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_scfa, aes(x = weeks, y = plasma_butyrate,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_butyrate_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = plasma_butyrate_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_butyrate_mean - (plasma_butyrate_sd/sqrt(plasma_butyrate_n)),
                          ymax = plasma_butyrate_mean + (plasma_butyrate_sd/sqrt(plasma_butyrate_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_butyrate_lm, y.position = 6, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,8), breaks = seq(from = 0, to = 8, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Butyric acid (µM)", title = "Plasma butyric acid", color = ""))

(plot_p_propionate <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 10),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_scfa, aes(x = weeks, y = plasma_propionate,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_scfa, aes(x = weeks, y = plasma_propionate,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = plasma_propionate_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = plasma_propionate_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = plasma_propionate_mean - (plasma_propionate_sd/sqrt(plasma_propionate_n)),
                          ymax = plasma_propionate_mean + (plasma_propionate_sd/sqrt(plasma_propionate_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(p_propionate_lm, y.position = 6, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,10), breaks = seq(from = 0, to = 10, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Propionic acid (µM)", title = "Plasma propionic acid", color = ""))

ggarrange(plot_p_acetate, plot_p_butyrate, plot_p_propionate,
          plot_acetate, plot_butyrate, plot_propionate,
          nrow = 2, ncol = 3,
          labels = c("A", "B", "C", "D", "E", "F"),
          common.legend = TRUE,
          legend = "bottom")
ggsave(filename = "results/plasmascfa/plasma_fecalscfa.pdf", width = 10, height = 7)
ggsave(filename = "results/plasmascfa/plasma_fecalscfa.svg", width = 10, height = 7)
ggsave(filename = "results/plasmascfa/plasma_fecalscfa.png", width = 10, height = 7)

