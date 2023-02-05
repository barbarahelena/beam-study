# Lab results analysis
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

linearmixed_lab <- function(data, var){
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

save_function_lab <- function(plot, name, a = 5, b = 4){
    ggsave(plot = plot, 
           filename = str_c("results/lab/", name, ".pdf"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/lab/", name, ".svg"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/lab/", name, ".png"), width = a, height = b)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
# bia <- readRDS("data/bia_data.RDS")
lab <- readRDS("data/lab_results.RDS")
#     select(ID, visit, GFR, LDL, CRP, UrineSodium, contains("DW"), contains("WW"))
urinedata <- readRDS("data/urinesamples.RDS") %>% select(ID, visit, Volume)
# dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium, Fibers)
covariates <- right_join(urinedata, df, by = "ID")

#### Dataset for SCFA analysis ####
df_lab <- right_join(lab, covariates, by = c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks)
           ) %>% 
    droplevels(.) 

df_means <- df_lab %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("CRP", "Hb", "Ht", "Trombo", "Leuko",
                       "Na", "K", "Kreat", "GFR", "TC", "HDL", "LDL", "TG"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

#### Lab plots with LMMs ####
crp_lm <- df_lab %>% linearmixed_bia(CRP)
hb_lm <- df_lab %>% linearmixed_bia(Hb)
ht_lm <- df_lab %>% linearmixed_bia(Ht)
trombo_lm <- df_lab %>% linearmixed_bia(Trombo)
leuko_lm <- df_lab %>% linearmixed_bia(Leuko)
na_lm <- df_lab %>% linearmixed_bia(Na)
k_lm <- df_lab %>% linearmixed_bia(K)
kreat_lm <- df_lab %>% linearmixed_bia(Kreat)
gfr_lm <- df_lab %>% linearmixed_bia(GFR)
tc_lm <- df_lab %>% linearmixed_bia(TC)
hdl_lm <- df_lab %>% linearmixed_bia(HDL)
ldl_lm <- df_lab %>% linearmixed_bia(LDL)
tg_lm <- df_lab %>% linearmixed_bia(TG)

(plot_tc <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 3, ymax = 8),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = TC_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_lab, aes(x = weeks, y = TC,
                                      color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = TC_mean - (TC_sd/sqrt(TC_n)),
                          ymax = TC_mean + (TC_sd/sqrt(TC_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(tc_lm, y.position = 6, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        #scale_y_continuous(limits = c(18,30), breaks = seq(from = 18, to = 30, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Total cholesterol (mmol/L)", title = "TC", color = ""))

(plot_ldl <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 1, ymax = 6),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = LDL_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_lab, aes(x = weeks, y = LDL,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = LDL_mean - (LDL_sd/sqrt(LDL_n)),
                          ymax = LDL_mean + (LDL_sd/sqrt(LDL_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(tc_lm, y.position = 5, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        #scale_y_continuous(limits = c(18,30), breaks = seq(from = 18, to = 30, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "LDL (mmol/L)", title = "LDL", color = ""))

(plot_hdl <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0.5, ymax = 3),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = HDL_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_lab, aes(x = weeks, y = HDL,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = HDL_mean - (HDL_sd/sqrt(HDL_n)),
                          ymax = HDL_mean + (HDL_sd/sqrt(HDL_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(tc_lm, y.position = 2.5, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        #scale_y_continuous(limits = c(18,30), breaks = seq(from = 18, to = 30, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "HDL (mmol/L)", title = "HDL", color = ""))

(plot_tg <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 2.5),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_means, aes(x = weeks, y = TG_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_lab, aes(x = weeks, y = TG,
                                     color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = df_means,
                      aes(ymin = TG_mean - (TG_sd/sqrt(TG_n)),
                          ymax = TG_mean + (TG_sd/sqrt(TG_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(tg_lm, y.position = 2.0, label = "p_signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        #scale_y_continuous(limits = c(18,30), breaks = seq(from = 18, to = 30, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Triglycerides (mmol/L)", title = "TG", color = ""))

(pl_lab <- ggarrange(plot_tc, plot_ldl, plot_hdl, plot_tg, 
                        labels = c("A", "B", "C", "D"),
                        nrow = 2, ncol = 2, common.legend = TRUE,
                     legend = "bottom")
    )
save_function_lab(pl_lab, "cholesterol_plots", a = 7, b = 6)

