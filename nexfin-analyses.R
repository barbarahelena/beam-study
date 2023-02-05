# Nexfin analysis
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

linearmixed_nexfin <- function(data, var){
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
    statres$pval <- case_when(statres$pval < 0.1 ~ paste0(str_c("p = ", statres$pval)),
                              statres$pval > 0.1 ~ paste0(""))
    return(statres)
}

save_function_nexfin <- function(plot, name, width = 5, height = 4){
    ggsave(plot = plot, 
           filename = str_c("results/nexfin/", name, ".pdf"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/nexfin/", name, ".svg"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/nexfin/", name, ".png"), width = width, height = height)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
nexfin <- readRDS("data/nexfin_data.RDS") 
bia <- readRDS("data/bia_data.RDS") %>% select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS") %>% 
    select(ID, visit, GFR, Na, K, Kreat, UrineSodium, UrineK, UrineCreat)
officebp <- readRDS("data/officebp_summary.RDS") %>% select(ID, visit, Systolic, Diastolic)
dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium, Fibers)
covariates <- right_join(bia, dietarydata, by = c("ID", "visit")) %>% 
    right_join(officebp, ., by = c("ID", "visit")) %>% 
    right_join(df, ., by = c("ID")) %>% 
    filter(!ID %in% c("BEAM_713", "BEAM_664", "BEAM_299")) %>% 
    filter(visit != "V3")
names(covariates)
tail(covariates)
nexfin_total <- right_join(nexfin, covariates, by = c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_713", "BEAM_664", "BEAM_299"))
head(nexfin_total)

#### Prepare data ####
nexfin_total <- nexfin_total %>% 
    mutate(weeks = case_when(
        visit == "V2" ~ paste0(0),
        visit == "V4" ~ paste0(4),
        visit == "V5" ~ paste0(5)
    ),
    weeks = as.numeric(weeks),
    before_after = case_when(
        visit == "V2" ~ paste0("Before"),
        visit == "V4" ~ paste0("Treatment"),
        visit == "V5" ~ paste0("After")
    ),
    before_after = fct_relevel(before_after, "Before", after = 0L),
    before_after = fct_relevel(before_after, "After", after = 2L))

nexfin_means <- nexfin_total %>% 
    select(ID, MAP, SBP, DBP, AvgRR, SV, CO, SVR, meanBRS, SDNN, RMSDD, NN50, pNN50, dPdt,
           weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c(MAP, SBP, DBP, SV, CO, SVR, meanBRS, SDNN, RMSDD, NN50, pNN50, dPdt), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                          ),
                     .names = "{.col}_{.fn}"))


sbp_lmm <- linearmixed_nexfin(nexfin_total, SBP)
dbp_lmm <- linearmixed_nexfin(nexfin_total, DBP)
map_lmm <- linearmixed_nexfin(nexfin_total, MAP)
co_lmm <- linearmixed_nexfin(nexfin_total, CO)
dpdt_lmm <- linearmixed_nexfin(nexfin_total, dPdt)
sdnn_lmm <- linearmixed_nexfin(nexfin_total, SDNN)
pnn50_lmm <- linearmixed_nexfin(nexfin_total, pNN50)
dpdt_lmm <- linearmixed_nexfin(nexfin_total, dPdt)
rmsdd_lmm <- linearmixed_nexfin(nexfin_total, RMSDD)
brs_lmm <- linearmixed_nexfin(nexfin_total, meanBRS)

(plot_sbp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 200),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = SBP_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = SBP,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = SBP_mean - (SBP_sd/sqrt(SBP_n)),
                          ymax = SBP_mean + (SBP_sd/sqrt(SBP_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(sbp_lmm, y.position = 170, label =  "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(100, 200), breaks = seq(from = 100, to = 200, by = 10)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Systolic BP (mmHg)", title = "Systolic BP",
             color = ""))

(plot_dbp <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 55, ymax = 100),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = DBP_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = DBP,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = DBP_mean - (DBP_sd/sqrt(DBP_n)),
                          ymax = DBP_mean + (DBP_sd/sqrt(DBP_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(dbp_lmm, y.position = 90, label =  "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(55, 100), breaks = seq(from = 55, to = 100, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Diastolic BP (mmHg)", title = "Diastolic BP",
             color = ""))

(plot_co <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 3, ymax = 7.5),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = CO_mean, 
                                           color = Treatment_group, 
                                           group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = CO,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = CO_mean - (CO_sd/sqrt(CO_n)),
                          ymax = CO_mean + (CO_sd/sqrt(CO_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(co_lmm, y.position = 7.0, label = "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(3,7.5), breaks = seq(from = 3, to = 7.5, by = 0.5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Cardiac output (L/min)", title = "Cardiac output",
             color = ""))

(plot_dpdt<- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 3000),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = dPdt_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = dPdt,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = dPdt_mean - (dPdt_sd/sqrt(MAP_n)),
                          ymax = dPdt_mean + (dPdt_sd/sqrt(MAP_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(dpdt_lmm, y.position = 2000, label = "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,3000), breaks = seq(from = 0, to = 3000, by = 500)) +
        theme_Publication() +
        labs(x = "Weeks", y = "dPdt", title = "dPdt",
             color = ""))

(plot_pnn50<- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 0.3),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = pNN50_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = pNN50,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = pNN50_mean - (pNN50_sd/sqrt(pNN50_n)),
                          ymax = pNN50_mean + (pNN50_sd/sqrt(pNN50_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(pnn50_lmm, y.position = 0.2, label = "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,0.3), breaks = seq(from = 0, to = 0.3, by = 0.1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "pNN50 (proportion)", title = "pNN50",
             color = ""))

(plot_sdnn<- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 0.10),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = SDNN_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = SDNN,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = SDNN_mean - (SDNN_sd/sqrt(SDNN_n)),
                          ymax = SDNN_mean + (SDNN_sd/sqrt(SDNN_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(sdnn_lmm, y.position = 0.08, label = "{pval}",
                           remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,0.10), breaks = seq(from = 0, to = 0.10, by = 0.02)) +
        theme_Publication() +
        labs(x = "Weeks", y = "SDNN", title = "Heart rate variability",
             color = ""))

(plot_brs <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 16),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = nexfin_means, aes(x = weeks, y = meanBRS_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = nexfin_total, aes(x = weeks, y = meanBRS,
                                           color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = nexfin_means,
                      aes(ymin = meanBRS_mean - (meanBRS_sd/sqrt(meanBRS_n)),
                          ymax = meanBRS_mean + (meanBRS_sd/sqrt(meanBRS_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(brs_lmm, y.position = 13, label = "{pval}",
                           remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,16), breaks = seq(from = 0, to = 16, by = 2)) +
        theme_Publication() +
        labs(x = "Weeks", y = "BRS", title = "Baroreceptor sensitivity",
             color = ""))

(pl_nexfin <- ggarrange(plot_sbp, plot_dbp, plot_co, plot_dpdt, plot_brs, plot_sdnn, 
          labels = c("A", "B", "C", "D", "E", "F"),
          nrow = 2, ncol = 3,
          legend = "bottom",
          common.legend = TRUE))
save_function_nexfin(pl_nexfin, "nexfin_plots", width = 9, height = 6)

save_function_nexfin(plot_sbp, "sbp_nexfin")
save_function_nexfin(plot_dbp, "dbp_nexfin")
save_function_nexfin(plot_co, "co_nexfin")
save_function_nexfin(plot_co, "dpdt_nexfin")
save_function_nexfin(plot_co, "brs_nexfin")
save_function_nexfin(plot_sdnn, "sdnn_nexfin")
save_function_nexfin(plot_pnn50, "pnn50_nexfin")
save_function_nexfin(plot_rmsdd, "rmsdd_nexfin")
