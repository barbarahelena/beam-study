# Fractional excretion analyses
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

linearmixed_excr <- function(data, var){
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

linearmixed_excr_cov <- function(data, var){
    data1 <- data %>% filter(weeks %in% c(0,4)) %>% 
        mutate(var = {{ var }})
    model1_v4 <- lmer(var ~ Treatment_group*weeks + (1|ID) + Fibers + BMI + V1_Systolic, 
                      data = data1)
    res_v4 <- summary(model1_v4)
    print(res_v4)
    pval <- format(round(res_v4$coefficients[7,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    statres_line1 <- cbind(group1 = 0, group2 = 4, pval)
    
    data2 <- data %>%
        filter(weeks %in% c(4,5)) %>% 
        mutate(var = {{ var }})
    model1_v5 <- lmer(var ~ Treatment_group*weeks + (1|ID) + Fibers + BMI + V1_Systolic,
                      data = data2)
    res_v5 <- summary(model1_v5)
    print(res_v5)
    pval <- format(round(res_v5$coefficients[7,5], 3), nsmall = 3)
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

save_function <- function(plot, name, width = 5, height = 4){
    ggsave(plot = plot, 
           filename = str_c("results/frac_excretion/", name, ".pdf"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/frac_excretion/", name, ".svg"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/frac_excretion/", name, ".png"), width = width, height = height)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
bia <- readRDS("data/bia_data.RDS") %>% select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS")  %>% 
    dplyr::select(ID, visit, GFR, Na, K, Kreat, UrineSodium, UrineK, UrineCreat)
officebp <- readRDS("data/officebp_summary.RDS") %>% select(ID, visit, Systolic, Diastolic)
urinedata <- readRDS("data/urinesamples.RDS") %>% select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium, Fibers)
covariates <- right_join(bia, dietarydata, by = c("ID", "visit")) %>% 
    right_join(officebp, ., by = c("ID", "visit")) %>% 
    right_join(df, ., by = c("ID")) %>% 
    filter(!ID %in% c("BEAM_713", "BEAM_664", "BEAM_299")) %>% 
    filter(visit != "V3")
names(covariates)
tail(covariates)
urinelab <- right_join(lab, urinedata, by = c("ID", "visit"))
urine_total <- right_join(urinelab, covariates, by = c("ID", "visit")) %>% 
    filter(!ID %in% c("BEAM_713", "BEAM_664", "BEAM_299"))
# Replace missing volume for BEAM_862 by average of V4 and V5
urine_total$Volume[which(urine_total$ID == "BEAM_862") & urine_total$visit == "V2"] <- 3075
names(urine_total)
head(urine_total)

#### Calculations fractional excretion ####
# gfr = urine kreat x urine vol / plasma kreat
# fr na excr : Natrium-urine/Natrium-plasma x (Kreatinine-plasma /1000)/Kreatinine-urine x 100
# fr K excr : K-urine/K-plasma x (Kreatinine-plasma /1000)/Kreatinine-urine x 100
head(urinelab)
urine_total <- urine_total %>% mutate(
    BSA = 0.007184 *  (V1_Weight ^ 0.425) * (V1_Height ^ 0.725),
    Urine_Na_Total = UrineSodium * (Volume / 1000),
    Urine_K_Total = UrineK * (Volume / 1000),
    Urine_Kreat_Total = UrineCreat * (Volume / 1000),
    Urine_Na_L = UrineSodium / (Volume / 1000),
    Urine_K_L = UrineK / (Volume / 1000),
    Urine_Kreat_L = UrineCreat / (Volume / 1000),
    FENa = ((Urine_Na_L * (Kreat/1000)) / ((Urine_Kreat_L) * Na)) * 100,
    FEK = ((Urine_K_L * (Kreat/1000)) / ((Urine_Kreat_L)* K)) * 100,
    Creatclearance =  ((Urine_Kreat_L * Volume) * 1000 / (Kreat *1440)),
    visit = as.factor(visit), 
    weeks = case_when(
        visit == "V2" ~ paste0(0),
        visit == "V4" ~ paste0(4),
        visit == "V5" ~ paste0(5)
    ),
    weeks = as.numeric(weeks),
    Treatment_group = as.factor(Treatment_group)
) %>% 
    droplevels(.)

#### Data / calculation checks ####
hist(urine_total$BSA) # body surface area
plot(urine_total$BSA, urine_total$V1_Height)
plot(urine_total$BSA, urine_total$V1_Weight)
hist(urine_total$Urine_Na_L)
hist(urine_total$Urine_K_L)
hist(urine_total$Urine_Kreat_L)
hist(urine_total$FENa)
hist(urine_total$FEK)
hist(urine_total$Creatclearance)
plot(urine_total$Creatclearance, urine_total$GFR)
plot(urine_total$Sodium, urine_total$Urine_Na_Total)
plot(urine_total$FENa, urine_total$Urine_Na_Total)
plot(urine_total$FENa, urine_total$Creatclearance)

#### Calculating means ####
urine_means <- urine_total %>% 
    dplyr::select(ID, FENa, FEK, Creatclearance, Urine_Na_L, Urine_Kreat_L, GFR, Na,
           weeks, Treatment_group) %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c(FENa, FEK, Creatclearance, Urine_Na_L, Urine_Kreat_L, GFR, Na), 
                     list(mean = mean, sd = sd), .names = "{.col}_{.fn}"))
urine_means

#### LMM and plot ####
fena_lmm <- linearmixed_excr(urine_total, FENa)

(plot_fena <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 1.5),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = FENa_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = FENa,
                                      color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = FENa_mean - (FENa_sd/sqrt(11)),
                          ymax = FENa_mean + (FENa_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(fena_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,1.5), breaks = seq(from = 0, to = 1.5, by = 0.25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "FENa (%)", title = "Fractional sodium excretion",
             color = ""))

fek_lmm <- linearmixed_excr(urine_total, FEK)

(plot_fek <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 18),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = FEK_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = FEK,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = FEK_mean - (FEK_sd/sqrt(11)),
                          ymax = FEK_mean + (FEK_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(fek_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,18), breaks = seq(from = 0, to = 18, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "FEK (%)", title = "Fractional potassium excretion",
             color = ""))

urinena_lmm <- linearmixed_excr(urine_total, Urine_Na_Total)

(plot_na <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 300),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = Urine_Na_L_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = Urine_Na_L,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = Urine_Na_L_mean - (Urine_Na_L_sd/sqrt(11)),
                          ymax = Urine_Na_L_mean + (Urine_Na_L_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(urinena_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,300), breaks = seq(from = 0, to = 300, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Urine sodium (mmol/L)", title = "Urine sodium",
             color = ""))

urinekreat_lmm <- linearmixed_excr(urine_total, Urine_Kreat_L)

(plot_totalkreat <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 30),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = Urine_Kreat_L_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = Urine_Kreat_L,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = Urine_Kreat_L_mean - (Urine_Kreat_L_sd/sqrt(11)),
                          ymax = Urine_Kreat_L_mean + (Urine_Kreat_L_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(urinekreat_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,30), breaks = seq(from = 0, to = 30, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Urine creatinine (mmol/L)", title = "Total urine creatinine",
             color = ""))


urinegfr_lmm <- linearmixed_excr(urine_total, Creatclearance)

(plot_totalgfr <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 150),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = Creatclearance_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = Creatclearance,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = Creatclearance_mean - (Creatclearance_sd/sqrt(11)),
                          ymax = Creatclearance_mean + (Creatclearance_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(urinegfr_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,150), breaks = seq(from = 0, to = 150, by = 25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Creatinine clearance (ml/min)", title = "Creatinine clearance",
             color = ""))

gfr_lmm <- linearmixed_excr(urine_total, GFR)

(plot_gfr <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 60, ymax = 90),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = GFR_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = GFR,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = GFR_mean - (GFR_sd/sqrt(11)),
                          ymax = GFR_mean + (GFR_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(gfr_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(60,90), breaks = seq(from = 60, to = 90, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "eGFR (ml/min/1.73)", title = "eGFR",
             color = ""))


plasmana_lmm <- linearmixed_excr(urine_total, Na)

(plot_plasmana <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 135, ymax = 145),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = urine_means, aes(x = weeks, y = Na_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = urine_total, aes(x = weeks, y = Na,
                                          color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = urine_means,
                      aes(ymin = Na_mean - (Na_sd/sqrt(11)),
                          ymax = Na_mean + (Na_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(plasmana_lmm, y.position = 5, label = "p_signif", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(135, 145), breaks = seq(from = 135, to = 145, by = 1)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Plasma sodium (umol/L)", title = "Plasma sodium",
             color = ""))

pl_natrium <- ggarrange(plot_na, plot_plasmana, plot_fena, labels = c("A", "B", "C"),
          nrow = 1, ncol = 3,
          common.legend = TRUE,
          legend = "bottom")
save_function(pl_natrium, "sodium_plasmaurine", width = 10, height = 4)

save_function(plot_totalgfr, "creatclearance")
save_function(plot_fena, "fract_na_excr")
save_function(plot_fek, "fract_k_excr")
save_function(plot_na, "urine_sodium")


#### Boxplots before after intervention ####
urine_total %>% filter(visit %in% c("V2", "V4")) %>% 
    mutate(before_after = case_when(
    visit == "V2" ~ paste0("Before"),
    visit == "V4" ~ paste0("After")
),
    before_after = fct_relevel(before_after, "Before", after = 0L)) %>% 
    ggplot(., aes(x = before_after, y = FENa)) +
    geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
    geom_point(size = 0.75)+
    geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
    stat_compare_means(method = "wilcox.test") +
    facet_wrap(~Treatment_group) + 
    scale_fill_jama(guide = "none") +
    scale_color_jama(guide = "none") +
    labs(x = "", title = "Differences in FENa", y = "FENa (%)") +
    theme_Publication() 

urine_boxplots <- urine_total %>% 
    mutate(before_after = case_when(
        visit == "V2" ~ paste0("Before"),
        visit == "V4" ~ paste0("Treatment"),
        visit == "V5" ~ paste0("After")
    ),
    before_after = fct_relevel(before_after, "Before", after = 0L),
    before_after = fct_relevel(before_after, "After", after = 2L))

(boxplot_fena <- ggplot(urine_boxplots, aes(x = before_after, y = FENa)) +
    geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
    geom_point(size = 0.75)+
    geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
    stat_compare_means(method = "wilcox.test", label = "p.signif",
                       comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
    facet_wrap(~Treatment_group) + 
    scale_fill_jama(guide = "none") +
    scale_color_jama(guide = "none") +
    labs(x = "", title = "Effect on FENa", y = "FENa (%)") +
    theme_Publication() )
save_function(boxplot_fena, "boxplot_fena", height = 5)

(boxplot_fek <- ggplot(urine_boxplots, aes(x = before_after, y = FEK)) +
        geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        stat_compare_means(method = "wilcox.test", label = "p.signif",
                           comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
        facet_wrap(~Treatment_group) + 
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(x = "", title = "Effect on FEK", y = "FEK (%)") +
        theme_Publication() )
save_function(boxplot_fek, "boxplot_fek", height = 5)

(boxplot_cc <- ggplot(urine_boxplots, aes(x = before_after, y = Creatclearance)) +
        geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        stat_compare_means(method = "wilcox.test", label = "p.signif",
                           comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
        facet_wrap(~Treatment_group) + 
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(x = "", title = "Effect on creatinine clearance", y = "Creatinine clearance (ml/min)") +
        theme_Publication() )
save_function(boxplot_cc, "boxplot_cc", height = 5)
