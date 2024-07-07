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
    statres$p.signif <- case_when(
        statres$pval < 0.05 ~paste0("*"),
        statres$pval < 0.01 ~paste0("**"),
        statres$pval < 0.001 ~paste0("***"),
        statres$pval > 0.05 ~paste0("")
    )
    statres <- statres %>% filter(p.signif != "")
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
bia <- readRDS("data/bia_data.RDS") %>% dplyr::select(ID, visit, BMI)
lab <- readRDS("data/lab_results.RDS")  %>% 
    dplyr::select(ID, visit, GFR, Na, K, Kreat, UrineSodium, UrineK, UrineCreat)
officebp <- readRDS("data/officebp_summary.RDS") %>% 
    dplyr::select(ID, visit, Systolic, Diastolic)
urinedata <- readRDS("data/urinesamples.RDS") %>% 
    dplyr::select(ID, visit, Volume)
dietarydata <- readRDS("data/diet_summary.RDS") %>% 
    dplyr::select(ID, visit, Sodium, Fibers)
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
        geom_line(data = urine_total, aes(x = weeks, y = FENa,
                                      color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = urine_total, aes(x = weeks, y = FENa,
                                       color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = urine_means, aes(x = weeks, y = FENa_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = urine_means, aes(x = weeks, y = FENa_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = urine_means,
                      aes(ymin = FENa_mean - (FENa_sd/sqrt(11)),
                          ymax = FENa_mean + (FENa_sd/sqrt(11)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(fena_lmm, y.position = 5, label = "p.signif", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,1.5), breaks = seq(from = 0, to = 1.5, by = 0.25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "FENa (%)", title = "Fractional sodium excretion",
             color = ""))

