# Monocytes results analyses
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

linearmixed_mono_elisa <- function(data, var){
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

linearmixed_mono <- function(data, var){
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

function_plotmono <- function(data, data_means, 
                              var, var_mean, var_sd, var_n,
                              title, yname,
                              var_lm, positiony) {
    data1 <- data %>% mutate(var = {{ var }})
    max_var <- max(data1$var)
    min_var <- min(data1$var)
    data_means1 <- data_means %>% mutate(var_mean = {{ var_mean }},
                                         var_sd = {{ var_sd }},
                                         var_n = {{ var_n }})
    (ggplot() +
            geom_rect(aes(xmin = 0, xmax = 4, ymin = min_var, ymax = max_var),
                      fill = "#CDCDCD", alpha = 0.3) +
            geom_line(data = data1, aes(x = weeks, y = var,
                                          color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
            geom_point(data = data1, aes(x = weeks, y = var,
                                           color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
            geom_line(data = data_means1, aes(x = weeks, y = var_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
            geom_point(data = data_means1, aes(x = weeks, y = var_mean, 
                                            color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
            geom_errorbar(data = data_means1,
                          aes(ymin = var_mean - (var_sd/sqrt(var_n)),
                              ymax = var_mean + (var_sd/sqrt(var_n)),
                              x = weeks,
                              color = Treatment_group), width=0.1) +
            stat_pvalue_manual(var_lm, y.position = positiony, label = "p.signif",
                               remove.bracket = TRUE, bracket.size = 0) +
            scale_color_jama() + 
            scale_y_continuous(limits = c(min_var,max_var)) +
            theme_Publication() +
            labs(x = "Weeks", y = yname, title = title, color = ""))
}

function_plotelisa <- function(data, data_means, 
                               var, var_mean, var_sd, var_n,
                               title, yname,
                               var_lm, positiony) {
    data1 <- data %>% mutate(var = {{ var }})
    max_var <- max(data1$var)*1.2
    min_var <- min(data1$var)*1.2
    data_means1 <- data_means %>% mutate(var_mean = {{ var_mean }},
                                         var_sd = {{ var_sd }},
                                         var_n = {{ var_n }})
    (ggplot() +
            geom_rect(aes(xmin = 0, xmax = 4, ymin = min_var, ymax = max_var),
                      fill = "#CDCDCD", alpha = 0.3) +
            geom_line(data = data1, aes(x = weeks, y = var,
                                        color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
            geom_point(data = data1, aes(x = weeks, y = var,
                                         color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
            geom_line(data = data_means1, aes(x = weeks, y = var_mean, 
                                              color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
            geom_point(data = data_means1, aes(x = weeks, y = var_mean, 
                                               color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
            geom_errorbar(data = data_means1,
                          aes(ymin = var_mean - (var_sd/sqrt(var_n)),
                              ymax = var_mean + (var_sd/sqrt(var_n)),
                              x = weeks,
                              color = Treatment_group), width=0.1) +
            stat_pvalue_manual(var_lm, y.position = positiony, label = "p.signif",
                               remove.bracket = TRUE, bracket.size = 0) +
            scale_color_jama() + 
            scale_y_continuous(limits = c(min_var,max_var)) +
            theme_Publication() +
            labs(x = "Weeks", y = yname, title = title, color = ""))
}

save_function_mono <- function(plot, name, a = 5, b = 4){
    ggsave(plot = plot, 
           filename = str_c("results/monocytes/", name, ".pdf"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/monocytes/", name, ".svg"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/monocytes/", name, ".png"), width = a, height = b)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS") %>% select(ID, Sex, Age, Treatment_group)
compliance <- readRDS("data/compliance_incl_pharmacy.RDS") %>% select(1:7, 31, 34)
df_cov <- right_join(df, compliance, by = c("ID"))
mono <- rio::import("data/monocytes_n_tidy.xlsx") %>% 
    mutate(`CD11b+` = `CD11b+ (x10^4)` * 10000)
tnfa <- rio::import("data/monocytes_tnfa_tidy.xlsx")
il6 <- rio::import("data/monocytes_il6_tidy.xlsx")
elisa <- dplyr::bind_rows(tnfa, il6)

#### Datasets for ELISA analysis ####
df_elisa <- right_join(df_cov, elisa, by = c("ID")) %>% 
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

df_elisa_wide <- df_elisa %>% pivot_wider(.,
            names_from = "measurement", values_from = "concentration") %>% 
    right_join(., mono[,1:3])

# 1st stim raw values adjusted for monocyte count, then control substr
# 2st stim fold change compared to 1st
df_elisa_new <- df_elisa_wide %>% 
    mutate(
        tnf_mono_adj = TNFa / `CD11b+ (x10^4)`,
        il6_mono_adj = IL6 / `CD11b+ (x10^4)`
    )

df_elisa_prop <- df_elisa_new %>% group_by(ID, visit) %>% 
    mutate(#il6_fc_ctrl = log2(IL6/ IL6[which(stim == "ctr")]), 
           #tnfa_fc_ctrl = log2(TNFa / TNFa[which(stim == "ctr")]),
           il6_substr_ctrl = il6_mono_adj - il6_mono_adj[which(stim == "ctr")],
           tnfa_substr_ctrl = tnf_mono_adj - tnf_mono_adj[which(stim == "ctr")]
           )

df_elisa_prop2 <- df_elisa_prop %>% group_by(ID, visit, stim) %>% 
    mutate(il6_fc_train = il6_mono_adj / il6_mono_adj[which(trained == 0)], 
           tnfa_fc_train = tnf_mono_adj / tnf_mono_adj[which(trained == 0)]
    )

df_palm <- df_elisa_prop2 %>% filter(stim == "Palmitate 100uM" & trained == 0) 
df_palm_mean <- df_palm %>%   group_by(Treatment_group, weeks) %>% 
    summarise(across(c("tnfa_substr_ctrl", "il6_substr_ctrl"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

df_lps <- df_elisa_prop2 %>% filter(stim == "LPS 1ng/ml" & trained == 0)
df_lps_mean <- df_lps %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("tnfa_substr_ctrl", "il6_substr_ctrl"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

df_palm_trained <- df_elisa_prop2 %>% filter(stim == "Palmitate 100uM" & trained == 1) 
df_palm_trained_mean <- df_palm_trained  %>% group_by(Treatment_group, weeks) %>% 
    summarise(across(c("il6_fc_train", "tnfa_fc_train"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

df_lps_trained <- df_elisa_prop2 %>% filter(stim == "LPS 1ng/ml" & trained == 1) 
df_lps_trained_mean <- df_lps_trained %>% group_by(Treatment_group, weeks) %>% 
    summarise(across(c("il6_fc_train", "tnfa_fc_train"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))


#### ELISA LMMs ####
palm_lm_tnf <- df_palm %>% linearmixed_mono_elisa(tnfa_substr_ctrl)
palm_lm_il6 <- df_palm %>% linearmixed_mono_elisa(il6_substr_ctrl)

lps_lm_tnf <- df_lps %>% linearmixed_mono_elisa(tnfa_substr_ctrl)
lps_lm_il6 <- df_lps %>% linearmixed_mono_elisa(il6_substr_ctrl)

palmt_lm_il6 <- df_palm_trained %>% linearmixed_mono_elisa(il6_fc_train)
palmt_lm_tnf <- df_palm_trained %>% linearmixed_mono_elisa(tnfa_fc_train)

lpst_lm_il6 <- df_lps_trained %>% linearmixed_mono_elisa(il6_fc_train)
lpst_lm_tnf <- df_lps_trained %>% linearmixed_mono_elisa(tnfa_fc_train)

#### ELISA plots ####
(firststim_tnf_palm <- function_plotelisa(df_palm, df_palm_mean,
                                          tnfa_substr_ctrl, tnfa_substr_ctrl_mean, tnfa_substr_ctrl_sd, 
                                          tnfa_substr_ctrl_n,
                                         "TNFa (palmitate 100uM)", "TNFa pg/ml/10.000",
                                         palm_lm_tnf, 800))

(firststim_il6_palm <- function_plotelisa(df_palm, df_palm_mean,
                                          il6_substr_ctrl, il6_substr_ctrl_mean, il6_substr_ctrl_sd, 
                                          il6_substr_ctrl_n,
                                          "IL6 (palmitate 100uM)", "IL6 pg/ml/10.000",
                                          palm_lm_il6, 800))

(firststim_tnf_lps <- function_plotelisa(df_lps, df_lps_mean,
                                          tnfa_substr_ctrl, tnfa_substr_ctrl_mean, tnfa_substr_ctrl_sd, 
                                          tnfa_substr_ctrl_n,
                                         "TNFa (LPS 1ng/ml)", "TNFa pg/ml/10.000",
                                         lps_lm_tnf, 5000))

(firststim_il6_lps <- function_plotelisa(df_lps, df_lps_mean,
                                          il6_substr_ctrl, il6_substr_ctrl_mean, il6_substr_ctrl_sd, 
                                          il6_substr_ctrl_n,
                                          "IL6 (LPS 1ng/ml)", "IL6 pg/ml/10.000",
                                          lps_lm_il6, 8000))


(trained_tnf_palm <- function_plotelisa(df_palm_trained, df_palm_trained_mean,
                   tnfa_fc_train, tnfa_fc_train_mean, tnfa_fc_train_sd, tnfa_fc_train_n,
                   "TNFa (palmitate + LPS)", "fold change TNFa",
                   palmt_lm_tnf, 4))

(trained_il6_palm <- function_plotelisa(df_palm_trained, df_palm_trained_mean,
                   il6_fc_train, il6_fc_train_mean, il6_fc_train_sd, il6_fc_train_n,
                   "IL6 (palmitate + LPS)", "fold change IL6",
                   palmt_lm_il6, 5))

(trained_tnf_lps <- function_plotelisa(df_lps_trained, df_lps_trained_mean,
                   tnfa_fc_train, tnfa_fc_train_mean, tnfa_fc_train_sd, tnfa_fc_train_n,
                   "TNFa (LPS + LPS)", "fold change TNFa",
                   lpst_lm_tnf, 1))

(trained_il6_lps <- function_plotelisa(df_lps_trained, df_lps_trained_mean,
                   il6_fc_train, il6_fc_train_mean, il6_fc_train_sd, il6_fc_train_n,
                   "IL6 (LPS + LPS)", "fold change IL6",
                   lpst_lm_il6, 1))

palm_plot <- ggarrange(firststim_il6_palm, firststim_tnf_palm,
          trained_il6_palm, trained_tnf_palm,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom")

lps_plot <- ggarrange(firststim_il6_lps, firststim_tnf_lps,
          trained_il6_lps, trained_tnf_palm,
          labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2,
          common.legend = TRUE, legend = "bottom")

firststim_plot <- ggarrange(firststim_il6_palm, firststim_tnf_palm,
                            firststim_il6_lps, firststim_tnf_lps,
                           labels = c("A", "B", "C", "D"),
                           ncol = 2, nrow = 2, 
                           common.legend = TRUE, legend = "bottom")

trained_plot <- ggarrange( trained_il6_palm, trained_tnf_palm,
                           trained_il6_lps, trained_tnf_palm,
                          labels = c("A", "B", "C", "D"),
                          ncol = 2, nrow = 2, 
                          common.legend = TRUE, legend = "bottom")

total_plot <- ggarrange(firststim_il6_palm, firststim_tnf_palm,
                        firststim_il6_lps, firststim_tnf_lps,
                        trained_il6_palm, trained_tnf_palm,
                        trained_il6_lps, trained_tnf_lps,
                        nrow = 4, ncol = 2,
                        common.legend = TRUE, legend = "bottom",
                        labels = c("A", "", "B", "",
                                   "C", "", "D", ""))

save_function_mono(palm_plot, "palm_plots", a = 7, b = 8)
save_function_mono(lps_plot, "lps_plots", a = 7, b = 8)
save_function_mono(firststim_plot, "firststim", a = 7, b = 8)
save_function_mono(trained_plot, "trained_plot", a = 7, b = 8)
save_function_mono(total_plot, "total_plot", a = 7, b = 14)




#### Monocyte count data prep ####
df_mono <- right_join(df, mono, by = c("ID")) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           nonclass = `non-classical (%CD16+CD14-)`*100,
           class = `classical (%CD16-CD14+)`*100,
           im = `intermediate (%CD16+CD14+)`*100,
           nonclass_log = log(nonclass),
           im_log = log(im),
           class_log = log(class)
    ) %>% 
    droplevels(.) 

df_means_mono <- df_mono %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c( "nonclass", "im", "class",
                        "nonclass_log", "im_log", "class_log"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))


#### Monocyte LMMs ####
nonclass_lm <- df_mono %>% linearmixed_mono(nonclass)
class_lm <- df_mono %>% linearmixed_mono(class)
im_lm <- df_mono %>% linearmixed_mono(im)

#### Monocyte plots ####
(nonclassmono <- function_plotmono(df_mono, df_means_mono,
                                  nonclass, nonclass_mean, nonclass_sd, nonclass_n,
                                   "Non-classical monocytes", "proportion (%)",
                                  nonclass_lm, 20))

(classmono <- function_plotmono(df_mono, df_means_mono,
                                class, class_mean, class_sd, class_n,
                                "Classical monocytes", "proportion (%)",
                                class_lm, 80))

(immono <- function_plotmono(df_mono, df_means_mono,
                                im, im_mean, im_sd, im_n,
                                "Intermediate monocytes", "proportion (%)",
                                im_lm, 4))

mono_plot <- ggarrange(classmono, immono, nonclassmono,
          common.legend = TRUE, labels = c("A", "B", "C"), legend = "bottom",
          nrow = 2, ncol = 2)
save_function_mono(mono_plot, "monocyte_proportions", a = 8, b = 8)
