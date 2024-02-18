# Renin / aldosterone analysis
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

linearmixed_reninaldo <- function(data, var){
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

save_function_reninaldo <- function(plot, name, a = 5, b = 4){
    ggsave(plot = plot, 
           filename = str_c("results/reninaldo/", name, ".pdf"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/reninaldo/", name, ".svg"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/reninaldo/", name, ".png"), width = a, height = b)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS")
reninaldo <- readRDS("data/reninaldo.RDS")
lab <- readRDS("data/lab_results.RDS") %>% 
     select(ID, visit, Na, K, GFR, LDL, CRP, UrineSodium, contains("DW"), contains("WW"))
# urinedata <- readRDS("data/urinesamples.RDS") %>% select(ID, visit, Volume)
# dietarydata <- readRDS("data/diet_summary.RDS") %>% select(ID, visit, Sodium, Fibers)
covariates <- right_join(urinedata, df, by = "ID") %>% right_join(lab, ., by = c("ID", "visit"))
abpm <- readRDS("data/abpm_total.RDS")

#### Dataset for renin and aldosterone analysis ####
df_ra <- right_join(reninaldo, covariates, by = c("ID", "visit")) %>% 
    # filter(Renin < 15) %>% 
    filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
    mutate(Treatment_group = as.factor(Treatment_group),
           visit = as.factor(visit), 
           weeks = case_when(
               visit == "V2" ~ paste0(0),
               visit == "V4" ~ paste0(4),
               visit == "V5" ~ paste0(5)
           ),
           weeks = as.numeric(weeks),
           logRenin = log10(Renin),
           logAldo = log10(Aldosterone),
           ARR = Aldosterone / Renin,
           logARR = log10(ARR),
           boxcoxRenin = car::boxCoxVariable(Renin),
           boxcoxAldo = car::boxCoxVariable(Aldosterone),
           boxcoxARR = car::boxCoxVariable(ARR)) %>% 
    droplevels(.) 

df_means <- df_ra %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(across(c("Renin", "Aldosterone", "ARR", "logRenin", "logARR", "logAldo",
                       "boxcoxRenin", "boxcoxAldo", "boxcoxARR"), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

#### Renin and aldosterone plots with LMMs ####
renin_lm <- df_ra %>% linearmixed_reninaldo(logRenin)
aldo_lm <- df_ra %>% linearmixed_reninaldo(Aldosterone)
arr_lm <- df_ra %>% linearmixed_reninaldo(logARR)

(plot_renin <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = 22),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_ra, aes(x = weeks, y = Renin,
                                          color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_ra, aes(x = weeks, y = Renin,
                                           color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Renin_mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Renin_mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Renin_mean - (Renin_sd/sqrt(Renin_n)),
                          ymax = Renin_mean + (Renin_sd/sqrt(Renin_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(renin_lm, y.position = 1.0, label = "p.signif",
                           remove.bracket = TRUE, bracket.size = 0, size = 5) +
        scale_color_jama(guide = "none") +
        scale_y_continuous(limits = c(0,22), breaks = seq(from = 0, to = 20, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Renin", title = "Renin"))

(plot_logrenin <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0.25, ymax = 1.50),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_ra, aes(x = weeks, y = logRenin,
                                    color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_ra, aes(x = weeks, y = logRenin,
                                     color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = logRenin_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = logRenin_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = logRenin_mean - (logRenin_sd/sqrt(logRenin_n)),
                          ymax = logRenin_mean + (logRenin_sd/sqrt(logRenin_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(renin_lm, y.position = 1.25, label = "p.signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0.25,1.50), breaks = seq(from = 0.25, to = 1.50, by = 0.25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Log10 renin (pg/ml)", title = "Renin", color = ""))

(plot_logarr <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0.50, ymax = 2.0),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_ra, aes(x = weeks, y = logARR,
                                    color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_ra, aes(x = weeks, y = logARR,
                                     color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = logARR_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = logARR_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = logARR_mean - (logARR_sd/sqrt(logARR_n)),
                          ymax = logARR_mean + (logARR_sd/sqrt(logARR_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(arr_lm, y.position = 1.5, label = "p.signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0.5,2.0), breaks = seq(from = 0.50, to = 2.0, by = 0.25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Log10 ARR", title = "Aldosterone-renin ratio", color = ""))

(plot_aldo <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 45, ymax = 320),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_ra, aes(x = weeks, y = Aldosterone,
                                    color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_ra, aes(x = weeks, y = Aldosterone,
                                     color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = Aldosterone_mean, 
                                       color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = Aldosterone_mean, 
                                        color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = Aldosterone_mean - (Aldosterone_sd/sqrt(Aldosterone_n)),
                          ymax = Aldosterone_mean + (Aldosterone_sd/sqrt(Aldosterone_n)),
                          x = weeks,
                          color = Treatment_group), width=0.1, linewidth = 0.75) +
        stat_pvalue_manual(aldo_lm, y.position = 250, label = "p.signif",
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(45,320), breaks = seq(from = 50, to = 300, by = 50)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Aldosterone (pg/ml)", title = "Aldosterone", color = ""))


(pl_reninaldo <- ggarrange(plot_logrenin, plot_aldo,
                     labels = c("A", "B"),
                     nrow = 1, ncol = 2,
                     common.legend = TRUE,
                     legend = "bottom"))
save_function_reninaldo(pl_reninaldo, "reninaldo_plots", a = 9, b = 3)

save_function_reninaldo(plot_logrenin, "logrenin")
save_function_reninaldo(plot_aldo, "aldo")
save_function_reninaldo(plot_logarr, "logarr")

pl_mechanisms <- ggarrange(plot_brs, plot_sdnn, plot_dpdt,
                           plot_logrenin, plot_aldo,
                           labels = c("A", "B", "C", "D", "E"), 
                           nrow = 3,
                           ncol = 2,
                           common.legend = TRUE,
                           legend = "bottom")
pl_mechanisms
save_function_reninaldo(pl_mechanisms, "nexfin_aldorenin_fena_bia", a = 8, b = 10)

### Boxplots renin aldo ##
reninaldo_boxplots <- df_ra %>% 
    mutate(before_after = case_when(
        visit == "V2" ~ paste0("Before"),
        visit == "V4" ~ paste0("Treatment"),
        visit == "V5" ~ paste0("After")
    ),
    before_after = fct_relevel(before_after, "Before", after = 0L),
    before_after = fct_relevel(before_after, "After", after = 2L))

(boxplot_aldo <- ggplot(reninaldo_boxplots, aes(x = before_after, y = Aldosterone,
                                          fill = Treatment_group)) +
        geom_boxplot(outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        facet_wrap(~Treatment_group) +
        stat_compare_means(method = "t.test", label = "p.format",
                           comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
        scale_fill_jama() +
        scale_color_jama(guide = "none") +
        labs(x = "", title = "Aldosterone", y = "Levels") +
        theme_Publication())
save_function_bp(boxplot_sbp, "abpm", "boxplot_sbp", height = 5)

(boxplot_dbp <- ggplot(abpm_boxplots, aes(x = before_after, y = Awake_diastolic_Mean)) +
        geom_boxplot(aes(fill = Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        stat_compare_means(method = "t.test", label = "p.value",
                           comparisons = list(c("Before", "Treatment"), c("Treatment", "After"))) +
        facet_wrap(~Treatment_group) + 
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(x = "", title = "Daytime diastolic BP", y = "Diastolic BP (mmHg)") +
        theme_Publication() )
save_function_bp(boxplot_officesbp, "abpm", "boxplot_dbp", height = 5)

