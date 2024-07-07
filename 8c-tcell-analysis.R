# T cell results analyses
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

linearmixed_tcell <- function(data, var){
    library(lme4)
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

function_plottcell <- function(data, data_means, 
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
                               tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
            scale_color_jama() + 
            scale_y_continuous(limits = c(min_var,max_var)) +
            theme_Publication() +
            labs(x = "Weeks", y = yname, title = title, color = ""))
}

save_function_tcell <- function(plot, name, a = 5, b = 4){
    ggsave(plot = plot, 
           filename = str_c("results/tcells/", name, ".pdf"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/tcells/", name, ".svg"), width = a, height = b)
    ggsave(plot = plot, 
           filename = str_c("results/tcells/", name, ".png"), width = a, height = b)
}

#### Data ####
df <- readRDS("data/demographics_BEAM.RDS") %>% select(ID, Sex, Age, Treatment_group)
compliance <- readRDS("data/compliance_incl_pharmacy.RDS") %>% select(1:7, 31, 34)
df_cov <- right_join(df, compliance, by = c("ID"))
tcells <- rio::import("data/240503_TcellFACS_tidy.xlsx")
chrvars <- names(tcells)[str_detect(names(tcells), "%") | str_detect(names(tcells), "Freq.")]
tcells <- tcells %>% mutate(
    ID = str_extract(Idlong, "[A-Z]+_[0-9]+"),
    visit = str_extract(Idlong, "V[0-9]")
) %>% filter(!is.na(ID)) %>% 
    mutate(across(all_of(chrvars), ~str_replace(str_remove(str_remove_all(.x, " "), "%"), ",", "."))) %>% 
    mutate(across(all_of(chrvars), ~as.numeric(.x)))

#### Datasets for ELISA analysis ####
df_tcell <- right_join(df_cov, tcells, by = c("ID")) %>% 
    # filter(!ID %in% c("BEAM_299", "BEAM_664", "BEAM_713")) %>% 
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

varnames <- names(df_tcell)[14:45]

df_tcell_mean <- df_tcell %>% group_by(Treatment_group, weeks) %>% 
    summarise(across(all_of(varnames), 
                     list(mean = ~mean(.x, na.rm = TRUE),
                          sd = ~sd(.x, na.rm = TRUE),
                          n = ~length(.x)
                     ), .names = "{.col}_{.fn}"))

## Varnames
# [1] "Singlets Count"                               "population Count"                            
# [3] "CD3 Count"                                    "CD4 Count"                                   
# [5] "Singlets Freq. of Parent"                     "Singlets FSC-A, SSC-A subset Freq. of Parent"
# [7] "%CD3+ T cells"                                "%CD8+ cells of CD3+ T cells"                 
# [9] "%CD4+ cells of CD3+ T cells"                  "%conv T cells of CD4 T cells"                
# [11] "%memory T cells of conv-CD4 T cells"          "%Th17 of memory CD4 T cells"                 
# [13] "%Th22 of memory CD4 T cells"                  "%Th2 of memory CD4 T cells"                  
# [15] "%CD294+ of Th2 cells"                         "%CD294- of Th2 cells"                        
# [17] "%CD294+ Th2  of memory CD4 T cells"           "%CD294- Th2 of memory CD4 T cells"           
# [19] "%ThG of memory CD4 T cells"                   "%Th1 of memory CD4 T cells"                  
# [21] "%Th9 of memory CD4 T cells"                   "%Tfh of conv-CD4 T cells"                    
# [23] "%Tfh1 of Tfh cells"                           "%Tfh2 of Tfh cells"                          
# [25] "%Tfh17 of Tfh cells"                          "%Treg of CD4 T cells"                        
# [27] "%Tfr of CD4 Treg"                             "%memory T cells of CD4 Treg cells"           
# [29] "%Th17-like of memory Treg"                    "%Th22-like of memory Treg"                   
# [31] "%Th1-like of memory Treg"                     "%Th2-like of memory Treg" 

#### LMMs ####
## block 1 ##
singletcounts <- df_tcell %>% linearmixed_tcell(`Singlets Count`)
popcounts <- df_tcell %>% linearmixed_tcell(`population Count`)
cd3counts <- df_tcell %>% linearmixed_tcell(`CD3 Count`)
cd4counts <- df_tcell %>% linearmixed_tcell(`CD4 Count`)
singletfreqparent <- df_tcell %>% linearmixed_tcell(`Singlets Freq. of Parent`)
singletfreqsubset <- df_tcell %>% linearmixed_tcell(`Singlets FSC-A, SSC-A subset Freq. of Parent`)
cd3perc <- df_tcell %>% linearmixed_tcell(`%CD3+ T cells`)
cd8perc <- df_tcell %>% linearmixed_tcell(`%CD8+ cells of CD3+ T cells`)
cd4perc <- df_tcell %>% linearmixed_tcell(`%CD4+ cells of CD3+ T cells`)
convT <- df_tcell %>% linearmixed_tcell(`%conv T cells of CD4 T cells`)

## block 2 ##
memT <- df_tcell %>% linearmixed_tcell(`%memory T cells of conv-CD4 T cells`)
th17 <- df_tcell %>% linearmixed_tcell(`%Th17 of memory CD4 T cells`)
th22 <- df_tcell %>% linearmixed_tcell(`%Th22 of memory CD4 T cells`)
th2 <- df_tcell %>% linearmixed_tcell(`%Th2 of memory CD4 T cells`)
cd294 <- df_tcell %>% linearmixed_tcell(`%CD294+ of Th2 cells`)
cd294min <- df_tcell %>% linearmixed_tcell(`%CD294- of Th2 cells`)
cd294mem <- df_tcell %>% linearmixed_tcell(`%CD294+ Th2  of memory CD4 T cells`)
cd294memmin <- df_tcell %>% linearmixed_tcell(`%CD294- Th2 of memory CD4 T cells`)
thGmem <- df_tcell %>%  linearmixed_tcell(`%ThG of memory CD4 T cells`)
th1 <- df_tcell %>% linearmixed_tcell(`%Th1 of memory CD4 T cells`)
th9mem <- df_tcell %>% linearmixed_tcell(`%Th9 of memory CD4 T cells`)

## block 3 ##
tfhmem <- df_tcell %>% linearmixed_tcell(`%Tfh of conv-CD4 T cells`)
tfh1 <- df_tcell %>% linearmixed_tcell(`%Tfh1 of Tfh cells`)
tfh2 <- df_tcell %>% linearmixed_tcell(`%Tfh2 of Tfh cells`)
tfh17 <- df_tcell %>% linearmixed_tcell(`%Tfh17 of Tfh cells`)
treg <- df_tcell %>% linearmixed_tcell(`%Treg of CD4 T cells`)
tfr <- df_tcell %>% linearmixed_tcell(`%Tfr of CD4 Treg`)
tmem <- df_tcell %>% linearmixed_tcell(`%memory T cells of CD4 Treg cells`)
th17memreg <- df_tcell %>% linearmixed_tcell(`%Th17-like of memory Treg`)
th22memreg <- df_tcell %>% linearmixed_tcell(`%Th22-like of memory Treg`)
th1memreg <- df_tcell %>% linearmixed_tcell(`%Th1-like of memory Treg`)
th2memreg <- df_tcell %>% linearmixed_tcell(`%Th2-like of memory Treg`)


#### plots ####
## block 1 ##
(singcountplot <- function_plottcell(df_tcell, df_tcell_mean,
                                  `Singlets Count`, `Singlets Count_mean`,
                                  `Singlets Count_sd`, 
                                  `Singlets Count_n`,
                                  "Singlets Count", "Count",
                                  singletcounts, 40000))

(popcountplot <- function_plottcell(df_tcell, df_tcell_mean,
                                 `population Count`, `population Count_mean`,
                                 `population Count_sd`, 
                                 `population Count_n`,
                                 "Population Count", "Count",
                                 popcounts, 16000))

(cd3countplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `CD3 Count`, `CD3 Count_mean`,
                                `CD3 Count_sd`, 
                                `CD3 Count_n`,
                                "CD3 Count", "Count",
                                cd3counts, 9000))

(cd4countplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `CD4 Count`, `CD4 Count_mean`,
                                `CD4 Count_sd`, 
                                `CD4 Count_n`,
                                "CD4 Count", "Count",
                                cd4counts, 6500))

(cd3percplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `%CD3+ T cells`, `%CD3+ T cells_mean`,
                                `%CD3+ T cells_sd`, 
                                `%CD3+ T cells_n`,
                                "CD3+ T cells", "Percentage (%)",
                                cd3perc, 80))

(cd4percplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `%CD4+ cells of CD3+ T cells`, `%CD4+ cells of CD3+ T cells_mean`,
                                `%CD4+ cells of CD3+ T cells_sd`, 
                                `%CD4+ cells of CD3+ T cells_n`,
                                "CD4+ cells of CD3+ T cells", "Percentage (%)",
                                cd4perc, 80))

(cd8percplot <- function_plottcell(df_tcell, df_tcell_mean,
                                    `%CD8+ cells of CD3+ T cells`, `%CD8+ cells of CD3+ T cells_mean`,
                                    `%CD8+ cells of CD3+ T cells_sd`, 
                                    `%CD8+ cells of CD3+ T cells_n`,
                                    "CD8+ cells of CD3+ T cells", "Percentage (%)",
                                    cd8perc, 55))

(convtcellsplot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%conv T cells of CD4 T cells`, `%conv T cells of CD4 T cells_mean`,
                                      `%conv T cells of CD4 T cells_sd`, 
                                      `%conv T cells of CD4 T cells_n`,
                                      "conventional T cells of CD4+ T cells", "Percentage (%)",
                                      convT, 95))

plotarr1 <- ggarrange(singcountplot, popcountplot, cd3countplot, 
                      cd4countplot, cd3percplot, cd4percplot, cd8percplot,
                      convtcellsplot,
                        nrow = 4, ncol = 2)
save_function_tcell(plotarr1, "tcells_results1", a = 7, b = 14)

## block 2 ##
(memtplot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%memory T cells of conv-CD4 T cells`, 
                                      `%memory T cells of conv-CD4 T cells_mean`,
                                      `%memory T cells of conv-CD4 T cells_sd`, 
                                      `%memory T cells of conv-CD4 T cells_n`,
                                      "memory T cells of conv-CD4 T cells", "Percentage (%)",
                                      memT, 60))

(th17memt <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%Th17 of memory CD4 T cells`, `%Th17 of memory CD4 T cells_mean`,
                                      `%Th17 of memory CD4 T cells_sd`, 
                                      `%Th17 of memory CD4 T cells_n`,
                                      "Th17 of memory CD4+ T cells", "Percentage (%)",
                                      th17, 16))

(th22memt <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%Th22 of memory CD4 T cells`, `%Th22 of memory CD4 T cells_mean`,
                                      `%Th22 of memory CD4 T cells_sd`, 
                                      `%Th22 of memory CD4 T cells_n`,
                                      "Th22 of memory CD4 T cells", "Percentage (%)",
                                      th22, 7.5))

(th2memt <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%Th2 of memory CD4 T cells`, `%Th2 of memory CD4 T cells_mean`,
                                      `%Th2 of memory CD4 T cells_sd`, 
                                      `%Th2 of memory CD4 T cells_n`,
                                      "Th2 of memory CD4 T cells", "Percentage (%)",
                                      th2, 30))

(cd294plot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%CD294+ of Th2 cells`, `%CD294+ of Th2 cells_mean`,
                                      `%CD294+ of Th2 cells_sd`, 
                                      `%CD294+ of Th2 cells_n`,
                                      "CD294+ of Th2 cells", "Percentage (%)",
                                      cd294, 30))

(cd294minplot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%CD294- of Th2 cells`, `%CD294- of Th2 cells_mean`,
                                      `%CD294- of Th2 cells_sd`, 
                                      `%CD294- of Th2 cells_n`,
                                      "CD294- of Th2 cells", "Percentage (%)",
                                      cd294min, 90))

(cd294memplot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%CD294+ Th2  of memory CD4 T cells`, 
                                      `%CD294+ Th2  of memory CD4 T cells_mean`,
                                      `%CD294+ Th2  of memory CD4 T cells_sd`, 
                                      `%CD294+ Th2  of memory CD4 T cells_n`,
                                      "CD294+ Th2  of memory CD4 T cells", "Percentage (%)",
                                      cd294mem, 6))

(cd294memminplot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%CD294- Th2 of memory CD4 T cells`, 
                                      `%CD294- Th2 of memory CD4 T cells_mean`,
                                      `%CD294- Th2 of memory CD4 T cells_sd`, 
                                      `%CD294- Th2 of memory CD4 T cells_n`,
                                      "CD294- Th2 of memory CD4 T cells", "Percentage (%)",
                                      cd294memmin, 25))

(thgmemplot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%ThG of memory CD4 T cells`, `%ThG of memory CD4 T cells_mean`,
                                      `%ThG of memory CD4 T cells_sd`, 
                                      `%ThG of memory CD4 T cells_n`,
                                      "ThG of memory CD4 T cells", "Percentage (%)",
                                      thGmem, 7))

(th1plot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%Th1 of memory CD4 T cells`, `%Th1 of memory CD4 T cells_mean`,
                                      `%Th1 of memory CD4 T cells_sd`, 
                                      `%Th1 of memory CD4 T cells_n`,
                                      "Th1 of memory CD4 T cells", "Percentage (%)",
                                      th1, 30))

(th9plot <- function_plottcell(df_tcell, df_tcell_mean,
                                      `%Th9 of memory CD4 T cells`, `%Th9 of memory CD4 T cells_mean`,
                                      `%Th9 of memory CD4 T cells_sd`, 
                                      `%Th9 of memory CD4 T cells_n`,
                                      "Th9 of memory CD4 T cells", "Percentage (%)",
                                      th9mem, 24))
plotarr2 <- ggarrange(memtplot, th17memt, th22memt, th2memt, 
                      cd294plot, cd294minplot, cd294memplot, cd294memminplot,
                      thgmemplot, th1plot, th9plot,
                      nrow = 4, ncol = 3)
save_function_tcell(plotarr2, "tcells_results2", a = 10, b = 14)


## block 3 ##
(tfhmemplot <- function_plottcell(df_tcell, df_tcell_mean,
                                 `%Tfh of conv-CD4 T cells`, `%Tfh of conv-CD4 T cells_mean`,
                                 `%Tfh of conv-CD4 T cells_sd`, 
                                 `%Tfh of conv-CD4 T cells_n`,
                                 "Tfh of conv-CD4 T cells", "Percentage %",
                                 tfhmem, 15))

(tfh1plot <- function_plottcell(df_tcell, df_tcell_mean,
                                 `%Tfh1 of Tfh cells`, `%Tfh1 of Tfh cells_mean`,
                                 `%Tfh1 of Tfh cells_sd`, 
                                 `%Tfh1 of Tfh cells_n`,
                                 "Tfh1 of Tfh cells", "Percentage %",
                                 tfh1, 40))

(tfh2plot <- function_plottcell(df_tcell, df_tcell_mean,
                                `%Tfh2 of Tfh cells`, `%Tfh2 of Tfh cells_mean`,
                                `%Tfh2 of Tfh cells_sd`, 
                                `%Tfh2 of Tfh cells_n`,
                                "Tfh2 of Tfh cells", "Percentage %",
                                tfh2, 40))

(tfh17plot <- function_plottcell(df_tcell, df_tcell_mean,
                                  `%Tfh17 of Tfh cells`, `%Tfh17 of Tfh cells_mean`,
                                  `%Tfh17 of Tfh cells_sd`, 
                                  `%Tfh17 of Tfh cells_n`,
                                  "Tfh17 of Tfh cells", "Percentage %",
                                  tfh17, 45))

(tregplot <- function_plottcell(df_tcell, df_tcell_mean,
                             `%Treg of CD4 T cells`, `%Treg of CD4 T cells_mean`,
                             `%Treg of CD4 T cells_sd`, 
                             `%Treg of CD4 T cells_n`,
                             "Treg of CD4 T cells", "Percentage %",
                             treg, 9))

(tfrplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `%Tfr of CD4 Treg`, `%Tfr of CD4 Treg_mean`,
                                `%Tfr of CD4 Treg_sd`, 
                                `%Tfr of CD4 Treg_n`,
                                "Tfr of CD4 Treg", "Percentage %",
                                tfr, 20))

(tmemplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `%memory T cells of CD4 Treg cells`, 
                                `%memory T cells of CD4 Treg cells_mean`,
                                `%memory T cells of CD4 Treg cells_sd`, 
                                `%memory T cells of CD4 Treg cells_n`,
                               "memory T cells of CD4 Treg cells", "Percentage %",
                               tmem, 80))

(t17memregplot <- function_plottcell(df_tcell, df_tcell_mean,
                                `%Th17-like of memory Treg`, 
                                `%Th17-like of memory Treg_mean`,
                                `%Th17-like of memory Treg_sd`, 
                                `%Th17-like of memory Treg_n`,
                                "Th17-like of memory Treg", "Percentage %",
                                th17memreg, 45))


(t22memregplot <- function_plottcell(df_tcell, df_tcell_mean,
                                     `%Th22-like of memory Treg`, 
                                     `%Th22-like of memory Treg_mean`,
                                     `%Th22-like of memory Treg_sd`, 
                                     `%Th22-like of memory Treg_n`,
                                     "Th22-like of memory Treg", "Percentage %",
                                     th22memreg, 40))

(t2memregplot <- function_plottcell(df_tcell, df_tcell_mean,
                                     `%Th2-like of memory Treg`, 
                                     `%Th2-like of memory Treg_mean`,
                                     `%Th2-like of memory Treg_sd`, 
                                     `%Th2-like of memory Treg_n`,
                                     "Th2-like of memory Treg", "Percentage %",
                                     th2memreg, 30))

(t1memregplot <- function_plottcell(df_tcell, df_tcell_mean,
                                    `%Th1-like of memory Treg`, 
                                    `%Th1-like of memory Treg_mean`,
                                    `%Th1-like of memory Treg_sd`, 
                                    `%Th1-like of memory Treg_n`,
                                    "Th1-like of memory Treg", "Percentage %",
                                    th1memreg, 30))

plotarr3 <- ggarrange(tfhmemplot, tfh1plot, tfh2plot, tfh17plot,
                      tregplot, tfrplot, tmemplot, t17memregplot, 
                      t22memregplot, t1memregplot, t2memregplot,
                      nrow = 4, ncol = 3)
save_function_tcell(plotarr3, "tcells_results3", a = 10, b = 14)


plotsuppl1 <- ggarrange(plot_ifng, cd3percplot, cd8percplot, cd4percplot, 
                       memtplot, th1plot, th2memt, th9plot,
                       th17memt, th22memt, thgmemplot, tfhmemplot,
                       tregplot, nrow = 5, ncol = 3, labels = LETTERS[1:13],
                       common.legend = TRUE, legend = "bottom")
save_function_tcell(plotsuppl, "tcells_suppl10", a = 10, b = 16)

plotsuppl2 <- ggarrange(cd294plot, cd294minplot, tfh1plot, tfh2plot,
                        tfh17plot, tfrplot, t1memregplot, t2memregplot,
                        t17memregplot, t22memregplot,
                       nrow = 4, ncol = 3, labels = LETTERS[1:10],
                       common.legend = TRUE, legend = "bottom")
save_function_tcell(plotsuppl2, "tcells_suppl11", a = 10, b = 14)
