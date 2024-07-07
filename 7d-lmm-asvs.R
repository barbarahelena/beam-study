# Linear mixed for 16S
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(phyloseq)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(lme4)
library(afex)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

# Data
## Open phyloseq object
asvtab <- readRDS("data/phyloseq_BEAM_rarefied.RDS")
tax <- readRDS("data/tax_table.RDS")
asvtab_sel <- as.data.frame(as(asvtab@otu_table, "matrix"))
asvtab_rel <- (asvtab_sel / 14000) * 100
asvtab_rel_v2 <- asvtab_rel[str_detect(rownames(asvtab_rel), "V2"),]
tkb <- apply(asvtab_rel_v2, 2, function(x) sum(x > 0.05) > (0.3*length(x)))
asvtab_rel <- asvtab_rel[,tkb]
asvtab_log <- log2(asvtab_rel+1)
asvtab_log$ID <- str_extract(rownames(asvtab_log), "[A-Z]*_[0-9]*")
asvtab_log$visit <- str_extract(rownames(asvtab_log), "V[0-9]")
asvno <- ncol(asvtab_rel)

# Metadata
clindata <- readRDS("data/demographics_BEAM.RDS")
sampleinfo <- rio::import("data/fecessamples.csv")
clindf <- left_join(sampleinfo, clindata, by = "ID")
names(clindf)
head(clindf)
clindf <- clindf %>% mutate(
    sampleID = str_c(ID, "_", visit),
    weeks = case_when(
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
    before_after = fct_relevel(before_after, "After", after = 2L)
) %>% 
    filter(ID %in% asvtab_log$ID)
clindf <- clindf %>% arrange(sampleID)

df_tot <- left_join(asvtab_log, clindf, by = c("ID", "visit"))

statres <- c()
for(i in c(1:asvno)) {
    df_tot1 <- df_tot %>% filter(visit %in% c("V2", "V4"))
    df_tot1$ASV <- df_tot1[,i]
    asvname <- colnames(df_tot1)[i]
    model1 <- lmer(ASV ~ Treatment_group*weeks + (1|ID), 
                      data = df_tot1)
    res <- summary(model1)
    confint_model1 <- confint(model1)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model1[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model1[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(asvname, group1 = 0, group2 = 4, pval, 
                          sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
    
    df_tot2 <- df_tot %>% filter(visit %in% c("V4", "V5"))
    df_tot2$ASV <- df_tot2[,i]
    asvname <- colnames(df_tot2)[i]
    model2 <- lmer(ASV ~ Treatment_group*weeks + (1|ID), 
                   data = df_tot2)
    res <- summary(model2)
    confint_model2 <- confint(model2)
    estimate <- as.numeric(format(round(res$coefficients[4,1], 3), nsmall = 3))
    conflow <- as.numeric(format(round(confint_model2[6,1], 3), nsmall = 3))
    confhigh <- as.numeric(format(round(confint_model2[6,2], 3), nsmall = 3))
    pval <- format(round(res$coefficients[4,5], 3), nsmall = 3)
    pval <- as.numeric(pval)
    sig <- case_when(
        pval <= 0.05 ~paste0("*"),
        pval < 0.01 ~paste0("**"),
        pval < 0.001 ~paste0("***"),
        pval > 0.05 ~paste0("")
    )
    statres_line <- cbind(asvname, group1 = 4, group2 = 5, pval, sig, estimate, conflow, confhigh)
    statres <- rbind(statres, statres_line)
}

statres <- as.data.frame(statres)
statres$tax <- tax$Tax[match(statres$asvname, tax$ASV)]
statres <- statres %>% arrange(pval, group1) %>% 
    mutate(qval = p.adjust(pval, method = "fdr"),
           group1 = as.numeric(group1),
           group2 = as.numeric(group2))
head(statres, n = 20)
maxsig <- statres %>% filter(pval <= 0.05) %>% filter(!duplicated(asvname))

plist <- list()
for(i in 1:nrow(maxsig)){
    taxname <- maxsig$tax[i]
    pval <- maxsig$pval[i]
    asvnumber <- maxsig$asvname[i]
    df_tot$asv <- df_tot[,maxsig$asvname[i]]
    df_means <- df_tot %>% group_by(Treatment_group, weeks) %>% 
        summarise(mean = mean(asv), sd = sd(asv), n = length(asv))
    res_lmm <- statres %>% filter(tax == taxname) %>% dplyr::select(-asvname) %>% filter(sig != "")
    asvmax <- max(df_tot$asv*1.1)
    pl2 <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = asvmax),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = df_tot, aes(x = weeks, y = asv,
                                    color = Treatment_group, group = ID), alpha = 0.2, linewidth = 0.5) +
        geom_point(data = df_tot, aes(x = weeks, y = asv,
                                     color = Treatment_group, group = Treatment_group), alpha = 0.2, size = 0.8) +
        geom_line(data = df_means, aes(x = weeks, y = mean, 
                                          color = Treatment_group, group = Treatment_group), alpha = 1, linewidth = 0.8) +
        geom_point(data = df_means, aes(x = weeks, y = mean, 
                                           color = Treatment_group, group = Treatment_group), alpha = 1, size = 1.3) +
        geom_errorbar(data = df_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(res_lmm, y.position = asvmax*0.8, label = "{sig}", 
                           tip.length = 0, bracket.shorten = 0.1, size = 5, hide.ns = TRUE) +
        scale_color_jama() + 
        coord_cartesian(ylim = c(0,asvmax)) +
        theme_Publication() +
        labs(x = "Weeks", y = "log2(counts)", title = taxname,
             color = "")
        plist[[i]] <- pl2
}

(plots <- ggarrange(plotlist = plist, common.legend = TRUE, legend = "bottom",
          labels = LETTERS[1:12],
          nrow = 4, ncol = 3))
ggsave(plots, filename = "results/16S/lmer/lmer_plots.pdf", width = 10, height = 13)
ggsave(plots, filename = "results/16S/lmer/lmer_plots.png", width = 10, height = 14)
ggsave(plots, filename = "results/16S/lmer/lmer_plots.svg", width = 10, height = 14)
write.csv2(statres, "results/16S/lmer/lmm_results.csv")
