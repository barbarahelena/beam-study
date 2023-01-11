## Maaslin2
## b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(ggsci)
library(ggpubr)
library(Maaslin2)

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

## Open phyloseq object
asvtab <- readRDS("data/phyloseq_BEAM_rarefied.RDS")
asvtab_v2v4 <- subset_samples(asvtab, Time_Point == 0 | Time_Point == 1)
tax <- readRDS("data/tax_table.RDS")

# Metadata
clindata <- readRDS("data/demographics_BEAM.RDS")
sampleinfo <- rio::import("data/fecessamples.csv")
clindf <- left_join(sampleinfo, clindata, by = "ID")
clindata <- readRDS("data/clinicaldata.RDS")
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
)
clindf <- clindf %>% arrange(sampleID)


# Maaslin2
input_data <- as.data.frame(as(asvtab_v2v4@otu_table, "matrix"))
input_metadata <- clindf %>% dplyr::select(sampleID, ID, Treatment_group, weeks) %>% 
    mutate(Treatment_week = as.numeric(Treatment_group) * (weeks +1))
rownames(input_metadata) <- input_metadata$sampleID

fit_data = Maaslin2(
    input_data = input_data, 
    input_metadata = input_metadata, 
    min_abundance = 0.1,
    min_prevalence = 0.3,
    output = "results/16S/maaslin2", 
    fixed_effects = c("Treatment_group", "weeks", "Treatment_week"),
    normalization = "TSS",
    analysis_method = "LM",
    transform = "LOG",
    correction = "BH",
    random_effects = c("ID"))

res <- fit_data$results
res2 <- res %>% filter(metadata == "Treatment_week") %>% 
    mutate(qval2 = p.adjust(pval, method = "fdr")) %>% 
    filter(pval < 0.05) %>% 
    mutate(taxonomy = tax$Tax[match(feature, tax$ASV)])

asvtab_df <- as.data.frame(as(asvtab@otu_table, "matrix"))
asvtab_df <- asvtab_df[,str_c(res2$feature)]
asvtab_df$sampleID <- rownames(asvtab_df)
asvtab_df <- left_join(asvtab_df, clindf)

for(i in 1:8){
    taxname <- res2$taxonomy[i]
    asvtab_df$asv <- (asvtab_df[,i] / 14000) * 100
    comp <- list(c("Before", "Treatment"), c("Treatment", "After"))
    pl <- ggplot(data = asvtab_df, aes(x=before_after, y=asv)) +
        geom_boxplot(aes(fill=Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(title = taxname, x = "", y = "Relative abundance (%)") +
        facet_wrap(~Treatment_group) +
        stat_compare_means(method = "wilcox.test", label = "p.format",
                           comparisons = comp) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(pl, filename = str_c("results/16S/maaslin2/boxplot_asv_", i, ".pdf"),
           device = "pdf", width = 4, height = 5)
    
    pval <- res2$pval[i]
    asvtab_df$asv <- log(asvtab_df[,i] + 1)
    asvtab_mean <- asvtab_df %>% group_by(Treatment_group, weeks) %>% 
        summarise(mean = mean(asv), sd = sd(asv), n = length(asv))
    
    linearmixed_maaslin <- function(data, var){
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
    
    res_lmm <- linearmixed_maaslin(asvtab_df, asv)
    ymax <- max(asvtab_df$asv)*1.1
    pl2 <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 0, ymax = ymax),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = asvtab_mean, aes(x = weeks, y = mean, 
                                color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = asvtab_df, aes(x = weeks, y = asv,
                                color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = asvtab_mean,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(res_lmm, y.position = ymax*0.67, label = "p = {pval}",
                           remove.bracket = TRUE) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(0,ymax)) +
        theme_Publication() +
        labs(x = "Weeks", y = "log(counts)", title = taxname,
             color = "")
    ggsave(pl2, filename = str_c("results/16S/maaslin2/lm_asv_", i, ".pdf"),
           device = "pdf", width = 4, height = 5)
}
