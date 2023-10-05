## Differences microbiota between groups; descriptive
## Barbara Verhaar

## Libraries
library(phyloseq)
library(tidyverse)
library(ggsci)
library(ggpubr)

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

(rarefaction_level <- sample_sums(asvtab)[1]) # rarefied to 14000
phy <- tax_glom(asvtab, taxrank="Phylum")
fam <- tax_glom(asvtab, taxrank="Family")
gen <- tax_glom(asvtab, taxrank="Genus")

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
)
clindf <- clindf %>% arrange(sampleID)

## Family level
fam2 <- filter_taxa(fam, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/16S/family_boxplots")

for(i in 1:nrow(fam2@tax_table)){
    fa_name <- fam2@tax_table[,"Family"][[i]]
    famsub <- subset_taxa(fam2, Family == fa_name)
    df <- as.data.frame(as(famsub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- rownames(df)
    dfa <- left_join(df, clindf, by="sampleID")
    dfa$family <- dfa[,1]
    comp <- list(c("Before", "Treatment"), c("Treatment", "After"))
    pl <- ggplot(data = dfa, aes(x=before_after, y=family)) +
        geom_boxplot(aes(fill=Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(title = fa_name, x = "", y = "Relative abundance (%)") +
        facet_wrap(~Treatment_group) +
        stat_compare_means(method = "wilcox.test", label = "p.format",
            comparisons = comp) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(pl, filename = str_c("results/16S/family_boxplots/fam_", str_to_lower(fa_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}


## Phylum level
phy2 = filter_taxa(phy, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/16S/phylum_boxplots")

for(i in 1:nrow(phy2@tax_table)){
    ph_name <- phy2@tax_table[,"Phylum"][[i]]
    physub <- subset_taxa(phy2, Phylum == ph_name)
    df <- as.data.frame(as(physub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- rownames(df)
    dfa <- left_join(df, clindf, by="sampleID")
    dfa$phylum <- dfa[,1]
    comp <- list(c("Before", "Treatment"), c("Treatment", "After"))
    pl <- ggplot(data = dfa, aes(x=before_after, y=phylum)) +
        geom_boxplot(aes(fill=Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(title = ph_name, x = "", y = "Relative abundance (%)") +
        facet_wrap(~Treatment_group) +
        stat_compare_means(method = "wilcox.test", label = "p.format",
                           comparisons = comp) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(pl, filename = str_c("results/16S/phylum_boxplots/phy_", str_to_lower(ph_name), ".pdf"), 
           device = "pdf", width = 4, height = 5)
}

## Genus level
gen2 = filter_taxa(gen, function(x) sum(x > 3) > (0.3*length(x)), TRUE)
dir.create("results/16S/genus_boxplots")

for(i in 1:nrow(gen2@tax_table)){
    gen_name <- gen2@tax_table[,"Genus"][[i]]
    gensub <- subset_taxa(gen2, Genus == gen_name)
    df <- as.data.frame(as(gensub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- rownames(df)
    dfa <- left_join(df, clindf, by="sampleID")
    dfa$genus <- dfa[,1]
    comp <- list(c("Before", "Treatment"), c("Treatment", "After"))
    pl <- ggplot(data = dfa, aes(x=before_after, y=genus)) +
        geom_boxplot(aes(fill=Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        geom_line(aes(color = Treatment_group, group = ID), alpha = 0.2) +
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(title = gen_name, x = "", y = "Relative abundance (%)") +
        facet_wrap(~Treatment_group) +
        stat_compare_means(method = "wilcox.test", label = "p.format",
                           comparisons = comp) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(pl, filename = str_c("results/16S/genus_boxplots/gen_", str_to_lower(gen_name), ".pdf"), 
           device = "pdf", width = 5, height = 5)
}

for(i in 1:nrow(gen2@tax_table)){
    gen_name <- gen2@tax_table[,"Genus"][[i]]
    gensub <- subset_taxa(gen2, Genus == gen_name)
    df <- as.data.frame(as(gensub@otu_table@.Data, "matrix"))
    df <- (df / rarefaction_level) * 100
    df$sampleID <- rownames(df)
    dfa <- left_join(df, clindf, by="sampleID")
    dfa$genus <- dfa[,1]
    comp <- list(c("Before", "Treatment"), c("Treatment", "After"))
    pl <- ggplot(data = dfa, aes(x=Treatment_group, y=genus)) +
        geom_boxplot(aes(fill=Treatment_group), outlier.shape = NA) +
        geom_point(size = 0.75)+
        scale_fill_jama(guide = "none") +
        scale_color_jama(guide = "none") +
        labs(title = gen_name, x = "", y = "Relative abundance (%)") +
        facet_wrap(~before_after) +
        stat_compare_means(method = "wilcox.test", label = "p.format",
                           comparisons = comp) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(pl, filename = str_c("results/16S/genus_boxplots/gen_visit_", str_to_lower(gen_name), ".pdf"), 
           device = "pdf", width = 5, height = 5)
}

