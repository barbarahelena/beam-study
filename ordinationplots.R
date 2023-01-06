## 16S analyses
## Barbara Verhaar

## libraries
library(phyloseq)
library(tidyverse)
library(rio)
library(ggsci)
library(mixOmics)
library(ggpubr)
library(vegan)

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

braypcoa <- function(bray, clinicaldata = clindf, visitname){
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
    pcoord <- ape::pcoa(bray, correction = "cailliez")
    expl_variance <- pcoord$values$Relative_eig * 100
    x_comp <- paste0('Axis.',1)
    y_comp <- paste0('Axis.',2)
    dbray <- pcoord$vectors[, c(x_comp, y_comp)]
    dbray <- as.data.frame(dbray)
    
    # add metadata / covariates
    dbray$sampleID <- rownames(dbray)
    dbray <- left_join(dbray, clindf, by = 'sampleID')
    head(dbray)
    
    plotbray <- dbray %>% 
        ggplot(aes(Axis.1, Axis.2)) +
        scale_fill_nejm() +
        geom_point(aes(color = Treatment_group), size = 2) +
        xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
        theme_Publication() +
        scale_color_lancet() +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        stat_ellipse(aes(color = Treatment_group), type = "norm", level = 0.8) +
        ggtitle(str_c(visitname))
    return(plotbray)
}

linearmixed_div <- function(data, var){
    library(lme4)
    library(afex)
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

save_function <- function(plot, name, width = 5, height = 4){
    ggsave(plot = plot, 
           filename = str_c("results/16S/", name, ".pdf"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/16S/", name, ".svg"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/16S/", name, ".png"), width = width, height = height)
}

## Data
asvtab <- readRDS("data/phyloseq_BEAM_rarefied.RDS")
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

## Ordination plots
# calculate Weighted Unifrac
wunifrac <- UniFrac(asvtab, normalized = T, weighted = T)

# PCoA
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)

# get PCoA coordinates
dfpc <- pcoord$vectors[, c(x_comp, y_comp)]
dfpc <- as.data.frame(dfpc)

# add metadata / covariates
dfpc$sampleID <- rownames(dfpc)
dfpc <- left_join(dfpc, clindf, by = 'sampleID')

pl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = Treatment_group, shape = visit), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) 
    #stat_ellipse(aes(color = Treatment_group), type = "norm")
pl
ggsave("results/16S/PCoA_WeightedUnifrac_total.pdf", device = "pdf", width = 6, height = 5)

pl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = Treatment_group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = Treatment_group), type = "norm") +
    facet_wrap(~visit)

pl
ggsave("results/16S/PCoA_WeightedUnifrac_Visit.pdf", device = "pdf", width = 6, height = 5)

pl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = visit), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = visit), type = "norm") +
    facet_wrap(~Treatment_group) + 
    ggtitle("Beta diversity: weighted UniFrac")
#guides(shape = guide_legend(override.aes = list(size = 4)))
pl
ggsave("results/16S/PCoA_WeightedUnifrac_Treatment.pdf", device = "pdf", width = 6, height = 5)

## for Bray-Curtis PCoA
names(asvtab@sam_data)
asvtab_v2 <- subset_samples(asvtab, Time_Point == 0)
asvtab_v4 <- subset_samples(asvtab, Time_Point == 1)
asvtab_v5 <- subset_samples(asvtab, Time_Point == 2)
mat <- as(asvtab@otu_table, 'matrix')
bray <- vegan::vegdist(mat, method = 'bray')

mat2 <- as(asvtab_v2@otu_table, 'matrix')
brayv2 <- vegan::vegdist(mat2, method = 'bray')

mat4 <- as(asvtab_v4@otu_table, 'matrix')
brayv4 <- vegan::vegdist(mat4, method = 'bray')

mat5 <- as(asvtab_v5@otu_table, 'matrix')
brayv5 <- vegan::vegdist(mat5, method = 'bray')

# PCoA
pcoord <- ape::pcoa(bray, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)
dbray <- pcoord$vectors[, c(x_comp, y_comp)]
dbray <- as.data.frame(dbray)

# add metadata / covariates
dbray$sampleID <- rownames(dbray)
dbray <- left_join(dbray, clindf, by = 'sampleID')
head(dbray)

pl2 <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = Treatment_group), size = 2) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = Treatment_group), type = "norm") +
    facet_wrap(~visit)
pl2
ggsave("results/PCA_BrayCurtis_Visit.pdf", device = "pdf", width = 6, height = 5)

pl2 <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = visit), size = 2) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = visit), type = "norm") +
    facet_wrap(~Treatment_group)
pl2
ggsave("results/PCA_BrayCurtis_Visit.pdf", device = "pdf", width = 6, height = 5)

# PCoA per visit, see function at the top
braypcoa(brayv2, clindf, visitname = "visit 2")
braypcoa(brayv4, clindf, visitname = "visit 4")
braypcoa(brayv5, clindf, visitname = "visit 5")
ggsave(str_c("results/PCA_BrayCurtis_", visit, ".pdf"), 
       device = "pdf", width = 6, height = 5)

## CLR-transformed PCA
pseudocount <- min(mat[mat != 0]) / 2
pc <- mixOmics::pca(mat + pseudocount, center = T, 
                    scale = F, logratio = 'CLR', ) # scale should be F when CLR is used
expl_variance <- pc$prop_expl_var$X * 100
dp <- as.data.frame(pc$variates$X)

# add metadata
rownames(dp)
dp$Treatment_group <- clindf$Treatment_group[match(rownames(dp), clindf$sampleID)]
dp$visit <- clindf$visit[match(rownames(dp), clindf$sampleID)]
#dp$Subject <- m$Subject[match(rownames(dp), m$Sample)]

# plot PCA CLR-transformed
pl <- dp %>% 
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(color = Treatment_group), size = 2) +
    xlab(paste0('PC1 (', round(expl_variance[1], digits = 1),'%)')) + # PC, not PCo
    ylab(paste0('PC2 (', round(expl_variance[2], digits = 1),'%)')) +
    #scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    ggtitle('PCA CLR-transformed') + 
    stat_ellipse(aes(color = Treatment_group), level = 5, type = "euclid") + 
    facet_wrap(~visit)
pl
ggsave("results/PCA_CLR.pdf", device = "pdf", width = 6, height = 5)

## PERMANOVA / adonis
set.seed(1234)
# distance matrix and metadata (df with the outcome / covariates) must have the same sample order as bray distance matrix / distance object
# you can use bray distance of weighted unifrac distance
all(clindf$sampleID == sample_names(asvtab)) # FALSE
clindf2 <- clindf %>%
    slice(match(sample_names(asvtab), sampleID)) %>% 
    select(sampleID, ID, visit, Treatment_group)
all(clindf2$sampleID == sample_names(asvtab)) # TRUE
dim(clindf2)

nestdf <- clindf2 %>% group_by(visit) %>% nest()
r2 <- adonis2(brayv2 ~ Treatment_group, data = nestdf$data[[2]]) 
r4 <- adonis2(brayv4 ~ Treatment_group, data = nestdf$data[[3]]) 
r5 <- adonis2(brayv5 ~ Treatment_group, data = nestdf$data[[1]]) 

# names(clindf2)
# res <- adonis(bray ~ Age + Sex + Treatment_group, data = clindf2) 

## Bray curtis with PERMANOVA annotation
pl <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    ggtitle("Bray curtis distance PCoA") +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
pl <- pl + annotate("text", x = 0.2, y = 0.4, label = paste0("PERMANOVA p = ", res$aov.tab[3,6]))

ggsave("results/PCA_BrayCurtis_permanova.pdf", device = "pdf", width = 6, height = 5)

## Alpha diversity
# Shannon
shannon <- vegan::diversity(asvtab@otu_table, index = 'shannon')
df_shan <- data.frame(sampleID = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, clindf, by = "sampleID")

(pl4 <- ggplot(data = df_shan, aes(x = Treatment_group, y = shannon, 
                                  fill = Treatment_group)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() + 
    scale_fill_jama(guide = FALSE) + 
    labs(title = "Shannon alpha diversity", y = "Shannon index", x="") +
    stat_compare_means(method = "wilcox.test", label = "p.format", hide.ns = TRUE) +
    facet_wrap(~before_after))
save_function(pl4, "shannondiversity", width = 6, height = 5)

shan_lmm <- linearmixed_div(df_shan, shannon)
shan_means <- df_shan %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(mean = mean(shannon), sd = sd(shannon), n = length(shannon))

(plot_shan <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 3.0, ymax = 5.0),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = shan_means, aes(x = weeks, y = mean, 
                                         color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = df_shan, aes(x = weeks, y = shannon,
                                    color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = shan_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(shan_lmm, y.position = 4.7, label = "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(3.0,5.0), breaks = seq(from = 3.0, to = 5.0, by = 0.5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Shannon index", title = "Shannon diversity",
             color = ""))
save_function(plot_shan, "shannon_lmm")

# Richness (ASV / Species)
richness <- vegan::specnumber(asvtab@otu_table)
dfveg <- data.frame(sampleID = names(richness), richness = richness)
dfveg <- left_join(dfveg, clindf, by = "sampleID")

(richnesspl <- ggplot(data = dfveg, aes(x = Treatment_group, y = richness, fill = Treatment_group)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_jama(guide = FALSE) + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(~before_after))
save_function(richnesspl, "richness", width = 6, height = 5)

rich_lmm <- linearmixed_div(dfveg, richness)
rich_means <- dfveg %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(mean = mean(richness), sd = sd(richness), n = length(richness))

(plot_rich <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 100, ymax = 305),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = rich_means, aes(x = weeks, y = mean, 
                                         color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = dfveg, aes(x = weeks, y = richness,
                                    color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = rich_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(rich_lmm, y.position = 275, label = "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(100,305), breaks = seq(from = 100, to = 300, by = 25)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Richness", title = "Richness",
             color = ""))
save_function(plot_rich, "richness_lmm")

# Faith's PD
faith <- picante::pd(asvtab@otu_table, tree = asvtab@phy_tree)
dffai <- as.data.frame(faith)
dffai$sampleID <- rownames(faith)
dffai <- left_join(dffai, clindf, by = "sampleID")

(faithviolin <- ggplot(data = dffai, aes(x = Treatment_group, y = PD, fill = Treatment_group)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_jama(guide = "none") + 
    labs(title = "Faith's PD", y = "Faith's PD", x = "") +
    stat_compare_means(method = "wilcox.test", label = "p.format") +
    facet_wrap(~before_after))
save_function(faithviolin, "faiths_pd", width = 6, height = 5)

div_lmm <- linearmixed_div(dffai, PD)
faith_means <- dffai %>% 
    group_by(Treatment_group, weeks) %>% 
    summarise(mean = mean(PD), sd = sd(PD), n = length(PD))

(plot_faith <- ggplot() +
        geom_rect(aes(xmin = 0, xmax = 4, ymin = 10, ymax = 35),
                  fill = "#CDCDCD", alpha = 0.3) +
        geom_line(data = faith_means, aes(x = weeks, y = mean, 
                                color = Treatment_group, group = Treatment_group), alpha = 1) +
        geom_line(data = dffai, aes(x = weeks, y = PD,
                                color = Treatment_group, group = ID), alpha = 0.2) +
        geom_errorbar(data = faith_means,
                      aes(ymin = mean - (sd/sqrt(n)),
                          ymax = mean + (sd/sqrt(n)),
                          x = weeks,
                          color = Treatment_group), width=0.1) +
        stat_pvalue_manual(div_lmm, y.position = 30, label = "p = {pval}", 
                           remove.bracket = TRUE, bracket.size = 0) +
        scale_color_jama() + 
        scale_y_continuous(limits = c(10,35), breaks = seq(from = 10, to = 35, by = 5)) +
        theme_Publication() +
        labs(x = "Weeks", y = "Faith's PD", title = "Faith's PD",
             color = ""))

save_function(plot_faith, "faith_pd_lmm")

ggarrange(p3, pl, pl4, labels = c("A", "B", "C"))
ggsave("results/descriptives.pdf", device = "pdf", width = 12, height = 10)
ggsave("results/descriptives.svg", device = "svg", width = 12, height = 10)
