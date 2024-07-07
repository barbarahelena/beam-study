## Composition plot at genus level
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(phyloseq)
library(tidyverse)
library(rio)
library(ggsci)
library(mixOmics)
library(ggpubr)
library(vegan)
library(breakerofchains) # put mouse cursor on the %>%  and press Ctrl + Shift + B 

## Theme
theme_composition <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
}

save_function <- function(plot, name, width = 5, height = 4){
    ggsave(plot = plot, 
           filename = str_c("results/16S/", name, ".pdf"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/16S/", name, ".svg"), width = width, height = height)
    ggsave(plot = plot, 
           filename = str_c("results/16S/", name, ".png"), width = width, height = height)
}

cols <- c("darkgreen", 'firebrick', "navy", "dodgerblue",  "goldenrod2", "chartreuse4", "darkorange2", "rosybrown1", "darkred", "lightskyblue",
          "seagreen", "gold1", "olivedrab", "royalblue", "linen", "maroon4", "mediumturquoise", "plum2", "darkslateblue", "sienna", "grey70", "grey90")


## Data
asvtab <- readRDS("data/phyloseq_BEAM_rarefied.RDS")
tab <- as.data.frame(as(asvtab@otu_table, 'matrix'))
dim(tab)
(rarefaction_level <- sample_sums(asvtab)[1]) # rarefied to 14000
tax <- readRDS("data/tax_table.RDS")
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

# convert to rel. abundance %
tab <- (tab / rarefaction_level) * 100
rowSums(tab) # samples should all sum up to 100%


#### Genus-level ####
N <- 20
# convert to long format
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample) %>% 
    left_join(., clindf)
# add species taxonomy (including ambiguous)
d$Genus <- tax$Genus[match(d$ASV, tax$ASV)]
# remove amiguous species
d$Genus[str_detect(d$Genus, 'spp.')] <- 'ambiguous'

top_taxa <- d %>%
    group_by(Genus) %>%
    summarise(Abund = sum(Abundance)) %>%
    arrange(-Abund) %>%
    dplyr::select(Genus) %>%
    filter(Genus != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_taxa

d %>% group_by(Sample) %>% summarise(x = sum(Abundance))

# summarize abundance per group for this tax level
dx <- d %>% 
    group_by(Treatment_group, before_after, ASV) %>% 
    summarise(Abundance = mean(Abundance))
dx$Genus <- d$Genus[match(dx$ASV, d$ASV)] # add curated taxonomy

dx <- dx %>% 
    group_by(Genus, Treatment_group, before_after) %>% 
    summarise(Abundance = sum(Abundance))
dx

# check
dx %>% group_by(Treatment_group, before_after) %>% summarise(x = sum(Abundance))

dx$Genus2 <- case_when(
    dx$Genus %in% top_taxa ~ paste(dx$Genus),
    is.na(dx$Genus) ~ paste("Unknown"),
    !(dx$Genus %in% top_taxa) ~ paste("Other genera")
)

dx <- dx %>% mutate(
    Genus2 = as.factor(Genus2),
    Genus2 = fct_relevel(Genus2, "Other genera", after = 0L),
    Genus2 = fct_relevel(Genus2, "Unknown", after = 0L),
    # Genus2 = fct_recode(Genus2, `Oscillospiraceae UCG-002` = "UCG-002",
    #                     `Oscillospiraceae NK4A214 group` = "NK4A214 group")
)

lev <- levels(dx$Genus2)
lev
dx <- dx %>% group_by(Treatment_group, before_after, Genus2) %>% summarise(Abundance2 = sum(Abundance))

# library(ggtext)
comp_genus <- dx %>% 
    #filter(Genus %in% top_taxa) %>% 
    #mutate(Tax = factor(Genus, levels = rev(make.unique(top_taxa)))) %>% 
    ggplot(aes(x = before_after, y = Abundance2, fill = Genus2)) +
    geom_bar(stat = "identity", color = 'black') +
    scale_fill_manual(values = rev(cols), labels = lev) +
    guides(fill = guide_legend(title = "Genus", ncol = 1)) +
    labs(y="Composition (%)", x = "", title = "Composition (genus level)") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition() +
    theme(
        legend.text = element_markdown(),
        axis.text.x =  element_text(angle = 45, hjust = 1)
        )+
    facet_wrap(~Treatment_group)
comp_genus
save_function(comp_genus, "composition_genus", width = 7, height = 7)

