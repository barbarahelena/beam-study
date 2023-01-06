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
    mutate(Treatment_week = as.numeric(Treatment_group) * weeks)
rownames(input_metadata) <- input_metadata$sampleID

fit_data = Maaslin2(
    input_data = input_data, 
    input_metadata = input_metadata, 
    min_abundance = 0.25,
    min_prevalence = 0.3,
    output = "results/16S/maaslin2", 
    fixed_effects = c("Treatment_group", "weeks", "Treatment_week"),
    normalization = "TSS",
    analysis_method = "LM",
    transform = "LOG",
    correction = "BH",
    random_effects = c("ID"))

