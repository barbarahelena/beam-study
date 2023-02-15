## Adverse effects and compliance
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(lme4)
library(afex)
library(ggpubr)
library(ggsci)
library(patchwork)

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

## Open data
df_med <- rio::import("data/BEAM_randomisatie_medicatie_retour.xlsx") %>% 
    select(ID = `Subject ID`, everything())
df <- readRDS("data/intervention.RDS")
metadata <- readRDS("data/demographics_BEAM.RDS")
head(df)
names(df)
df$Total_NonCompliantDays
df2 <- left_join(df, metadata, by = "ID")
head(df2)

table(df2$Treatment_group, df2$Total_NonCompliantDays)
df2 %>% group_by(Treatment_group) %>% summarise(mean(Total_NonCompliantDays),
                                                median(Total_NonCompliantDays),
                                                quantile(Total_NonCompliantDays, 1/4),
                                                quantile(Total_NonCompliantDays, 3/4))


head(df_med)
df_med <- df_med %>% filter(!ID %in% c("BEAM_713", "BEAM_664", "BEAM_299")) %>% 
    mutate(Capsules_left = as.numeric(`Capsules left`),
            Perc_pills_taken = ((788-Capsules_left)/788) *100) %>% 
    select(ID, Capsules_left, Perc_pills_taken, Group)

sum <- df_med %>% group_by(Group) %>% summarise(mean = mean(Perc_pills_taken),
                                         sd = sd(Perc_pills_taken),
                                         median = median(Perc_pills_taken),
                                         lowestq = quantile(Perc_pills_taken, 1/4),
                                         highestq = quantile(Perc_pills_taken, 3/4))

df3 <- left_join(df2, df_med, by = "ID") %>% 
    mutate(Calculated_perc = 100-(((Total_NonCompliantDays * 26)/788)*100))

df3 %>% select(Calculated_perc, Perc_pills_taken) %>% filter(Calculated_perc != 100 |
                                                                 Perc_pills_taken != 100)
df3 %>% select(Calculated_perc, Perc_pills_taken, Group) %>% filter(Calculated_perc == 100 &
                                                                 Perc_pills_taken != 100)
df3 %>% select(Calculated_perc, Perc_pills_taken, Group) %>% filter(Calculated_perc != 100 &
                                                                        Perc_pills_taken == 100)
df3 %>% select(Calculated_perc, Perc_pills_taken, Group) %>% filter(Calculated_perc >                  Perc_pills_taken )

saveRDS(df3, "data/compliance_incl_pharmacy.RDS")

## Adverse events
df_ae <- rio::import("data/BEAM_Adverse_event_export_20221219.csv") %>% 
    select(ID = "Participant Id", everything())
head(df_ae)
df_ae <- left_join(df_ae, metadata, by = "ID")
df_ae
table(df_ae$AE_Description, df_ae$Treatment_group, df_ae$ID)
