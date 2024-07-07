# BEAM trial analyses

# Project aims
In this randomized, placebo-controlled double-blind trial, we aimed to investigate the impact of oral butyrate treatment on BP in patients with grade I hypertension. The primary outcome was daytime systolic blood pressure as measured with an ambulatory blood pressure device.

# Paper
A paper with more details on the methodology can be found here: [insert link].

# Trial design
We included Dutch participants, aged 40 to 65, with untreated grade I hypertension (office systolic BP between 140 and 160 mmHg or diastolic BP between 90 and 100 mmHg), or the use of one BP lowering drug, and BMI lower than 27 kg/m2. Included participants were randomized using a stratified alternating block randomization (1:1 ratio) performed in the data collection system Castor EDC. Block sizes were 4, 6, or 8 with two strata for age (=<50 years and >50 years) and sex, to ensure that the two groups had similar baseline risk and were balanced in numbers. Physicians and study participants were blinded to treatment allocation, until the statistical analyses were completed in a blinded fashion.

# RStudio and renv
The analyses were performed in RStudio (v.2023.9.1.494) using R (v.4.2.1). We used renv and uploaded a lockfile in this repository to reconstruct the renv.

# ToC: scripts
The scripts shared in this repository are numbered, and are discussed below in the same order:

1. Data cleaning (`1_data_cleaning.R`). 

2. Tables (`2_tables.R`): for table 1, and some descriptives for adverse events and compliance numbers. Results (in csv) can be found in results folder. 

3. Dietary data: linear mixed models in `3-dietary-analysis.R` (plot can be found in Supplements)

4. Linear mixed models of blood pressure outcomes:
- Plots and linear mixed models of office and ambulatory blood pressure in `4a-bp-analyses.R` (Figure 1)
- Adjusted linear mixed models for ambulatory blood pressure in `4b-bp-lmm-adjusted.R` (Figure 2)

5. Neurohumoral balance:
- Nexfin analyses (heart rate variability, baroreceptor sensitivity, etc) in `5a-nexfin-analyses.R`
- Body impedance analysis (BIA), including total body water, in `5b-bia-analyses.R`
- Fractional excretion of sodium in `5c-frac-excretion.R`
- Renin and aldosterone plasma concentrations in `5d-reninaldo.R`. This script also includes a ggarrange function that uses outputs of scripts 5a-5c above and results in Figure 3 in the manuscript.

6. SCFA and serotonin analyses:
- Plasma and fecal SCFA analyses, plots and adjusted LMMs in `6a-scfa-analyses` which is Figure 4 in the manuscript.
- Serotonin plasma level analysis in `6b-serotonin-analysis`. This plot can be found in the Supplements.

7. Gut microbiota composition:
- Preparation of phyloseq object in `7a-processing-phyloseq.R` - make, root and add tree after multiple sequence alignment, make nice taxonomy.
- Compositional plot at genus level in `7b-compositional-plots.R`
- Ordination and diversity plots in `7c-ordination-plots.R` - this script also includes a ggarrange function that uses the compositional plot of 7b.
- Linear mixed models for separate ASVs in `7d-lmm-asvs.R` (plot can be found in Supplements)

8. Monocyte and inflammatory markers analyses:
- Analyses of monocyte FACS data and ex vivo stimulations (ELISA data IL6 and TNFa) in script `8a-monocyte-analysis.R`
- Analyses of circulating inflammatory markers in serum (ELISA data IL6 and IFNg) in script `8b-inflammatory-markers.R`
- Analyses of T-cell FACS data in script `8c-tcell-analysis.R`

