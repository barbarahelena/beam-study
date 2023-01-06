## 16S analyses
## Barbara Verhaar

## libraries
library(phyloseq)
library(dplyr)

## Data
asvtab <- readRDS("data/ps.2022_28_BEAM.global.2022-12-22.curated.RDS")

sum(asvtab@otu_table[1,]) 
sum(asvtab@otu_table[2,]) # check that samples are not rarefied
asvtab@otu_table[1,1]
sample_names(asvtab)
asvtab

# make multiple sequence alignment with MAFFT binary available in the 'bin' folder
# READ: https://mafft.cbrc.jp/alignment/software/
# install mafft for mac here: https://mafft.cbrc.jp/alignment/software/macstandard.html
# version 7.5.1.1 (mafft-7.511-signed.pkg)
make_multiple_alignment <- "mafft --auto --thread 2 data/ps.2022_28_BEAM.global.2022-12-22.curated.ASV.fasta > data/ASVs.msa"
system(make_multiple_alignment)

# make phylogenetic tree with FastTreeDbl (Double Precision) --> needs FastTreeDbl available in the 'bin' folder
# READ: http://www.microbesonline.org/fasttree/#BranchLen (read for MAC OS)
# check if FasttreeDbl is available
# FastTreeDbl doesnt work on mac, therefore iqtree
# download iqtree and install
iqtree <- '/Users/barbaraverhaar/Documents/VUmc/BEAM/beam-study/data/iqtree-2.2.0-MacOSX/bin/iqtree2'
make_tree <- paste0(iqtree, " -s data/ASVs.msa -m MFP")
make_tree
system(make_tree)

# import and root tree
tree <- ape::read.tree('data/ASVs.msa.treefile')
tree_rooted <- phytools::midpoint.root(tree)

asvtab@phy_tree <- tree_rooted

sample_sums(asvtab)
asvtab
saveRDS(asvtab, 'results/phyloseq_BEAM.RDS')

## Inspect phyloseq object
ntaxa(asvtab)
nsamples(asvtab)
sample_names(asvtab)

## Rarefaction
asv_rare <- rarefy_even_depth(asvtab, sample.size = 14000, 
                              rngseed = 4321, replace = F, 
                              trimOTUs = T, verbose = T)

# Remove constant/empty ASVs
asvtab2 <- prune_taxa(taxa_sums(asv_rare) > 0, asv_rare)
nsamples(asvtab2)
ntaxa(asvtab2)
sample_names(asvtab2)

### fix ASV taxonomy
tax <- as.data.frame(as(asvtab2@tax_table, 'matrix'))
head(tax)
sum(!is.na(tax$Species)) / nrow(tax) * 100 # only 9.15 % of all ASVs have species level
sum(!is.na(tax$Genus)) / nrow(tax) * 100 # 82.49 % of all ASVs have genus level
sum(!is.na(tax$Family)) / nrow(tax) * 100 # 98.21 % of all ASVs have family level
sum(!is.na(tax$Phylum)) / nrow(tax) * 100 # 99.96 % of all ASVs have phylum level
nrow(tax)

# get top 300 ASV by abundance 
ss <- taxa_sums(asvtab2)
ss <- ss[order(ss, decreasing = T)]
ss <- ss[1:300]
top300 <- names(ss)
tax300 <- tax[rownames(tax) %in% top300, ]

sum(!is.na(tax300$Species)) / nrow(tax300) * 100 # 24.0 % of top 300 ASVs have species level
sum(!is.na(tax300$Genus)) / nrow(tax300) * 100 # 87.7 % of top 300 ASVs  have genus level
sum(!is.na(tax300$Family)) / nrow(tax300) * 100 # 100 % of top 300 ASVs have family level
sum(!is.na(tax300$Phylum)) / nrow(tax300) * 100 # 100 % of top 300 ASVs have phylum level
nsamples(asvtab2) # 62

# get 'nice' taxonomy for ASVs
tax <- tax %>% 
    mutate(Tax = case_when(
        !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species),
        !is.na(Genus) & is.na(Species) ~ paste(Genus, 'spp.'),
        !is.na(Family) & is.na(Genus) ~ paste(Family, 'spp.'),
        !is.na(Order) & is.na(Family) ~ paste(Order, 'spp.'),
        !is.na(Class) & is.na(Order) ~ paste(Class, 'spp.'),
        !is.na(Phylum) & is.na(Class) ~ paste(Phylum, 'spp.'),
        !is.na(Kingdom) & is.na(Phylum) ~ paste(Kingdom, 'spp.'),
        is.na(Kingdom) ~ 'unclassified'),
        ASV = rownames(.)
    )
unique(tax$Tax)    

saveRDS(asvtab2, file = "data/phyloseq_BEAM_rarefied.RDS")
saveRDS(tax, file = "data/tax_table.RDS")
