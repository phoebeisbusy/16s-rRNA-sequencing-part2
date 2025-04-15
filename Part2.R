getwd()
setwd("/Users/phoebemac/EAM_Lacto2_16S")

suppressMessages({
  library(tidyverse)
  library(phyloseq)
  library(Biostrings)
  library(gridExtra)
})

load("./EAM_Lacto2_nochim.rda")

seqfile <- "dada2_Lacto2_nochim.fa"
treefile <- "dada2_Lacto2_nochim.tree"
ls()

metadata_file <- "metadata_file.tsv"
metadata_df = read.csv(
  file = metadata_file,
  strip.white = TRUE,
  sep = '\t',
  header = TRUE,
  colClasses = c("FileFwd" = "character", "FileRev" = "character")
)

seqtab.nochim.0 <- seqtab.nochim
taxa.0 <- taxa
metadata_df.0 <- metadata_df
track.0 = track

colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa) <- paste0("ASV_", seq(1:nrow(taxa)))

rownames(metadata_df) <- metadata_df$SampleID
seqtab.nochim[1:3,1:3]
taxa[1:3,]
metadata_df[1:3,]

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(metadata_df),
               tax_table(taxa))
ps

# Add rep seqs
repseqs <- readDNAStringSet(seqfile)
# repseqs %>% head()
ps <- merge_phyloseq(ps, repseqs)
ps

# Add tree
tree <- read_tree(treefile)
ps <- merge_phyloseq(ps, tree)
ps

otu_table(ps)[1:5,1:3]
sample_data(ps) %>% head()
tax_table(ps)[1:5,]

nsamples(ps)
ntaxa(ps)
sample_names(ps) %>% head()
taxa_names(ps) %>% head(n=10)
rank_names(ps)

ps.clean <- subset_taxa(ps, Kingdom == "Bacteria" &
                          !is.na(Phylum) &
                          Class != "Chloroplast" & 
                          Family != "Mitochondria")
ps.clean

#any(sample_sums(ps.clean) < 1000)
#prune_samples(sample_sums(ps.clean) < 30000, ps.clean)
#sample_sums(ps.clean) < 30000
#any(taxa_sums(ps.clean) == 0)
#sum(taxa_sums(ps.clean) == 0)
#prune_taxa(taxa_sums(ps.clean) < 100, ps.clean)
#prune_taxa(taxa_sums(ps.clean) < 1500, ps.clean)

ps.clean.re <- transform_sample_counts(ps.clean, function(x) x / sum(x) *100 )

#####RA
library(RColorBrewer)
n <- 8
m<- 125
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert=T)]
color2<-sample(color,n)
color1 = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert=T)]

color3<-sample(color1,m)

micro_type <- c("LF", "LF_Lacto")
ps.clean.re.LF1 <- subset_samples(ps.clean.re, Microbiome %in% micro_type)
ps.clean.re.LF <- subset_samples(ps.clean.re.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))




###Family Level
relative_abundance_Family <- plot_bar(ps.clean.re.LF, x = "SampleID", fill = "Family")+
  geom_bar(aes(), stat = "identity", position = "stack", width = 0.9) + 
  scale_fill_manual(values=c(color3)) +
  theme(axis.text.x = element_blank(),    # Remove x-axis text
        axis.ticks.x = element_blank())   # Remove x-axis ticks
#+ geom_text(aes(label = Genus), size = 3,  position = position_stack(vjust = 0.5))
relative_abundance_Family

A = c("LF_Control_Day0", "LF_Control_Day7","LF_Control_Day35",
      "LF_Lacto_Control_Day0",
      "LF_Lacto_Control_Day7",
      "LF_Lacto_Control_Day35",
      "LF_EAM_Day0", "LF_EAM_PreEAM","LF_EAM_PostEAM",
      "LF_Lacto_EAM_Day0",
      "LF_Lacto_EAM_PreEAM",
      "LF_Lacto_EAM_PostEAM",
      "LFR_EAM_Day0",
      "LFR_EAM_PreEAM",
      "LFR_EAM_PostEAM")

relative_abundance_Family$data$Group <- as.character(relative_abundance_Family$data$Group)
relative_abundance_Family$data$Group <- factor(relative_abundance_Family$data$Group, levels=A)
print(relative_abundance_Family) ->RA2

RA3 = RA2 + facet_wrap(~Group, scales = "free_x", shrink = TRUE, ncol = 6)
RA3
ggsave("Gut_barplot_by_family_pub.pdf", width = 20, height = 8, RA3)

###Genus Level
topN <- 150
topN <- names(sort(taxa_sums(ps.clean.re.LF), decreasing=TRUE))[1:topN]
ps.clean.topN <- prune_taxa(topN, ps.clean.re.LF)
ps.clean.topN

relative_abundance_Family <- plot_bar(ps.clean.topN, x = "SampleID", fill = "Genus")+
  geom_bar(aes(), stat = "identity", position = "stack") + 
  scale_fill_manual(values=c(color3)) + 
  theme(axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) 

#+ geom_text(aes(label = Genus), size = 3,  position = position_stack(vjust = 0.5))
relative_abundance_Family

A = c("LF_Control_Day0", "LF_Control_Day7","LF_Control_Day35",
      "LF_Lacto_Control_Day0",
      "LF_Lacto_Control_Day7",
      "LF_Lacto_Control_Day35",
      "LF_EAM_Day0", "LF_EAM_PreEAM","LF_EAM_PostEAM",
      "LF_Lacto_EAM_Day0",
      "LF_Lacto_EAM_PreEAM",
      "LF_Lacto_EAM_PostEAM",
      "LFR_EAM_Day0",
      "LFR_EAM_PreEAM",
      "LFR_EAM_PostEAM")

relative_abundance_Family$data$Group <- as.character(relative_abundance_Family$data$Group)
relative_abundance_Family$data$Group <- factor(relative_abundance_Family$data$Group, levels=A)
print(relative_abundance_Family) ->RA2

RA3 = RA2 + facet_wrap(~Group, scales = "free_x", shrink = TRUE, ncol = 6)
RA3
ggsave("Gut_barplot_by_genus150_pub.pdf", width = 20, height = 8, RA3)





###Species Level
topN <- 150
topN <- names(sort(taxa_sums(ps.clean.re.LF), decreasing=TRUE))[1:topN]
ps.clean.topN <- prune_taxa(topN, ps.clean.re.LF)
ps.clean.topN
relative_abundance_Family <- plot_bar(ps.clean.topN, x = "SampleID", fill = "Species")+
  geom_bar(aes(), stat = "identity", position = "stack", width = 0.9) + 
  scale_fill_manual(values=c(color3)) +
  theme(axis.text.x = element_blank(),    # Remove x-axis text
        axis.ticks.x = element_blank())   # Remove x-axis ticks
#+ geom_text(aes(label = Genus), size = 3,  position = position_stack(vjust = 0.5))
relative_abundance_Family

A = c("LF_Control_Day0", "LF_Control_Day7","LF_Control_Day35",
      "LF_Lacto_Control_Day0",
      "LF_Lacto_Control_Day7",
      "LF_Lacto_Control_Day35",
      "LF_EAM_Day0", "LF_EAM_PreEAM","LF_EAM_PostEAM",
      "LF_Lacto_EAM_Day0",
      "LF_Lacto_EAM_PreEAM",
      "LF_Lacto_EAM_PostEAM",
      "LFR_EAM_Day0",
      "LFR_EAM_PreEAM",
      "LFR_EAM_PostEAM")

relative_abundance_Family$data$Group <- as.character(relative_abundance_Family$data$Group)
relative_abundance_Family$data$Group <- factor(relative_abundance_Family$data$Group, levels=A)
print(relative_abundance_Family) ->RA2

RA3 = RA2 + facet_wrap(~Group, scales = "free_x", shrink = TRUE, ncol = 6)
RA3
ggsave("Gut_barplot_by_species_pub.pdf", width = 16, height = 8, RA3)





####For EAM figures
###Family level only showing EAM mice 
#########Relative abundance
n <- 6
m<- 200
color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert=T)]
color2<-sample(color,n)
color1 = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert=T)]
color3<-sample(color1,m)

ps.clean.re.1 <- subset_samples(ps.clean.re.LF, Treatment == "EAM")
relative_abundance_Family <- plot_bar(ps.clean.re.1, x = "SampleID", fill = "Phylum")+
  geom_bar(aes(), stat = "identity", position = "stack") + 
  scale_fill_manual(values=c(color3)) +
  theme(axis.text.x = element_blank(),    # Remove x-axis text
        axis.ticks.x = element_blank())   # Remove x-axis ticks


relative_abundance_Family

A = c("LF_EAM_Day0", "LF_EAM_PreEAM","LF_EAM_PostEAM",
      "LF_Lacto_EAM_Day0","LF_Lacto_EAM_PreEAM","LF_Lacto_EAM_PostEAM")

relative_abundance_Family$data$Group <- as.character(relative_abundance_Family$data$Group)
relative_abundance_Family$data$Group <- factor(relative_abundance_Family$data$Group, levels=A)
print(relative_abundance_Family) ->RA2

RA3 = RA2 + facet_wrap(~Group, labeller = labeller(Group = 
                                                     c("LF_EAM_Day0" = "LF: Baseline",
                                                       "LF_EAM_PreEAM" = "LF: PreEAM",
                                                       "LF_EAM_PostEAM" = "LF: PostEAM",
                                                       "LF_Lacto_EAM_Day0" = "Lacto: PreLacto",
                                                       "LF_Lacto_EAM_PreEAM" = "Lacto: PostLacto & PreEAM",
                                                       "LF_Lacto_EAM_PostEAM" = "Lacto: PostEAM")),
                       scales = "free_x", shrink = TRUE, ncol = 3)
RA3
#ggsave("Gut_barplot_by_family(onlyEAM).pdf", width = 14, height = 6, RA3)



######For EAM Story Figure
getwd()
setwd("/Users/phoebemac/EAM_Lacto2_16s/Relative_Abundance")
write.csv(taxa, "taxa.csv")
#import taxa_1.csv
taxa_1_df <- as.data.frame(taxa_1)

colnames(seqtab.nochim) <- paste0("ASV_", seq(1:ncol(seqtab.nochim)))
rownames(taxa_1_df) <- paste0("ASV_", seq(1:nrow(taxa_1_df)))
taxa_1_matrix <- as.matrix(taxa_1_df)

ps_1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(metadata_df),
               tax_table(taxa_1_matrix))
ps_1

ps.clean_1 <- subset_taxa(ps_1, Kingdom == "Bacteria" &
                          !is.na(Phylum) &
                          Class != "Chloroplast" & 
                          Family != "Mitochondria")
ps.clean_1

micro_type <- c("LF", "LF_Lacto")
group_type <- c("LF_EAM_Day0", "LF_EAM_PreEAM","LF_EAM_PostEAM",
                "LF_Lacto_EAM_Day0","LF_Lacto_EAM_PreEAM","LF_Lacto_EAM_PostEAM", "LF_Control_Day35","LF_Lacto_Control_Day35")
p1 <- subset_samples(ps.clean_1, Microbiome %in% micro_type
                     & Group %in% group_type)
p1
p2 <- subset_samples(p1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))
# Assume 'ps.clean.re.1' is a phyloseq object with your data
# First, calculate the relative abundance if not already done
ps.clean.re.1.rel <- transform_sample_counts(p2, function(x) x / sum(x) * 100)  # Convert to percentage

# Summarize the relative abundance per Family
family_abundance <- psmelt(ps.clean.re.1.rel) %>%
  group_by(Family) %>%
  summarize(mean_abundance = mean(Abundance))  # Calculate mean abundance across samples

# Identify families with >1% relative abundance
high_abundance_families <- family_abundance %>% 
  filter(mean_abundance > 0.008) %>%
  pull(Family)  # Extract the names of families to be shown in the legend


# Now plot, keeping all families in the plot but filtering the legend
relative_abundance_Family <- plot_bar(ps.clean.re.1.rel, x = "SampleID", fill = "Family") +
  geom_bar(aes(), stat = "identity", position = "stack") + 
  scale_fill_manual(
    values = color3,  # Assuming color3 is a vector of colors
    breaks = high_abundance_families  # Only show families with >1% abundance in the legend
  ) +
  theme(
    axis.text.x = element_blank(),    # Remove x-axis text (sample names)
    axis.ticks.x = element_blank()    # Remove x-axis ticks
  )
relative_abundance_Family

A = c("LF_EAM_Day0", "LF_EAM_PreEAM","LF_EAM_PostEAM", "LF_Control_Day35",
      "LF_Lacto_EAM_Day0","LF_Lacto_EAM_PreEAM","LF_Lacto_EAM_PostEAM", "LF_Lacto_Control_Day35")

relative_abundance_Family$data$Group <- as.character(relative_abundance_Family$data$Group)
relative_abundance_Family$data$Group <- factor(relative_abundance_Family$data$Group, levels=A)
print(relative_abundance_Family) ->RA2

RA3 = RA2 + facet_wrap(~Group, labeller = labeller(Group = 
                                                     c("LF_EAM_Day0" = "LF: Baseline",
                                                       "LF_EAM_PreEAM" = "LF: PreEAM",
                                                       "LF_EAM_PostEAM" = "LF: PostEAM",
                                                       "LF_Control_Day35" = "LF: Control Day21",
                                                       "LF_Lacto_EAM_Day0" = "Lacto: PreLacto",
                                                       "LF_Lacto_EAM_PreEAM" = "Lacto: PostLacto & PreEAM",
                                                       "LF_Lacto_EAM_PostEAM" = "Lacto: PostEAM",
                                                       "LF_Lacto_Control_Day35" = "Lacto: Control Day21")),
                       scales = "free_x", shrink = TRUE, ncol = 4)
RA3
ggsave("Gut_barplot_by_family(onlyEAM)_1.pdf", width = 10, height = 6, RA3)







#####Alpha diversity 
library(ggpubr)
library(circlize)
library(devtools)

ps.clean.LF1 <- subset_samples(ps.clean, Microbiome %in% c("LF", "LF_Lacto"))
ps.clean.LF <- subset_samples(ps.clean.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))

a1 <- estimate_richness(ps.clean.LF, measures = c("Observed","Chao1", "Shannon"))
a1 <- cbind(a1, sample_data(ps.clean.LF))

###compare between LF versus Lacto at different timepoint
p1 <- ggplot(data = a1,aes(x = Microbiome, y = Shannon)) +
  geom_boxplot(aes(fill = Microbiome)) + 
  geom_jitter(width = 0.2,alpha= 0.3) +
  facet_wrap(~Timepoint, ncol = 4,
             labeller = labeller(Timepoint = 
                                   c("Day0" = "Baseline",
                                     "Day7" = "PostLacto",
                                     "Day35" = "Control: Day21",
                                     "PostEAM" = "EAM: PostEAM"))) +
  geom_signif(comparisons = list(c("LF","LF_Lacto")), 
              test = "wilcox.test", 
              test.args=list(alternative = "two.sided", 
                             var.equal = FALSE, paired=FALSE),
              map_signif_level=TRUE, size = 1,
              y_position = 3)+
  theme_bw()+
  scale_fill_manual(
    breaks = c("LF","LF_Lacto"),
    labels = c("LF", "LF+Lacto"),
    values = c("#FF7F50", "#FFC300")) +
  scale_x_discrete(labels = c("LF" = "LF", 
                              "LF_Lacto" = "LF+Lacto"))

Time_order1= c("Day0", "Day7", "PostEAM", "Day35")
p1$data$Timepoint <- as.character(p1$data$Timepoint)
p1$data$Timepoint <- factor(p1$data$Timepoint, levels=Time_order1)
print(p1) -> a1
ggsave("Shannon.pdf", width = 6, height = 4, a1)

##Chao1
a1 <- estimate_richness(ps.clean.LF, measures = c("Observed","Chao1", "Shannon"))
a1 <- cbind(a1, sample_data(ps.clean.LF))

p1 <- ggplot(data = a1,aes(x = Microbiome, y = Chao1)) +
  geom_boxplot(aes(fill = Microbiome)) + 
  geom_jitter(width = 0.2,alpha= 0.3) +
  facet_wrap(~Timepoint, ncol = 4,
             labeller = labeller(Timepoint = 
                                   c("Day0" = "Baseline",
                                     "Day7" = "PostLacto",
                                     "Day35" = "Control: Day21",
                                     "PostEAM" = "EAM: PostEAM"))) +
  geom_signif(comparisons = list(c("LF","LF_Lacto")), 
              test = "wilcox.test", 
              test.args=list(alternative = "two.sided", 
                             var.equal = FALSE, paired=FALSE),
              map_signif_level=TRUE, size = 1,
              y_position = 160)+
  theme_bw()+
  scale_fill_manual(
    breaks = c("LF","LF_Lacto"),
    labels = c("LF", "LF+Lacto"),
    values = c("#FF7F50", "#FFC300")) +
  scale_x_discrete(labels = c("LF" = "LF", 
                              "LF_Lacto" = "LF+Lacto"))
p1

Time_order1= c("Day0", "Day7", "PostEAM", "Day35")
p1$data$Timepoint <- as.character(p1$data$Timepoint)
p1$data$Timepoint <- factor(p1$data$Timepoint, levels=Time_order1)
print(p1) -> a1

ggsave("Chao1.pdf", width = 6, height = 4, a1)










#########beta-diversity
####Bray-curtis
#Make comparision between LF and Lacto mice
micro_type <- c("LF", "LF_Lacto")
ps.clean.re.LF1 <- subset_samples(ps.clean.re, Microbiome %in% micro_type)
ps.clean.re.LF <- subset_samples(ps.clean.re.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))

ps.clean.re.bra <- ordinate(ps.clean.re.LF, method = "PCoA", distance = "bray")


pH <- plot_ordination(ps.clean.re.LF, ps.clean.re.bra, color = "Microbiome") +
  geom_point(alpha = 0.7, size = 3) +
  theme_bw() +
  scale_color_manual(
    breaks = c("LF","LF_Lacto", "LF_R"),
    labels = c("LF", "LF+Lacto", "LF+R"),
    values = c("#FF7F50", "#FFC300", "#6495ED")) +
  facet_wrap(~Timepoint, ncol = 4,
             labeller = labeller(Timepoint = 
                                   c("Day0" = "Day-14 Baseline",
                                     "Day7" = "Day-7 PostLacto",
                                     "Day35" = "Day21 Control",
                                     "PostEAM" = "Day21 EAM"))) +
  theme_minimal() + stat_ellipse(aes(group= Microbiome),type = "norm",level = 0.95, size=0.5, color= "Dark gray") +
  labs(title = "PCoA Bray Curtis") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))+    # Center title position and size
  theme(legend.title = element_text(colour="grey20", size=10, face="plain"),
        legend.text = element_text(colour="grey20", size=10, face="plain"),
        legend.position = "right",legend.background = element_rect(fill=NA, color = NA),
        legend.key=element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 45, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) +
  xlab("PC1 (44.5%)") + ylab("PC2 (18.2%)")

pH

Time_order1= c("Day0", "Day7", "PostEAM", "Day35")
pH$data$Timepoint <- as.character(pH$data$Timepoint)
pH$data$Timepoint <- factor(pH$data$Timepoint, levels=Time_order1)
print(pH) -> H1

ggsave("beta_diversity_bray.pdf", width = 8, height = 3, H1)
getwd()






###To compare how timepoint differ 
##LF microbiome
ps.clean.re.LF1 <- subset_samples(ps.clean.re, Microbiome == "LF")
ps.clean.re.LF <- subset_samples(ps.clean.re.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))
ps.clean.re.LF2 <- subset_samples(ps.clean.re.LF, Treatment == "EAM") 


ps.clean.re.bra <- ordinate(ps.clean.re.LF2, method = "PCoA", distance = "bray")
pH <- plot_ordination(ps.clean.re.LF2, ps.clean.re.bra, color = "Microbiome", shape = "Timepoint") +
  geom_point(alpha = 0.7, size = 3) +
  theme_bw() +
  scale_color_manual(
    breaks = c("LF","LF_Lacto"),
    values = c("#FF7F50", "#FFC300")) +
  scale_shape_manual(
    breaks = c("Day0", "Day7", "PostEAM"),
    values = c(1,16,17),
    labels = c("Baseline", "PreEAM", "PostEAM")) +
  theme_minimal() + stat_ellipse(aes(group= Group),type = "norm",level = 0.95,size=0.5,color= "Dark gray") +
  labs(title = "PCoA Bray Curtis | LF") +
  theme(plot.title = element_text(hjust = 0.5, size = 10))+    # Center title position and size
  theme(legend.title = element_text(colour="grey20", size=10, face="plain"),
        legend.text = element_text(colour="grey20", size=10, face="plain"),
        legend.position = "right",legend.background = element_rect(fill=NA, color = NA),
        legend.key=element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 45, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) + 
  xlab("PC1 (56.5%)") + ylab("PC2 (15.9%)")

pH

Time_order1= c("Day0", "Day7", "PostEAM")
pH$data$Timepoint <- as.character(pH$data$Timepoint)
pH$data$Timepoint <- factor(pH$data$Timepoint, levels=Time_order1)
print(pH) -> H1
ggsave("beta_diversity_bray_LF.pdf", width = 4, height = 4, H1)



###Focuse at only EAM or Control mice
time_type <- c("Day7", "PostEAM")
ps.clean.re.LF1 <- subset_samples(ps.clean.re, Timepoint %in% time_type)
ps.clean.re.LF <- subset_samples(ps.clean.re.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))
ps.clean.re.LF3 <- subset_samples(ps.clean.re.LF, Treatment == "EAM")
ps.clean.re.LF2 <- subset_samples(ps.clean.re.LF3, Microbiome %in% c("LF", "LF_Lacto")) 

ps.clean.re.bra <- ordinate(ps.clean.re.LF2, method = "PCoA", distance = "bray")
pH <- plot_ordination(ps.clean.re.LF2, ps.clean.re.bra, color = "Microbiome", shape = "Timepoint") +
  geom_point(alpha = 0.7, size = 3) +
  theme_bw() +
  scale_color_manual(
    breaks = c("LF","LF_Lacto"),
    labels = c("LF", "LF+Lacto"),
    values = c("#FF7F50", "#FFC300")) +
  scale_shape_manual(
    breaks = c("Day7", "PostEAM"),
    values = c(16, 1),
    labels = c("PreEAM", "PostEAM")) +
  theme_minimal() + stat_ellipse(aes(group= Group),type = "norm",level = 0.95,size=0.5,color= "Dark gray") +
  labs(title = "PCoA Bray Curtis | EAM Mice") +
  facet_wrap(~Microbiome, ncol = 1,
             labeller = labeller(Microbiome = c("LF" = "LF",
                                     "LF_Lacto" = "LF+Lacto")))+
  theme(plot.title = element_text(hjust = 0.5, size = 10))+    # Center title position and size
  theme(legend.title = element_text(colour="grey20", size=10, face="plain"),
        legend.text = element_text(colour="grey20", size=10, face="plain"),
        legend.position = "right",legend.background = element_rect(fill=NA, color = NA),
        legend.key=element_blank())+
  theme(panel.grid.minor = element_blank())+
  theme_bw()+
  theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 45, hjust = .5, vjust = .5, face = "plain"),
        axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),  
        axis.title.x = element_text(color = "grey20", size = 10, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 10, angle = 90, hjust = .5, vjust = .5, face = "plain")) + 
  xlab("PC1 (38.5%)") + ylab("PC2 (27.3%)") + theme(strip.text = element_blank())

pH

ggsave("beta_diversity_bray_EAM_separate_1.pdf", width = 4, height = 4, pH)

#####PERMANOVA ANALYSIS
library(ggplot2)
library(vegan)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

getwd()
setwd("/Users/phoebemac/EAM_Lacto2_16s/PERMANOVA")

####Manage Abundance file
write.csv(seqtab.nochim, "seqtab.nochim.csv")
write.csv(metadata_df, "metadata.csv")
#Combine them in excel outside of R and preseve as "Abundance.csv"

#Determine LF or Lacto preEAM vs PostEAM
alllatee<-read.csv("Abundance.csv") 
alllatee1 <- subset(alllatee, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))
alllate <- subset(alllatee1, Microbiome == "LF_Lacto" 
                  & Treatment == "EAM" 
                  & Timepoint %in% c("Day7", "PostEAM")
                  & Study == "1")
alllate

#load ASV abundance table and in Bray
sum(is.na(alllate))
alllate$Timepoint<-factor(alllate$Group)
factor_key=TRUE
alllate2<-alllate[,-c(1,2,3,4,5,6,7,8,9,10)]

set.seed(4490)
mod1<-vegdist(alllate2, method="bray", binary=FALSE, na.rm=TRUE)
set.seed(1041)
mod2=adonis(mod1~Group, data=alllate, permutations=999)
mod2

pairwise.adonis(mod1, alllate$Group)
write.csv(pairwise.adonis(mod1, alllate$Group), "EAM_Lacto_result_bray.csv") # write it to a table







#Determine LF versus Lacto at each timepoint
alllatee<-read.csv("Abundance.csv") 
# Create a new column by combining Timepoint and Microbiome
alllatee$Timepoint_Microbiome <- paste(alllatee$Timepoint, alllatee$Microbiome, sep = "_")

alllatee1 <- subset(alllatee, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))

alllate <- subset(alllatee1, Microbiome %in% c("LF", "LF_Lacto"))
alllate

#load ASV abundance table and in Bray
sum(is.na(alllate))
alllate$Timepoint<-factor(alllate$Timepoint_Microbiome)
factor_key=TRUE
alllate2<-alllate[,-c(1,2,3,4,5,6,7,8,9,10)]

set.seed(4490)
mod1<-vegdist(alllate2, method="bray", binary=FALSE, na.rm=TRUE)
set.seed(1041)
mod2=adonis(mod1~Timepoint_Microbiome, data=alllate, permutations=999)
mod2

pairwise.adonis(mod1, alllate$Timepoint_Microbiome)
write.csv(pairwise.adonis(mod1, alllate$Timepoint_Microbiome), "Timepoint_comparison_result_bray.csv") # write it to a table










#####Log2 FoldChange
suppressMessages({
  library(tidyverse)
  library(phyloseq)
  library(Biostrings)
  library(gridExtra)
  library(DESeq2)
})

########Log2 FoldChange Analyis
set.seed(100) 
foldC <- 1.0       # log2 fold change [ default ] 
alpha <- 0.05      # FDR-corrected p-value (q-value)  [ default ] 
baseM <- 10      # minimum median count required to be included [ default ] 

getwd()
setwd("/Users/phoebemac/EAM_Lacto2_16s/Log2")

###############Compare LF versus Lacto
######at each timepoint
micro_type <- c("LF", "LF_Lacto")
ps.clean.LF1 <- subset_samples(ps.clean, Microbiome %in% micro_type 
                               & Timepoint == "Day7"
                               & Treatment == "EAM")
ps.clean.LF <- subset_samples(ps.clean.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))

ps.clean.PostLacto <- ps.clean.LF

colnames(sample_data(ps.clean.PostLacto)) 
unique(sample_data(ps.clean.PostLacto)$Microbiome) 
group_name.PostLacto <- "Microbiome" 

control.PostLacto <- "LF" 
treatment.PostLacto <- "LF_Lacto" 

dds_init_PostLacto <- phyloseq_to_deseq2(ps.clean.PostLacto, ~Microbiome) 
dds_init_PostLacto  

dds_PostLacto <- DESeq(dds_init_PostLacto) # creating models 
resultsNames(dds_PostLacto) 

res.PostLacto <- results(dds_PostLacto, contrast=c("Microbiome", "LF_Lacto", "LF")) 
res.PostLacto.0 <- res.PostLacto  # keep the original structure for later 
res.PostLacto 

res.PostLacto %>% nrow() 
ncol(res.PostLacto) 
res.PostLacto <- cbind(as(res.PostLacto, "data.frame"), as(tax_table(ps.clean)[rownames(res.PostLacto), ], "matrix")) 
res.PostLacto$abs_log2FoldChange = abs(res.PostLacto$log2FoldChange) 
res.PostLacto %>% arrange(abs_log2FoldChange) 
res.PostLacto %>% arrange(desc(abs_log2FoldChange)) 

res.PostLacto <- cbind(rownames(res.PostLacto), res.PostLacto) 
colnames(res.PostLacto)[1] <- "ID" 
res.PostLacto
write.table(as.data.frame(res.PostLacto), file=paste0('PostLacto_LF_vs_Lacto.csv'),sep = "\t", row.names = FALSE) 

res.PostLacto <- res.PostLacto.0 # recover original DESeq2 structure 

sigtab.PostLacto <- res.PostLacto[which(res.PostLacto$padj <= alpha & res.PostLacto$baseMean >= baseM & abs(res.PostLacto$log2FoldChange) >= foldC) , ] 

sigtab.PostLacto %>% nrow() 
sigtab.PostLacto <- cbind(as(sigtab.PostLacto, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab.PostLacto), ], "matrix")) 
sigtab.PostLacto %>% arrange(log2FoldChange) 
sigtab.PostLacto$abs_log2FoldChange = abs(sigtab.PostLacto$log2FoldChange) 

sigtab.PostLacto <- cbind(rownames(sigtab.PostLacto), sigtab.PostLacto); 
colnames(sigtab.PostLacto)[1] <- "ID" 

write.table(as.data.frame(sigtab.PostLacto),
            file=paste0('Sig_PostLacto_LF_vs_Lacto.tsv'),sep = "\t", row.names = FALSE)

plot <- ggplot(sigtab.PostLacto, aes(y = Genus, x = log2FoldChange, color = Family)) + 
  geom_vline(xintercept = 0, color = "gray", size = 0.5) + 
  geom_point(size=2) + 
  ggtitle(paste0("PostLacto: LF vs Lacto")) +
  theme(plot.title = element_text(hjust=0.5))
plot


###############Compare between timepoint
#####PreLacto vs PostLacto
ps.clean.LF1 <- subset_samples(ps.clean, Timepoint %in% c("Day7", "PostEAM")
                               & Microbiome == "LF_Lacto"
                               & Treatment == "EAM")
ps.clean.LF <- subset_samples(ps.clean.LF1, !(SampleID %in% c("2_10A", "2_7B", "2_30C")))

ps.clean.PostLacto <- ps.clean.LF

colnames(sample_data(ps.clean.PostLacto)) 
unique(sample_data(ps.clean.PostLacto)$Timepoint) 
group_name.PostLacto <- "Timepoint" 

control.PostLacto <- "Day7" 
treatment.PostLacto <- "PostEAM" 

dds_init_PostLacto <- phyloseq_to_deseq2(ps.clean.PostLacto, ~Timepoint) 
dds_init_PostLacto  

dds_PostLacto <- DESeq(dds_init_PostLacto) # creating models 
resultsNames(dds_PostLacto) 

res.PostLacto <- results(dds_PostLacto, contrast=c("Timepoint", "PostEAM", "Day7")) 
res.PostLacto.0 <- res.PostLacto  # keep the original structure for later 
res.PostLacto 

res.PostLacto %>% nrow() 
ncol(res.PostLacto) 
res.PostLacto <- cbind(as(res.PostLacto, "data.frame"), as(tax_table(ps.clean)[rownames(res.PostLacto), ], "matrix")) 
res.PostLacto$abs_log2FoldChange = abs(res.PostLacto$log2FoldChange) 
res.PostLacto %>% arrange(abs_log2FoldChange) 
res.PostLacto %>% arrange(desc(abs_log2FoldChange)) 

res.PostLacto <- cbind(rownames(res.PostLacto), res.PostLacto) 
colnames(res.PostLacto)[1] <- "ID" 
res.PostLacto
#write.table(as.data.frame(res.PostLacto), file=paste0('LF_PreEAM_vs_PostEAM.csv'),sep = "\t", row.names = FALSE) 

res.PostLacto <- res.PostLacto.0 # recover original DESeq2 structure 

sigtab.PostLacto <- res.PostLacto[which(res.PostLacto$padj <= alpha & res.PostLacto$baseMean >= baseM & abs(res.PostLacto$log2FoldChange) >= foldC) , ] 

sigtab.PostLacto %>% nrow() 
sigtab.PostLacto <- cbind(as(sigtab.PostLacto, "data.frame"), as(tax_table(ps.clean)[rownames(sigtab.PostLacto), ], "matrix")) 
sigtab.PostLacto %>% arrange(log2FoldChange) 
sigtab.PostLacto$abs_log2FoldChange = abs(sigtab.PostLacto$log2FoldChange) 

sigtab.PostLacto <- cbind(rownames(sigtab.PostLacto), sigtab.PostLacto); 
colnames(sigtab.PostLacto)[1] <- "ID" 

#write.table(as.data.frame(sigtab.PostLacto),
#            file=paste0('Sig_LF_PreEAM_vs_PostEAM.tsv'),sep = "\t", row.names = FALSE)

plot <- ggplot(sigtab.PostLacto, aes(y = Genus, x = log2FoldChange, color = Family)) + 
  geom_vline(xintercept = 0, color = "gray", size = 0.5) + 
  geom_point(size=2) + 
  ggtitle(paste0("LF | PreEAM vs PostEAM")) +
  theme(plot.title = element_text(hjust=0.5))
plot



###Look at the relative abundance of Oscillibacter
suppressMessages({
  library(tidyverse)
  library(phyloseq)
  library(Biostrings)
  library(gridExtra)
})

library("dplyr")
library("ggpubr")
library("grid")


ps.clean.re1 <- subset_taxa(ps.clean.re, Genus == "Oscillibacter")
ps.clean.re2 <- subset_samples(ps.clean.re1, 
                               !(SampleID %in% c("2_10A", "2_7B", "2_30C")) &
                                 Microbiome %in% c("LF","LF_Lacto") &
                                 Treatment == "EAM")

Lacto <- data.frame(id = ps.clean.re2@sam_data[["SampleID"]],
                    group = ps.clean.re2@sam_data[["Group"]],
                    timepoint = ps.clean.re2@sam_data[["Timepoint"]],
                    microbiome = ps.clean.re2@sam_data[["Microbiome"]],
                    Lactobacillus = rowSums(ps.clean.re2@otu_table))


Lacto
ggboxplot(Lacto, x = "microbiome", y = "Lactobacillus",
          ylab = "Oscillibacter", xlab = "Microbiome") +
  facet_wrap(~timepoint) + 
  geom_signif(comparisons = list(c("LF", "LF_Lacto")), 
              test = "wilcox.test", 
              test.args=list(alternative = "two.sided", paired=FALSE),
              map_signif_level=TRUE, size = 0.4,
              y_position = 20, tip_length = 0.02, vjust = 0.2) 