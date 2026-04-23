# load libraries
library(csv)
library(tidyverse)
library(ggh4x)


# load data
taxQC <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/taxQC_info_HQ_MQ_only.csv") %>%
  mutate(Tax_Class = Class) %>%
  select(-Class)
head(taxQC)

COG_groups <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/COG_info.csv")
head(COG_groups)
#View(COG_groups)

#load data and combine
COG_test1 <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/eggnog-mapper/Y5_9_cat.annotations.csv")
COG_test2 <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/eggnog-mapper/Y5_8_cat.annotations.csv")
COG_test3 <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/eggnog-mapper/W_cat.annotations.csv")
COG_test4 <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/eggnog-mapper/H_cat.annotations.csv")
COG_test5 <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/eggnog-mapper/Y5_Co_cat.annotations.csv")
colnames <- c("gene_id", "seed_ortholog","evalue","score","eggNOG_OGs","max_annot_lvl", 
                        "COG_category","Description","Preferred_name","GOs","EC","KEGG_ko",       
                        "KEGG_Pathway","KEGG_Module","KEGG_Reaction","KEGG_rclass","BRITE","KEGG_TC",       
                        "CAZy","BiGG_Reaction","PFAMs")
colnames(COG_test1) <- colnames
colnames(COG_test2) <- colnames
colnames(COG_test3) <- colnames
colnames(COG_test4) <- colnames
colnames(COG_test5) <- colnames
COG_test <- rbind(COG_test1, COG_test2, COG_test3, COG_test4, COG_test5)
head(COG_test)
#View(COG_test)

info <- COG_test %>%
  select(gene_id, Description, COG_category) %>%
  filter(str_detect(gene_id, "##", negate = TRUE)) %>%
  separate(gene_id, into = c("Genome", "gene_identifier"), sep = ".faa_", remove = FALSE) %>%
  left_join(COG_groups) %>%
  left_join(taxQC) %>%
  mutate(COG_category = ifelse(COG_category == "-", "S", COG_category),
         Class2 = ifelse(is.na(Class2), "Function unknown", Class2),
         Class = ifelse(is.na(Class), "Poorly Characterized", Class))
head(info)
#View(info)
#info$Genome %in% taxQC$Genome

# to make the stacked bar chart
info2 <- info %>%
  group_by(Genome, Class, Class2, Description) %>%
  summarise(
    Frequency = n(),
    Tax_Label = paste(Genome, ifelse(Species=="Unclassified", paste(Genus, Species, sep = " "), Species)),
    Tax_Label = str_replace_all(Tax_Label, "_concat|_assembly|pcnt|copy|-bin", ""),
    Class2 = factor(Class2, levels=c("Carbohydrate transport and metabolism",
                              "Inorganic ion transport and metabolism",
                              "Amino acid transport and metabolism",
                              "Nucleotide transport and metabolism",
                              "Energy production and conversion",
                              "Coenzyme transport and metabolism",
                              "Lipid transport and metabolism",
                              "Secondary metabolites biosynthesis, transport and catabolism",
                              "Chromatin structure and dynamics",
                              "Transcription",
                              "Translation, ribosomal structure and biogenesis",
                              "Replication, recombination and repair",
                              "RNA processing and modification",
                              "Cell wall/membrane/envelope biogenesis",
                              "Defense mechanisms",
                              "Intracellular trafficking, secretion, and vesicular transport",
                              "Signal transduction mechanisms",
                              "Cell motility",
                              "Posttranslational modification, protein turnover, chaperones",
                              "Cell cycle control, cell division, chromosome partitioning",
                              "Function unknown"))) %>%
  distinct() %>%
  arrange(Class)
head(info2)
#View(info2)

unique(info2$Class2)

# colors <- c("maroon3", "red4", "olivedrab3", "orange1", "orange3", "gold2", "lavender", "orange",
#             "purple", "grey40", "blue", "maroon", "yellow", "purple3", "firebrick", "pink",
#             "grey40", "gold", "lightskyblue", "tan4", "purple4", "olivedrab1", "darkblue", "dodgerblue3")

colors <- c("red4","maroon","maroon3","pink",
            "lavender","purple","purple3","purple4",
            "darkblue","blue","dodgerblue","lightskyblue1",
            "tan4","orange3","orange1","gold","yellow","olivedrab1","olivedrab3","olivedrab",
            "grey40")

ggplot(info2, aes(x = Tax_Label, y = Frequency, fill = Class2))+
  geom_col(position = "fill", color = "black")+
  scale_fill_manual(values = colors)+
  theme_linedraw()+
  xlab(NULL)+
  ylab("% Frequency")+
  theme(axis.title.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1, vjust = 1, color = "black"),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 12, angle = 0, hjust = 0, vjust = 0.5, color = "black"),
        legend.position = "right",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        legend.title = element_blank(),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.background = element_rect("white"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0, color = "black"),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0, color = "black"))+
  guides(fill = guide_legend(ncol = 1, byrow = FALSE))

unique(info2$Tax_Label)

info3 <- filter(info2, Tax_Label %in% c("Y5_8_28_25_1_Sub.12 Lactiplantibacillus pentosus",  "Y5_9_9B_25_1_Sub.9 Caproicibacter Unclassified",  
                                        "Y5_8_14_25_1_Sub.20 Prevotella sp002409785",        "Y5_8_14_25_1_Sub.16 Prevotella Unclassified",      
                                        "Y5_8_14_25_1_Sub.31 Lactobacillus amylovorus",      "Y5_8_14_25_1_Sub.21 Caproicibacter Unclassified",  
                                        "Y5_8_28_25_1_Sub.1 Anaerococcus sp943914475",       "Y5_8_14_25_1_Sub.30 Prevotella sp022486815",       
                                        "Y5_8_14_25_1_Sub.12 Dialister hominis",             "H2_20WD_25_1_Sub.5 Lactobacillus amylovorus",      
                                        "W7_8_14A_25_1_Sub.13 Lactobacillus sp905214545",   
                                        "Y5_8_14_25_1_Sub.34 Pseudoramibacter Unclassified", "H1_20WD_25_1_Sub.16 Lactobacillus amylovorus",     
                                        "Y5_8_14_25_1_Sub.25 Parafannyhessea sp003862195",   "Y5_8_14_25_1_Sub.11 Xanthomonas_B Unclassified",   
                                        "Y5_8_14_25_1_Sub.29 Prevotella Unclassified",       "Y5_8_14_25_1_Sub.39 Prevotella sp002305235",       
                                        "Y5_8_28_25_1_Sub.8 JALCQC01 sp022650895"))
unique(info3$Tax_Label)



ggplot(info3, aes(x = Tax_Label, y = Frequency, fill = Class2))+
  geom_col(position = "fill", color = "black")+
  scale_fill_manual(values = colors)+
  theme_linedraw()+
  xlab(NULL)+
  ylab("% Frequency")+
  theme(axis.title.y = element_text(size = 20, color = "black"),
        axis.title.x = element_text(size = 20, color = "black"),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 12, angle = 0, hjust = 0, vjust = 0.5, color = "black"),
        #legend.position = "right",
        legend.position = c(1.3,.2),
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        legend.title = element_blank(),
        plot.margin=grid::unit(c(0.1,128,0.1,55), "mm"),
        strip.background = element_rect("white"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 0, color = "black"),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0, color = "black"))+
  guides(fill = guide_legend(ncol = 1, byrow = FALSE))

head(info3)

# to make the pie chart

# inner = Class
inner <- info %>%
  count(Class) %>%
  arrange(desc(n)) %>%
  mutate(frac = n / sum(n),
         ymax = cumsum(frac),
         ymin = lag(ymax, default = 0),
         xmin = 1, xmax = 2)

# outer = Class2 (ordered by Class so chunks are contiguous)
outer <- info %>%
  count(Class, Class2) %>%
  arrange(match(Class, inner$Class), desc(n)) %>%
  mutate(frac = n / sum(n),
         ymax = cumsum(frac),
         ymin = lag(ymax, default = 0),
         xmin = 2, xmax = 3)

# colors <- c("maroon3", "red4", "olivedrab3", "orange1", "orange3", "gold2", "lavender", "orange",
#             "purple", "grey40", "blue", "maroon", "yellow", "purple3", "firebrick", "pink",
#             "grey40", "gold", "lightskyblue", "tan4", "purple4", "olivedrab1", "darkblue", "dodgerblue3")



ggplot() +
  scale_fill_manual(values = colors)+
  geom_rect(data = inner, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Class), color = "white") +
  geom_rect(data = outer, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Class2), color = "white") +
  coord_polar(theta = "y") +
  xlim(0.5, 3.5) + theme_void()


unique(outer$Class2)
unique(outer$Class)

ggplot() +
  geom_rect(data = outer,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Class2),
            color = "white", show.legend = TRUE) +
  geom_rect(data = inner,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = Class),
            color = "white", show.legend = TRUE, inherit.aes = FALSE) +
  coord_polar(theta = "y") +
  xlim(0.5, 3.5) +
  theme_void() +
  scale_fill_manual(
    name = NULL,
    values = c(
      # inner ring colors
      "Metabolism" = "firebrick",
      "Information Storage and Processing" = "blue",
      "Cellular Processes and Signaling" = "gold2",
      "Poorly characterized" = "grey40",
      # outer ring colors
      "Carbohydrate transport and metabolism" = "red4",
      "Inorganic ion transport and metabolism" = "maroon",
      "Amino acid transport and metabolism" = "maroon3",
      "Nucleotide transport and metabolism" = "pink",
      "Energy production and conversion" = "purple",
      "Coenzyme transport and metabolism" = "lavender",
      "Lipid transport and metabolism" = "purple3",
      "Secondary metabolites biosynthesis, transport and catabolism" = "purple4",
      "Transcription" = "darkblue",
      "Translation, ribosomal structure and biogenesis" = "dodgerblue3",
      "Replication, recombination and repair" = "lightskyblue",
      "RNA processing and modification" = "tan4",
      "Cell wall/membrane/envelope biogenesis" = "orange3",
      "Defense mechanisms" = "orange1",
      "Intracellular trafficking, secretion, and vesicular transport" = "yellow",
      "Signal transduction mechanisms" = "olivedrab1",
      "Cell motility" = "white",
      "Posttranslational modification, protein turnover, chaperones" = "gold",
      "Cell cycle control, cell division, chromosome partitioning" = "olivedrab3",
      "Function unknown" = "grey40"
    ),
    # legend order: inner classes first, then outer subclasses
    breaks = c(
      "Metabolism",
      "Information Storage and Processing",
      "Cellular Processes and Signaling",
      "Poorly characterized",
      "Carbohydrate transport and metabolism",
      "Inorganic ion transport and metabolism",
      "Amino acid transport and metabolism",
      "Nucleotide transport and metabolism",
      "Energy production and conversion",
      "Coenzyme transport and metabolism",
      "Lipid transport and metabolism",
      "Secondary metabolites biosynthesis, transport and catabolism",
      "Transcription",
      "Translation, ribosomal structure and biogenesis",
      "Replication, recombination and repair",
      "RNA processing and modification",
      "Cell wall/membrane/envelope biogenesis",
      "Defense mechanisms",
      "Intracellular trafficking, secretion, and vesicular transport",
      "Signal transduction mechanisms",
      "Cell motility",
      "Posttranslational modification, protein turnover, chaperones",
      "Cell cycle control, cell division, chromosome partitioning",
      "Function unknown"
    )
  ) +
  guides(fill = guide_legend(title = "Category", ncol = 1))

# Now looking at amino acids
aminos <- COG_test %>%
  filter(COG_category == "E") %>%
  select(-seed_ortholog, -score, -eggNOG_OGs, -max_annot_lvl, -BRITE, -CAZy, -BiGG_Reaction, -PFAMs,
         -KEGG_Module, -KEGG_Reaction, -KEGG_rclass, -KEGG_TC) %>%
  filter(str_detect(gene_id, "##", negate = TRUE)) %>%
  separate(gene_id, into = c("Genome", "gene_identifier"), sep = ".faa_", remove = FALSE) %>%
  left_join(COG_groups) %>%
  left_join(taxQC) %>%
  mutate(
    Aspartate = case_when(
      str_detect(EC, "2.6.1.1") ~ "Enzyme1-1_oxaloacetate",
      str_detect(EC, "6.3.5.4|6.3.1.1|3.5.1.38|3.5.1.1") ~ "Enzyme1-1_asparagine",
      str_detect(EC, "2.1.3.2|3.5.1.7") ~ "Enzyme1-1_pyrimidines",
      str_detect(EC, "4.1.1.22") ~ "Enzyme1-7_histidine",
      str_detect(EC, "1.4.3.22") ~ "Enzyme2-7_histidine",
      str_detect(EC, "1.2.1.3") ~ "Enzyme3-7_histidine",
      str_detect(EC, "1.14.13.5") ~ "Enzyme4-7_histidine",
      str_detect(EC, "3.5.2.-") ~ "Enzyme5-7_histidine",
      str_detect(EC, "3.5.3.5") ~ "Enzyme6-7_histidine",
      str_detect(EC, "3.5.1.15|3.5.1.8") ~ "Enzyme7-7_histidine"),
    Asparagine = case_when(
      str_detect(EC, "6.3.5.4") ~ "Enzyme1-3_aspartate",
      str_detect(EC, "6.3.1.1") ~ "Enzyme2-3_aspartate",
      str_detect(EC, "3.5.1.1|3.5.1.38") ~ "Enzyme3-3_aspartate"),
    Alanine = case_when(
      str_detect(EC, "4.1.1.12") ~ "Enzyme1-1_aspartate",
      str_detect(EC, "2.6.1.2|2.6.1.12|2.6.1.44|1.4.1.1") ~ "Enzyme1-1_pyruvate"),
    Arginine = case_when(
      str_detect(EC, "6.3.4.5") ~ "Enzyme1-2_aspartate",
      str_detect(EC, "4.3.2.1") ~ "Enzyme2-2_aspartate"),
    Ornithine = case_when(
      str_detect(EC, "3.5.3.1|3.5.3.27") ~ "Enzyme1-1_arginine",
      str_detect(EC, "2.3.1.1") ~ "Enzyme1-5a_glutamate",
      str_detect(EC, "2.7.2.8") ~ "Enzyme2-5a_glutamate",
      str_detect(EC, "1.2.1.38") ~ "Enzyme3-5a_glutamate",
      str_detect(EC, "2.6.1.11") ~ "Enzyme4-5a_glutamate",
      str_detect(EC, "3.5.1.14|3.5.1.16|2.3.1.35") ~ "Enzyme5-5a_glutamate",
      str_detect(EC, "6.3.2.60") ~ "Enzyme1-5b_glutamate",
      str_detect(EC, "2.7.2.19") ~ "Enzyme2-5b_glutamate",
      str_detect(EC, "1.2.1.106") ~ "Enzyme3-5b_glutamate",
      str_detect(EC, "2.7.2.124") ~ "Enzyme4-5b_glutamate",
      str_detect(EC, "3.5.1.132") ~ "Enzyme5-5b_glutamate",
      str_detect(EC, "4.3.1.12") ~ "Enzyme1-1b_proline",
      str_detect(EC, "2.6.1.13") ~ "Enzyme1-2a_proline",
      str_detect(EC, "1.5.1.2") ~ "Enzyme2-2a_proline",),
    Proline = case_when(
      str_detect(EC, "2.6.1.13") ~ "Enzyme1-2a_ornithine",
      str_detect(EC, "1.5.1.2") ~ "Enzyme2-2a_ornithine",
      str_detect(EC, "3.4.11.5") ~ "Enzyme1-1_peptide",
      str_detect(EC, "4.3.1.12") ~ "Enzyme1-1b_ornithine"),
    Glutamate = case_when(
      str_detect(EC, "1.2.1.88") ~ "Enzyme1-1_l-1-pyrroline-5-carboxylate",
      str_detect(EC, "6.3.1.2|1.4.1.13|1.4.1.14|1.4.1.2|3.5.1.2") ~ "Enzyme1-1_ammonia",
      str_detect(EC, "3.5.1.38") ~ "Enzyme1-1_aspartate, Enzyme1-1_ammonia",
      str_detect(EC, "3.5.1.1") ~ "Enzyme1-1_aspartate",
      str_detect(EC, "2.6.1.1|2.6.1.2") ~ "Enzyme1-1_2-oxoglutarate"),
    Glutamine = case_when(
      str_detect(EC, "6.3.1.2") ~ "Enzyme1-1_glutamate, Enzyme2-2_2-oxoglutarate, Enzyme2-2_ammonia",
      str_detect(EC, "1.4.1.13|1.4.1.14") ~ "Enzyme1-2_2-oxoglutarate",
      str_detect(EC, "1.4.1.2|1.4.1.3|1.4.1.4") ~ "Enzyme2-2_ammonia"),
    Homoserine = case_when(
      str_detect(EC, "2.7.2.4") ~ "Enzyme1-3_aspartate",
      str_detect(EC, "1.2.1.11") ~ "Enzyme2-3_aspartate",
      str_detect(EC, "1.1.1.3") ~ "Enzyme3-3_aspartate"),
    Glycine = case_when(
      str_detect(EC, "2.6.1.45|2.1.2.1") ~ "Enzyme1-1_serine",
      str_detect(EC, "3.5.3.3") ~ "Enzyme1-2_creatine",
      str_detect(EC, "2.1.1.20|1.5.3.1|1.5.3.24|1.5.8.3|1.5.7.3") ~ "Enzyme2-2_creatine",
      str_detect(EC, "4.1.2.5") ~ "Enzyme1-1_threonine",
      str_detect(EC, "4.1.2.48") ~ "Enzyme1-1_threonine, Enzyme2-2b_l-2-aminoacetoacetate",
      str_detect(EC, "2.3.1.29") ~ "Enzyme1-1a_l-2-aminoacetoacetate",
      str_detect(EC, "1.1.3.8") ~ "Enzyme1-2b_l-2-aminoacetoacetate",
      str_detect(EC, "4.1.2.49") ~ "Enzyme2-2b_l-2-aminoacetoacetate",
      str_detect(EC, "2.6.1.4|2.6.1.44") ~ "Enzyme1-1_glyoxylate",),
    Tryptophan = case_when(
      str_detect(EC, "4.1.3.27") ~ "Enzyme1-5_chorismate",
      str_detect(EC, "2.4.2.18") ~ "Enzyme2-5_chorismate",
      str_detect(EC, "5.3.1.24") ~ "Enzyme3-5_chorismate",
      str_detect(EC, "4.1.1.48") ~ "Enzyme4-5_chorismate",
      str_detect(EC, "4.2.1.20") ~ "Enzyme5-5_chorismate"),
    Methionine = case_when(
      str_detect(EC, "4.4.1.1") ~ "Enzyme1-3_cysteine",
      str_detect(EC, "4.4.1.13|4.2.1.22|2.5.1.134") ~ "Enzyme2-3_cysteine, Enzyme6-7_aspartate, Enzyme3-4a_homoserine, Enzyme3-4b_homoserine, Enzyme3-4c_homoserine",
      str_detect(EC, "2.1.1.5|2.1.1.10|2.1.1.13|2.1.1.14") ~ "Enzyme3-3_cysteine, Enzyme1-1_homocysteine, Enzyme7-7_aspartate, Enzyme4-4a_homoserine, Enzyme4-4b_homoserine, Enzyme4-4c_homoserine, Enzyme3-3d_homoserine, Enzyme3-3e_homoserine, Enzyme3-3f_homoserine",
      str_detect(EC, "2.7.2.4") ~ "Enzyme1-7_aspartate",
      str_detect(EC, "1.2.1.11") ~ "Enzyme2-7_aspartate",
      str_detect(EC, "1.1.1.3") ~ "Enzyme3-7_aspartate, Enzyme1-3f_homoserine",
      str_detect(EC, "2.3.1.46") ~ "Enzyme4-7_aspartate, Enzyme1-4a_homoserine, Enzyme1-3d_homoserine",
      str_detect(EC, "2.5.1.48") ~ "Enzyme5-7_aspartate, Enzyme2-4a_homoserine, Enzyme2-4c_homoserine, Enzyme2-3d_homoserine, Enzyme2-3e_homoserine",
      str_detect(EC, "2.7.1.39") ~ "Enzyme1-4b_homoserine",
      str_detect(EC, "2.5.1.160") ~ "Enzyme2-4b_homoserine",
      str_detect(EC, "2.3.1.31") ~ "Enzyme1-4c_homoserine, Enzyme1-3e_homoserine",
      str_detect(EC, "2.5.1.-") ~ "Enzyme2-3d_homoserine",
      str_detect(EC, "2.5.1.49") ~ "Enzyme2-3e_homoserine",
      str_detect(EC, "2.8.1.16") ~ "Enzyme2-3f_homoserine"),
    Phenylalanine = case_when(
      str_detect(EC, "4.2.1.20") ~ "Enzyme1-8_tryptophan",
      str_detect(EC, "4.1.1.48") ~ "Enzyme2-8_tryptophan",
      str_detect(EC, "5.3.1.24") ~ "Enzyme3-8_tryptophan",
      str_detect(EC, "2.4.2.18") ~ "Enzyme4-8_tryptophan",
      str_detect(EC, "4.1.3.27") ~ "Enzyme5-8_tryptophan",
      str_detect(EC, "5.4.99.5") ~ "Enzyme6-8_tryptophan, Enzyme1-3a_chorismate, Enzyme1-3b_chorismate",
      str_detect(EC, "4.2.1.91|4.2.1.51") ~ "Enzyme7-8_tryptophan, Enzyme2-3a_chorismate, Enzyme3-3b_chorismate",
      str_detect(EC, "2.6.1.1|2.6.1.5|2.6.1.9|2.6.1.58|1.4.1.20") ~ "Enzyme8-8_tryptophan",
      str_detect(EC, "2.6.1.57") ~ "Enzyme3-3a_chorismate, Enzyme8-8_tryptophan",
      str_detect(EC, "2.6.1.78|2.6.1.79") ~ "Enzyme2-3b_chorismate"),
    Homoserine = case_when(
      str_detect(EC, "2.7.2.4") ~ "Enzyme1-1_serine",
      str_detect(EC, "1.2.1.11") ~ "Enzyme2-3_aspartate",
      str_detect(EC, "1.1.1.3") ~ "Enzyme3-3_aspartate"),
    Threonine = case_when(
      str_detect(EC, "4.1.2.48|4.1.2.5") ~ "Enzyme1-1_glycine",
      str_detect(EC, "2.7.1.39") ~ "Enzyme1-2_homoserine",
      str_detect(EC, "4.2.3.1") ~ "Enzyme2-2_homoserine",),
    Cysteine = case_when(
      str_detect(EC, "4.4.1.13") ~ "Enzyme1-1a_pyruvate",
      str_detect(EC, "4.2.1.22") ~ "Enzyme1-1a_serine, Enzyme1-2b_serine",
      str_detect(EC, "4.4.1.1") ~ "Enzyme2-2b_serine, Enzyme1-1b_pyruvate"),
    Serine = case_when(
      str_detect(EC, "2.6.1.45|2.6.1.51") ~ "Enzyme1-1_hydroxy-pyruvate",
      str_detect(EC, "2.6.1.45|2.1.2.1") ~ "Enzyme1-1_glycine",
      str_detect(EC, "1.1.1.95") ~ "Enzyme1-3_3p-d-glycerate",
      str_detect(EC, "2.6.1.52") ~ "Enzyme2-3_3p-d-glycerate",
      str_detect(EC, "3.1.3.3") ~ "Enzyme3-3_3p-d-glycerate"),
    Lysine = case_when(
      str_detect(EC, "2.7.2.4") ~ "Enzyme1-6a_aspartate, Enzyme1-9b_aspartate, Enzyme1-9c_aspartate",
      str_detect(EC, "1.2.1.11") ~ "Enzyme2-6a_aspartate, Enzyme2-9b_aspartate, Enzyme2-9c_aspartate",
      str_detect(EC, "4.3.3.7") ~ "Enzyme3-6a_aspartate. Enzyme3-9b_aspartate, Enzyme3-9c_aspartate",
      str_detect(EC, "1.17.1.8") ~ "Enzyme4-6a_aspartate, Enzyme4-9b_aspartate, Enzyme4-9c_aspartate",
      str_detect(EC, "1.4.1.16|2.6.1.83") ~ "Enzyme5-6a_aspartate",
      str_detect(EC, "4.1.1.20") ~ "Enzyme6-6a_aspartate, Enzyme9-9b_aspartate, Enzyme9-9c_aspartate",
      str_detect(EC, "2.3.1.89") ~ "Enzyme5-9b_aspartate",
      str_detect(EC, "2.6.1.-") ~ "Enzyme6-9b_aspartate",
      str_detect(EC, "3.5.1.47") ~ "Enzyme7-9b_aspartate",
      str_detect(EC, "5.1.1.7") ~ "Enzyme8-9b_aspartate, Enzyme8-9c_aspartate",
      str_detect(EC, "2.3.1.117") ~ "Enzyme5-9c_aspartate",
      str_detect(EC, "2.6.1.17") ~ "Enzyme6-9c_aspartate",
      str_detect(EC, "3.5.1.18") ~ "Enzyme7-9c_aspartate",
      str_detect(EC, "2.3.3.14") ~ "Enzyme1-10a_2-oxogluterate, Enzyme1-10b_2-oxogluterate",
      str_detect(EC, "4.2.1.-|4.2.1.114") ~ "Enzyme2-10a_2-oxogluterate, Enzyme2-10b_2-oxogluterate",
      str_detect(EC, "4.2.1.114|4.2.1.36") ~ "Enzyme3-10a_2-oxogluterate, Enzyme3-10b_2-oxogluterate",
      str_detect(EC, "1.1.87|1.1.1.286") ~ "Enzyme4-10a_2-oxogluterate, Enzyme4-10b_2-oxogluterate",
      str_detect(EC, "2.6.1.39|2.6.1.57") ~ "Enzyme5-10a_2-oxogluterate, Enzyme5-10b_2-oxogluterate",
      str_detect(EC, "1.2.1.95") ~ "Enzyme6-10a_2-oxogluterate, Enzyme7-10a_2-oxogluterate, Enzyme8-10a_2-oxogluterate",
      str_detect(EC, "1.5.1.10") ~ "Enzyme9-10a_2-oxogluterate",
      str_detect(EC, "1.5.1.7") ~ "Enzyme10-10a_2-oxogluterate",
      str_detect(EC, "6.3.2.43") ~ "Enzyme6-10b_2-oxogluterate",
      str_detect(EC, "2.7.2.17") ~ "Enzyme7-10b_2-oxogluterate",
      str_detect(EC, "1.2.1.103") ~ "Enzyme8-10b_2-oxogluterate",
      str_detect(EC, "2.6.1.118") ~ "Enzyme9-10b_2-oxogluterate",
      str_detect(EC, "3.5.1.103") ~ "Enzyme10-10b_2-oxogluterate"),
    Isoleucine = case_when(
      str_detect(EC, "2.3.3.21") ~ "Enzyme1-9_pyruvate",
      str_detect(EC, "4.2.1.35") ~ "Enzyme2-9_pyruvate",
      str_detect(EC, "4.2.1.35") ~ "Enzyme3-9_pyruvate",
      str_detect(EC, "1.1.1.85") ~ "Enzyme4-9_pyruvate",
      str_detect(EC, "2.2.1.6") ~ "Enzyme5-9_pyruvate",
      str_detect(EC, "1.1.1.86|5.4.99.3") ~ "Enzyme6-9_pyruvate",
      str_detect(EC, "1.1.1.86") ~ "Enzyme7-9_pyruvate",
      str_detect(EC, "4.2.1.9") ~ "Enzyme8-9_pyruvate",
      str_detect(EC, "2.6.1.42|1.4.1.9") ~ "Enzyme9-9_pyruvate",
      str_detect(EC, "4.3.1.19") ~ "Enzyme1-6_threonine",
      str_detect(EC, "2.2.1.6") ~ "Enzyme2-6_threonine",
      str_detect(EC, "1.1.1.86|5.4.99.3") ~ "Enzyme3-6_threonine",
      str_detect(EC, "1.1.1.86") ~ "Enzyme4-6_threonine",
      str_detect(EC, "4.2.1.9") ~ "Enzyme5-6_threonine",
      str_detect(EC, "2.6.1.42|1.4.1.9") ~ "Enzyme6-6_threonine"),
    Valine = case_when(
      str_detect(EC, "2.2.1.6") ~ "Enzyme1-5_pyruvate",
      str_detect(EC, "1.1.1.86|5.4.99.3") ~ "Enzyme2-5_pyruvate",
      str_detect(EC, "1.1.1.86") ~ "Enzyme3-5_pyruvate",
      str_detect(EC, "4.2.1.9") ~ "Enzyme4-5_pyruvate",
      str_detect(EC, "2.6.1.42|1.4.1.9|2.6.1.66") ~ "Enzyme5-5_pyruvate"),
    Leucine = case_when(
      str_detect(EC, "2.6.1.6") ~ "Enzyme1-6_valine",
      str_detect(EC, "2.3.3.13") ~ "Enzyme2-6_valine",
      str_detect(EC, "4.2.1.33") ~ "Enzyme3-6_valine",
      str_detect(EC, "4.2.1.33") ~ "Enzyme4-6_valine",
      str_detect(EC, "1.1.1.85") ~ "Enzyme5-6_valine",
      str_detect(EC, "2.6.1.6|2.6.1.42|1.4.1.9") ~ "Enzyme6-6_valine",
      str_detect(EC, "4.2.1.33") ~ "Enzyme1-4_pyruvate",
      str_detect(EC, "4.2.1.33") ~ "Enzyme2-4_pyruvate",
      str_detect(EC, "1.1.1.85") ~ "Enzyme3-4_pyruvate",
      str_detect(EC, "2.6.1.6|2.6.1.42|1.4.1.9") ~ "Enzyme4-4_pyruvate"),
    Histidine = case_when(
      str_detect(EC, "2.6.1.38") ~ "Enzyme1-1_imidazole-pyruvate",
      str_detect(EC, "2.4.7.17") ~ "Enzyme1-10_PRPP",
      str_detect(EC, "3.6.1.31") ~ "Enzyme2-10_PRPP",
      str_detect(EC, "3.5.4.19") ~ "Enzyme3-10_PRPP",
      str_detect(EC, "5.3.1.16") ~ "Enzyme4-10_PRPP",
      str_detect(EC, "4.3.2.10") ~ "Enzyme5-10_PRPP",
      str_detect(EC, "4.2.1.19") ~ "Enzyme6-10_PRPP",
      str_detect(EC, "2.6.1.9") ~ "Enzyme7-10_PRPP",
      str_detect(EC, "3.1.3.15") ~ "Enzyme9-10_PRPP",
      str_detect(EC, "1.1.1.23") ~ "Enzyme9-10_PRPP",
      str_detect(EC, "1.1.1.23") ~ "Enzyme10-10_PRPP"),
    Tyrosine = case_when(
      str_detect(EC, "1.14.16.1") ~ "Enzyme1-1_phenylalanine",
      str_detect(EC, "5.4.99.5") ~ "Enzyme1-3a_chorismate",
      str_detect(EC, "1.3.1.12|1.3.1.13") ~ "Enzyme2-3a_chorismate",
      str_detect(EC, "2.6.1.57") ~ "Enzyme3-3a_chorismate",
      str_detect(EC, "5.4.99.5") ~ "Enzyme1-3b_chorismate",
      str_detect(EC, "2.6.1.78|2.6.1.57") ~ "Enzyme2-3b_chorismate",
      str_detect(EC, "1.3.1.43|1.3.1.78") ~ "Enzyme3-3b_chorismate"),
    Chorismate = case_when(
      str_detect(EC, "2.5.1.54") ~ "Enzyme1-7_shikimate_pthwy",
      str_detect(EC, "4.2.3.4") ~ "Enzyme2-7_shikimate_pthwy",
      str_detect(EC, "4.2.1.10") ~ "Enzyme3-7_shikimate_pthwy",
      str_detect(EC, "1.1.1.25|1.1.1.24|1.1.5.8|1.1.1.282") ~ "Enzyme4-7_shikimate_pthwy",
      str_detect(EC, "2.7.1.71") ~ "Enzyme5-7_shikimate_pthwy",
      str_detect(EC, "2.5.1.19") ~ "Enzyme6-7_shikimate_pthwy",
      str_detect(EC, "4.2.3.5") ~ "Enzyme7-7_shikimate_pthwy")
  ) %>%
  filter(!is.na(Aspartate)|!is.na(Asparagine)|!is.na(Alanine)|!is.na(Arginine)|!is.na(Ornithine)|!is.na(Proline)|
           !is.na(Glutamate)|!is.na(Glutamine)|!is.na(Glycine)|!is.na(Tryptophan)|!is.na(Phenylalanine)|!is.na(Homoserine)|
           !is.na(Threonine)|!is.na(Cysteine)|!is.na(Serine)|!is.na(Lysine)|!is.na(Isoleucine)|!is.na(Valine)|
           !is.na(Leucine)|!is.na(Histidine)|!is.na(Tyrosine)|is.na(Methionine)|!is.na(Homoserine)) %>%
  pivot_longer(cols = c("Aspartate", "Asparagine", "Alanine", "Arginine", "Ornithine",
                        "Proline" , "Glutamate", "Glutamine", "Glycine", "Tryptophan", "Phenylalanine",
                        "Homoserine", "Threonine", "Cysteine", "Serine", "Lysine", "Isoleucine",
                        "Valine", "Leucine", "Histidine", "Tyrosine", "Methionine", "Homoserine"), 
               names_to = "AminoAcid", values_to = "Pathway_Code", values_drop_na = TRUE) %>%
  separate(Pathway_Code, into = c("Pathway_Info", "Syn_From"), sep = "_") %>%
  separate(Pathway_Info, into = c("Enzyme_Pos", "Number"), sep = "-") %>%
  group_by(Genome, AminoAcid, Syn_From, Number) %>%
  summarise(
    Tax_Label = paste(Genome, ifelse(Species=="Unclassified", paste(Genus, Species, sep = " "), Species)),
    Tax_Label = str_replace_all(Tax_Label, "_concat|_assembly|pcnt|copy|-bin", ""),
    Alt_Tax_Label = paste(ifelse(Species=="Unclassified", paste(Genus, Species, sep = " "), Species), Genome),
    Alt_Tax_Label = str_replace_all(Alt_Tax_Label, "_concat|_assembly|pcnt|copy|-bin", ""),
    Num_in_pathway = as.numeric(str_extract_all(Number, "\\d+")),
    Num_present = n_distinct(Enzyme_Pos),
    Percent_complete = Num_present/Num_in_pathway,
    Pathway_identifier = paste(AminoAcid, "from", Syn_From, ifelse(is.na(str_extract(Number, "[A-Za-z]+")), "", str_extract(Number, "[A-Za-z]+")), sep = " "),
    Syn_From2 = paste(Syn_From, ifelse(is.na(str_extract(Number, "[A-Za-z]+")), "", str_extract(Number, "[A-Za-z]+")), sep = " "),
    ECs = paste(unique(EC), collapse = ", ")) %>%
  distinct() %>%
  filter(Tax_Label %in% c("Y5_8_28_25_1_Sub.12 Lactiplantibacillus pentosus",  "Y5_9_9B_25_1_Sub.9 Caproicibacter Unclassified",  
                          "Y5_8_14_25_1_Sub.20 Prevotella sp002409785",        "Y5_8_14_25_1_Sub.16 Prevotella Unclassified",      
                          "Y5_8_14_25_1_Sub.31 Lactobacillus amylovorus",      "Y5_8_14_25_1_Sub.21 Caproicibacter Unclassified",  
                          "Y5_8_28_25_1_Sub.1 Anaerococcus sp943914475",       "Y5_8_14_25_1_Sub.30 Prevotella sp022486815",       
                          "Y5_8_14_25_1_Sub.12 Dialister hominis",             "H2_20WD_25_1_Sub.5 Lactobacillus amylovorus",      
                          "W7_8_14A_25_1_Sub.13 Lactobacillus sp905214545",   
                          "Y5_8_14_25_1_Sub.34 Pseudoramibacter Unclassified", "H1_20WD_25_1_Sub.16 Lactobacillus amylovorus",     
                          "Y5_8_14_25_1_Sub.25 Parafannyhessea sp003862195",   "Y5_8_14_25_1_Sub.11 Xanthomonas_B Unclassified",   
                          "Y5_8_14_25_1_Sub.29 Prevotella Unclassified",       "Y5_8_14_25_1_Sub.39 Prevotella sp002305235",       
                          "Y5_8_28_25_1_Sub.8 JALCQC01 sp022650895"))
head(aminos)
#View(aminos)
#Enzyme1-1_peptide


ggplot(aminos, aes(x = Syn_From2, y = Alt_Tax_Label))+
  facet_grid(cols = vars(AminoAcid), scales = "free", space = "free")+
  geom_tile(aes(fill = Percent_complete), color = "black")+
  #geom_text(aes(label = Type), size = 6)+
  theme_linedraw()+
  ylab(NULL)+
  scale_fill_gradient(low = "purple4", high = "orange")+
  #scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 1, hjust = 1, color = "black"),
        axis.title.y = element_text(size = 20, color = "black"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.position = c(-0.3,-0.25),
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        legend.title = element_blank(),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x = element_text(size = 18, face = "bold", angle = 90, color = "black"),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0, color = "black"))+
  guides(fill = guide_colorbar(barheight = unit(1, "cm"), barwidth = unit(10, "cm"), direction = "horizontal"))

