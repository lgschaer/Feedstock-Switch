# load libraries
library(csv)
library(tidyverse)
library(ggh4x)
library(forcats)
#install.packages("tidytext")
library(tidytext)
library(vegan)
library(ggpubr)


# load data
taxQC <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/taxQC_info_HQ_MQ_only.csv") %>%
  mutate(Tax_Class = Class) %>%
  select(-Class)
head(taxQC)

COG_groups <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/COG_info.csv")
head(COG_groups)
#View(COG_groups)

clean_mag_ids <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/clean_mag_ids.csv")
head(clean_mag_ids)

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
#write_csv(COG_test, "/Users/lauraschaerer/Desktop/feedstock_switch/data/merged_and_cleaned_eggnog_output_11192025.csv")


clean_genes <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/metaErg_clean_genes_out_09192025.csv") %>%
  separate(Contig_id, into = c("Temp_Contig", NA), sep = "_k", remove = FALSE) %>%
  unite(gene_id, Temp_Contig, Gene_id, sep = ".faa_")
head(clean_genes)
# dim(clean_genes)
# head(unique(clean_genes$gene_id))
# head(unique(COG_test$gene_id))
head(unique(clean_genes$MAG))
head(unique(taxQC$MAG))
# anyDuplicated(clean_genes$gene_id)
# anyDuplicated(COG_test$gene_id)
# unique(clean_genes$gene_id) %in% unique(COG_test$gene_id)

info <- COG_test %>%
  select(gene_id, Description, COG_category, GOs, EC) %>%
  filter(str_detect(gene_id, "##", negate = TRUE)) %>%
  separate(gene_id, into = c("Genome", "gene_identifier"), sep = ".faa_", remove = FALSE) %>%
  left_join(COG_groups) %>%
  left_join(taxQC) %>%
  left_join(clean_genes, by = "gene_id") %>%
  left_join(clean_mag_ids) %>%
  filter(!is.na(Contig_id)) %>%
  filter(EC.x != "-") %>%
  mutate(EC = EC.x,
         COG_category = ifelse(COG_category == "-", "S", COG_category),
         Class2 = ifelse(is.na(Class2), "Function unknown", Class2),
         Class = ifelse(is.na(Class), "Poorly Characterized", Class)) %>%
  #filter(Class2 == "Carbohydrate transport and metabolism") %>%
  group_by(Genome, Class, Class2, Description, GOs, Enzyme_description, EC, Alt_Tax_Label, Quality, Completeness, Contamination) %>%
  summarise(
    Frequency = n(),
    EC = str_replace_all(EC, ",", "x"), 
    EC = paste0(EC, "x")) %>%
  distinct() #%>%
#  filter(Alt_Tax_Label %in% c("Lactobacillus amylovorus MAG.6","Lactobacillus amylovorus MAG.17","Lactobacillus sp905214545 MAG.29","Xanthomonas_B Unclassified MAG.57",
 #                             "Dialister hominis MAG.58","Prevotella Unclassified MAG.60","Prevotella sp002409785 MAG.63","Caproicibacter Unclassified MAG.64",
  #                            "Parafannyhessea sp003862195 MAG.65","Prevotella Unclassified MAG.66","Prevotella sp022486815 MAG.67","Lactobacillus amylovorus MAG.68",
   #                           "Pseudoramibacter Unclassified MAG.70","Prevotella sp002305235 MAG.73","Anaerococcus sp943914475 MAG.81","Lactiplantibacillus pentosus MAG.84",
    #                          "JALCQC01 sp022650895 MAG.86","Caproicibacter Unclassified MAG.93")) 
head(info)
#View(info)
unique(info$Class2)

Carbohydrates <- c("Hemicellulose", "Cellulose", "Galacturonan", "Chito_oligosaccharides", 
                   "Peptidoglycan", "Fructans", "Glycogen", "Polysaccharides")

carbz <- info %>%
  mutate(
    Lactate = case_when(
      str_detect(EC, "1\\.1\\.1\\.(27x|28x|436x)|1\\.1\\.2\\.(3x|4x5x)|1\\.1\\.5\\.12x|1\\.1\\.99\\.(6x|7x|40x)|1\\.1\\.3\\.2x") ~ "From pyruvate",
      str_detect(EC, "1\\.13\\.12\\.4x") ~ "From acetate",
      str_detect(EC, "3\\.7\\.1\\.26x") ~ "From L-fucose intermediate",
      str_detect(EC, "3\\.1\\.6\\.17x") ~ "From sulfolacetate",
      str_detect(EC, "4\\.1\\.1\\.101x") ~ "From malate",
      str_detect(EC, "3\\.5\\.1\\.124x") ~ "From protein",
      str_detect(EC, "4\\.2\\.1\\.126x") ~ "From acetylmuramic acid",
      str_detect(EC, "4\\.2\\.1\\.130x") ~ "From methylglyoxal",
      str_detect(EC, "2\\.8\\.3\\.1x") ~ "From propanoate, lactoyl-coa",
      str_detect(EC, "3\\.1\\.2\\.6x") ~ "From lactoylglutathione",
      str_detect(EC, "1\\.2\\.1\\.(22x|23x)") ~ "From lactaldehyde",
      str_detect(EC, "4\\.1\\.2\\.36x") ~ "From formate, acetaldehyde"),
    Ethanol = case_when(
      str_detect(EC, "3\\.1\\.1\\.67x") ~ "From fatty acyl ethyl ester",
      str_detect(EC, "3\\.5\\.1\\.75x") ~ "From urethane",
      str_detect(EC, "2\\.3\\.1\\.268x|3\\.1\\.1\\.113x") ~ "From ethyl acetate",
      str_detect(EC, "1\\.1\\.2\\.(8x|7x)|1\\.1\\.5\\.5x|1\\.1\\.99\\.36x") ~ "From acetaldehyde",
      str_detect(EC, "2\\.3\\.1\\.152x") ~ "From glucose, ethyl cinnamate",
      str_detect(EC, "1\\.1\\.1\\.(1x|2x|17x)") ~ "alcohol dehydrogenase"),
    Acetate = case_when(
      str_detect(EC, "3\\.5\\.1\\.51x") ~ "From 4-acetamindobutyril-coa",
      str_detect(EC, "3\\.5\\.1\\.56x") ~ "From 4-methylumbelliferyl-acetate",
      str_detect(EC, "3\\.5\\.1\\.63x") ~ "From 5-acetamindopentanoate",
      str_detect(EC, "3\\.1\\.1\\.33x") ~ "From 6-acetyl-d-glucose",
      str_detect(EC, "1\\.2\\.99\\.6x|1\\.2\\.1\\.(3x|4x|5x)|1\\.2\\.5\\.2x") ~ "From acetaldehyde",
      str_detect(EC, "3\\.1\\.1\\.6x") ~ "From acetic ester",
      str_detect(EC, "3\\.1\\.1\\.(113x|112x)") ~ "From X acetate",
      str_detect(EC, "3\\.5\\.1\\.4x") ~ "From acetamide",
      str_detect(EC, "6\\.2\\.1\\.(1x|2x)|3\\.6\\.1\\.20x") ~ "From acetyl adenylate",
      str_detect(EC, "3\\.1\\.1\\.7x") ~ "From acetylcholine",
      str_detect(EC, "2\\.7\\.2\\.(1x|12x|15x)|3\\.6\\.1\\.7x") ~ "From acetyl phosphate",
      str_detect(EC, "3\\.7\\.1\\.6x") ~ "From acetylpyruvate",
      str_detect(EC, "3\\.1\\.2\\.1x|6\\.2\\.1\\.(13x|1x)|2\\.8\\.3\\.(11x|12x|14x|-x|19x)") ~ "From acetyl-coa",
      str_detect(EC, "6\\.2\\.1\\.16x|4\\.1\\.1\\.4x") ~ "From acetoacetate",
      str_detect(EC, "3\\.1\\.1\\.55x") ~ "From asprin",
      str_detect(EC, "3\\.5\\.1\\.21x") ~ "Beta-alanine metabolism",
      str_detect(EC, "2\\.8\\.3\\.8x") ~ "From butanoic acid, acetyl-coa, or acetoacetate",
      str_detect(EC, "4\\.1\\.3\\.6x|2\\.8\\.3\\.10x") ~ "From citrate",
      str_detect(EC, "3\\.5\\.1\\.(41x|136x)") ~ "From chitin, chitobiose, chitosan",
      str_detect(EC, "3\\.1\\.1\\.45x") ~ "From chlorhexene, chlorobenzene",
      str_detect(EC, "3\\.1\\.1\\.(47x|53x|54x|58x|66x|71x|80x|106x)|3\\.1\\.2\\.(16x|12x)") ~ "Act on ester and carboxylic-ester bonds",
      str_detect(EC, "6\\.2\\.1\\.22x|6\\.2\\.1\\.35x|6\\.2\\.1\\.-x") ~ "Act on acid-thiol bonds",
      str_detect(EC, "3\\.7\\.1\\.24x") ~ "Act on ketones",
      str_detect(EC, "3\\.5\\.1\\.108x") ~ "Act on lipopolysaccharides",
      str_detect(EC, "1\\.13\\.12\\.4x") ~ "From lactate",
      str_detect(EC, "2\\.8\\.3\\.1x") ~ "From lactate, acetyl-coa",
      str_detect(EC, "4\\.1\\.3\\.22x") ~ "From methylmalate",
      str_detect(EC, "4\\.2\\.99\\.24x") ~ "From morphine",
      str_detect(EC, "2\\.3\\.1\\.187x|4\\.1\\.1\\.88x") ~ "From malonate",
      str_detect(EC, "3\\.5\\.1\\.(33x|25x)") ~ "From n-acetyl-d-glucosamine",
      str_detect(EC, "3\\.5\\.1\\.(76x|48x|85x|103x|104x|135x)") ~ "Act on linear amides",
      str_detect(EC, "3\\.5\\.1\\.89x") ~ "From oligosaccharidesl",
      str_detect(EC, "3\\.7\\.1\\.1x") ~ "From oxaloacetate",
      str_detect(EC, "1\\.13\\.11\\.50x") ~ "From pentane-2,4-dione",
      str_detect(EC, "3\\.11\\.1\\.2x") ~ "From phosphonoacetate",
      str_detect(EC, "1\\.2\\.5\\.1x") ~ "From pyruvate",
      str_detect(EC, "2\\.8\\.3\\.18x") ~ "From succinate",
      str_detect(EC, "3\\.7\\.1\\.25x") ~ "From xylans",
      str_detect(EC, "3\\.1\\.1\\.(80x|105x)|4\\.1\\.2\\.65x|1\\.14\\.11\\.11x") ~ "Alkaloid metabolism",
      str_detect(EC, "3\\.5\\.1\\.(112x|113x|-x)") ~ "Neomycin, kanamycin, gentamicin metabolism",
      str_detect(EC, "3\\.1\\.1\\.114x") ~ "Butanoate metabolism",
      str_detect(EC, "1\\.14\\.14\\.191x") ~ "Diterpenoid metabolism",
      str_detect(EC, "1\\.3\\.1\\.128x") ~ "Oxidoreductase reactions",
      str_detect(EC, "3\\.1\\.1\\.2x|4\\.1\\.2\\.(66x|65x)") ~ "Phenol metabolism",
      str_detect(EC, "1\\.1\\.1\\.(319x|318x)|2\\.3\\.1\\.-x") ~ "Phenylpropanoid metabolism",
      str_detect(EC, "1\\.14\\.13\\.54x|1\\.14\\.14\\.32x") ~ "Organic-aerobic degradation",
      str_detect(EC, "3\\.5\\.1\\.(29x|66x)") ~ "Vitamin B6 metabolism",
      str_detect(EC, "3\\.5\\.1\\.(47x|17x)") ~ "Lysine metabolism",
      str_detect(EC, "4\\.1\\.3\\.(4x|26x)") ~ "Pinene, camphor, geraniol metabolism",
      str_detect(EC, "3\\.1\\.1\\.94x") ~ "Aflatoxin metabolism",
      str_detect(EC, "2\\.1\\.1\\.152x|2\\.5\\.1\\.(47x|48x|49x)") ~ "Cystein, methionine & serine metabolism",
      str_detect(EC, "2\\.5\\.1\\.(50x|51x|52x|53x|125x|144x|119x|118x)|2\\.6\\.99\\.3x") ~ "Glycine, serine, threonine metabolism",
      str_detect(EC, "3\\.5\\.1\\.(17x|15x)") ~ "Alanine, aspartate, glutamate metabolism",
      str_detect(EC, "2\\.5\\.1\\.(65x|134x)") ~ "Cysteine & methoinine metabolism",
      str_detect(EC, "3\\.5\\.1\\.(14x|16x|62x|63x)") ~ "Arginine & proline metabolism"),
    Formate = case_when(
      str_detect(EC, "4\\.1\\.2\\.36x") ~ "From lactate",
      str_detect(EC, "2\\.7\\.2\\.6x") ~ "From formyl phosphate",
      str_detect(EC, "2\\.3\\.1\\.54x") ~ "From acetyl-coa",
      str_detect(EC, "4\\.1\\.99\\.(17x|23x)|1\\.13\\.11\\.(53x|54x)") ~ "Amino acid metabolism",
      str_detect(EC, "1\\.14\\.12\\.11x") ~ "Alkane metabolism",
      str_detect(EC, "1\\.13\\.11\\.89x") ~ "Phosphate metabolism",
      str_detect(EC, "4\\.1\\.2\\.44x") ~ "Benzoate metabolism",
      str_detect(EC, "4\\.1\\.99\\.12x") ~ "Riboflavin metabolism",
      str_detect(EC, "3\\.5\\.4\\.39x") ~ "Folate metabolism",
      str_detect(EC, "1\\.14\\.14\\.14x") ~ "Hormone biosynthesis",
      str_detect(EC, "6\\.3\\.4\\.(3x|17x|21x)|6\\.3\\.1\\.21x|1\\.3\\.5\\.1x|1\\.3\\.2\\.4x") ~ "From orthophosphate",
      str_detect(EC, "1\\.17\\.98\\.(5x|4x|3x)|1\\.17\\.1\\.(9x|10x|11x)|1\\.2\\.2\\.1x|1\\.17\\.2\\.3x|1\\.17\\.5\\.3x|1\\.1\\.98\\.6x|1\\.8\\.98\\.6x") ~ "From carbon dioxide"),
    Succinate = case_when(
      str_detect(EC, "1\\.3\\.1\\.6x|1\\.3\\.99\\.-x|1\\.3\\.98\\.1x|1\\.3\\.5\\.(1x|4x)|1\\.3\\.4\\.1x") ~ "From fumerate",
      str_detect(EC, "2\\.8\\.3\\.(13x|2x|-x|18x|20x)") ~ "TCA Cycle",
      str_detect(EC, "6\\.2\\.1\\.(5x|4x)|1\\.14\\.11\\.(1x|10x|18x|24x|55x|59x|36x|35x|37x|38x|40x|70x|71x|72x|75x|80x|81x|82x|-x|73x|58x|77x|78x|79x|64x|61x|51x|49x|48x|45x|44x|43x|42x)|1\\.14\\.20\\.(9x|7x|3x|10x|11x|12x|14x|16x)") ~ "Other oxogluterate-succinate conversions",
      str_detect(EC, "2\\.8\\.3\\.22x|3\\.1\\.2\\.3x") ~ "From succinyl-coa",
      str_detect(EC, "2\\.8\\.3\\.5x|3\\.5\\.1\\.(18x|96x)|1\\.14\\.11\\.(6x|2x|57x|7x|3x|8x|16x|41x|28x|39x|46x|56x|74x|76x)|1\\.2\\.1\\.(16x|24x|79x)|2\\.5\\.1\\.(48x|-x)|3\\.7\\.1\\.14x|1\\.14\\.20\\.7x") ~ "Amino Acid metabolism",
      str_detect(EC, "4\\.1\\.3\\.1x|1\\.1\\.1\\.93x|2\\.8\\.3\\.26x") ~ "Glyoxylate and dicarboxylate metabolism",
      str_detect(EC, "3\\.1\\.2\\.13x") ~ "Acting on ester bonds",
      str_detect(EC, "1\\.14\\.11\\.(27x|65x|66x|67x|68x|69x|63x)|1\\.14\\.20\\.-x") ~ "Protein metabolism",
      str_detect(EC, "1\\.14\\.11\\.17x") ~ "Taurine and hypotaurine metabolism",
      str_detect(EC, "1\\.14\\.11\\.(26x|52x)|1\\.14\\.20\\.1x") ~ "Antimicrobial metabolism",
      str_detect(EC, "1\\.2\\.1\\.76x|1\\.1\\.1\\.-x") ~ "Butanoate metabolism",
      str_detect(EC, "1\\.14\\.11\\.(20x|31x|32x)|1\\.14\\.20\\.13x") ~ "Alkaloid biosynthesis",
      str_detect(EC, "2\\.8\\.3\\.22x|1\\.14\\.11\\.21x") ~ "Organic acid metabolism",
      str_detect(EC, "2\\.8\\.3\\.6x|2\\.8\\.3\\.15x|2\\.3\\.1\\.315x|2\\.8\\.3\\.28x") ~ "Aromatic metabolism",
      str_detect(EC, "1\\.14\\.11\\.(13x|15x|12x)") ~ "Diterpenoid biosynthesis",
      str_detect(EC, "1\\.14\\.11\\.9x|1\\.14\\.20\\.(5x|6x|4x)") ~ "Flavenoid biosynthesis",
      str_detect(EC, "2\\.8\\.3\\.27x") ~ "Propanoate metabolism"),
    Propionate = case_when(
      str_detect(EC, "2\\.8\\.3\\.1x") ~ "From acetate, lactate",
      str_detect(EC, "4\\.1\\.3\\.32x") ~ "Nicotinate and nicotinamide metabolism",
      str_detect(EC, "6\\.2\\.1\\.(13x|1x|17x)|2\\.7\\.2\\.(1x|15x)|2\\.8\\.3\\.27x") ~ "Propanoate metabolism",
      str_detect(EC, "3.7.1.-x") ~ "Aromatic metabolism"),
    Butyrate = case_when(
      str_detect(EC, "2\\.3\\.1\\.16x|1\\.1\\.1\\.(35x|211x)|4\\.2\\.1\\.(17x|74x)|1\\.3\\.1\\.(38x|8x|125x)") ~ "core rBOX",
      str_detect(EC, "3\\.6\\.1\\.20x") ~ "Purine metabolism",
      str_detect(EC, "2\\.8\\.3\\.9x") ~ "Amino acid metabolism",
      str_detect(EC, "6\\.2\\.1\\.2x|2\\.7\\.2\\.7x|2\\.8\\.3\\.8x|1\\.3\\.1\\.31x|3\\.7\\.1\\.7x|3\\.1\\.1\\.51x") ~ "Butanoate metabolism"),
    Cellulose = case_when(
      str_detect(EC, "3\\.2\\.1\\.4x|2\\.4\\.1\\.49x|1\\.14\\.99\\.54x|1\\.14\\.99\\.56x") ~ "Degrades cellulose",
      str_detect(EC, "3\\.2\\.1\\.91x") ~ "Hydrolyse beta-linkages in cellulose & cellobiohydrolase",
      str_detect(EC, "3\\.2\\.1\\.(21x|23x)") ~ "Hydrolysis beta-D-glucosyl residues to beta-D-glucose",
      str_detect(EC, "3\\.2\\.1\\.(74x|31x)") ~ "Hydrolysis of O- and S-glycosyl compounds"),
    Hemicellulose = case_when(
      str_detect(EC, "3\\.2\\.1\\.(55x|22x|139x|78x)") ~ "Hydrolyse O- and S-glycosyl compounds",
      str_detect(EC, "3\\.2\\.1\\.(8x|7x|32x|37x|72x|131x|136x|156x)") ~ "Hydrolyze xylan",
      str_detect(EC, "3\\.1\\.1\\.(72x|73x)") ~ "Hydrolyze ester and carboxylic-ester bonds"),
    Galacturonan = case_when(
      str_detect(EC, "3\\.2\\.1\\.15x") ~ "Hydrolyse O- and S-glycosyl compounds",
      str_detect(EC, "4\\.2\\.2\\.(2x|9x)") ~ "Galacturonan into pectate)",
      str_detect(EC, "4\\.2\\.2\\.10x|2\\.4\\.1\\.43x") ~ "Make pectin",
      str_detect(EC, "3\\.1\\.1\\.11x") ~ "Make pectic acid"),
    Chito_oligosaccharides = case_when(
      str_detect(EC, "3\\.2\\.1\\.(14x|17x)|1\\.14\\.99\\.53x") ~ "Hydrolyze chitin",
      str_detect(EC, "3\\.2\\.1\\.(200x|201x)") ~ "Hydrolyse chitin and chitodextrins",
      str_detect(EC, "3\\.2\\.1\\.202x") ~ "Hydrolyse chitodextrins",
      str_detect(EC, "3\\.2\\.1\\.(132x|165x)") ~ "Hydrolyze chitosan",
      str_detect(EC, "3\\.5\\.1\\.41x") ~ "Hydrolyze chitin, chitobiose",
      str_detect(EC, "3\\.5\\.1\\.(105x|136x)") ~ "Hydrolyze chitin, produce acetate",
      str_detect(EC, "3\\.2\\.1\\.52x|2\\.7\\.1\\.196x") ~ "Hydrolyze chitobiose"),
    Peptidoglycan = case_when(
      str_detect(EC, "3\\.5\\.1\\.(28x|104x)|3\\.2\\.1\\.17x|4\\.2\\.2\\.29x") ~ "Hydrolyse peptidoglycan"),
    Fructans = case_when(
      str_detect(EC, "3\\.2\\.1\\.80x") ~ "Fructan to fructose",
      str_detect(EC, "3\\.2\\.1\\.65x") ~ "Levan hydrolysis",
      str_detect(EC, "3\\.2\\.1\\.26x") ~ "Hydrolyze beta-D-fructofuranosides"),
    Polysaccharides = case_when(
      str_detect(EC, "2\\.4\\.1\\.(1x|18x)|3\\.2\\.1\\.(10x|11x)") ~ "Hydrolyze starch and glucans",
      str_detect(EC, "3\\.2\\.1\\.6x") ~ "Glucans",
      str_detect(EC, "3\\.2\\.1\\.7x") ~ "Inulin",
      str_detect(EC, "3\\.2\\.1\\.24x") ~ "Glycan",
      str_detect(EC, "3\\.2\\.1\\.68x") ~ "Hydrolyze glycogen, amylopectin and their beta-limit dextrins",
      str_detect(EC, "3\\.2\\.1\\.(1x|2x|3x|20x)") ~ "Starch",
      str_detect(EC, "3\\.2\\.1\\.41x") ~ "Hydrolyze pullulan, amylopectin and glycogen")) %>%
  filter(!is.na(Hemicellulose)|!is.na(Galacturonan)|!is.na(Chito_oligosaccharides)|!is.na(Peptidoglycan)|!is.na(Fructans)|!is.na(Polysaccharides)|!is.na(Cellulose)|
           !is.na(Lactate)|!is.na(Ethanol)|!is.na(Acetate)|!is.na(Formate)|!is.na(Succinate)|!is.na(Propionate)|!is.na(Butyrate)) %>%
  pivot_longer(cols = c("Cellulose", "Hemicellulose", "Galacturonan", "Chito_oligosaccharides", "Peptidoglycan", "Acetate",
                        "Fructans", "Polysaccharides", "Lactate", "Ethanol", "Formate", "Succinate", "Propionate", "Butyrate"), 
               names_to = "Product", values_to = "Action", values_drop_na = TRUE) %>%
  group_by(Genome, Product, Enzyme_description,
           Quality, Completeness, Contamination, Alt_Tax_Label) %>%
  summarise(
    Type = ifelse(Product %in% Carbohydrates, "Carbohydrates", "Fermentation"),
    Frequency = n(),
    ECs = paste(unique(EC), collapse = ", ")) %>%
  distinct()
head(carbz)
#View(carbz)

relAb <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/abundance_mapping/coverm_relative_abundance_HQ_MQ.csv") %>%
  filter(Genome != "unmapped")
head(relAb)

Day_45 <- relAb %>% select(Genome, Y5_8_8) %>% mutate(Switch = "Pre-Switch", Day = "Day_45") %>% arrange(-Y5_8_8) %>% slice_head(n = 5)
Day_51 <- relAb %>% select(Genome, Y5_8_14) %>% mutate(Switch = "Pre-Switch", Day = "Day_51") %>% arrange(-Y5_8_14) %>% slice_head(n = 5)
Day_57 <- relAb %>% select(Genome, Y5_8_19E) %>% mutate(Switch = "Pre-Switch", Day = "Day_56") %>% arrange(-Y5_8_19E) %>% slice_head(n = 5)
Day_70 <- relAb %>% select(Genome, Y5_8_28) %>% mutate(Switch = "Post-Switch", Day = "Day_70") %>% arrange(-Y5_8_28) %>% slice_head(n = 5)
Day_77 <- relAb %>% select(Genome, Y5_9_2B) %>% mutate(Switch = "Post-Switch", Day = "Day_77") %>% arrange(-Y5_9_2B) %>% slice_head(n = 5)
Day_81 <- relAb %>% select(Genome, Y5_9_9B) %>% mutate(Switch = "Post-Switch", Day = "Day_81") %>% arrange(-Y5_9_9B) %>% slice_head(n = 5)
genomeListA <- unique(c(Day_45$Genome, Day_51$Genome, Day_57$Genome))
genomeListA
genomeListB <- unique(c(Day_70$Genome, Day_77$Genome, Day_81$Genome))
genomeListB

Day_45 <- relAb %>% select(Genome, Y5_8_8) %>% mutate(Switch = "B", Day = "Day_45") %>% arrange(-Y5_8_8)
Day_51 <- relAb %>% select(Genome, Y5_8_14) %>% mutate(Switch = "B", Day = "Day_51") %>% arrange(-Y5_8_14)
Day_57 <- relAb %>% select(Genome, Y5_8_19E) %>% mutate(Switch = "B", Day = "Day_57") %>% arrange(-Y5_8_19E)
Day_70 <- relAb %>% select(Genome, Y5_8_28) %>% mutate(Switch = "A", Day = "Day_70") %>% arrange(-Y5_8_28)
Day_77 <- relAb %>% select(Genome, Y5_9_2B) %>% mutate(Switch = "A", Day = "Day_77") %>% arrange(-Y5_9_2B)
Day_81 <- relAb %>% select(Genome, Y5_9_9B) %>% mutate(Switch = "A", Day = "Day_81") %>% arrange(-Y5_9_9B)


dayz <- Day_45 %>%
  full_join(Day_51) %>%
  full_join(Day_57) %>%
  full_join(Day_70) %>%
  full_join(Day_77) %>%
  full_join(Day_81) %>%
  pivot_longer(cols = c("Y5_8_8", "Y5_8_14", "Y5_8_19E", "Y5_8_28", "Y5_9_2B", "Y5_9_9B"), names_to = "Reactor_Name", values_to = "Relative_Abundance", values_drop_na = TRUE)
head(dayz)

carbz3 <- carbz %>%
  right_join(dayz, relationship = "many-to-many") %>%
  filter(Relative_Abundance > 0) %>%
  left_join(taxQC) %>%
  distinct() %>%
  mutate(
    Day = case_when(Day == "Day_45" ~ "45",
                    Day == "Day_51" ~ "51",
                    Day == "Day_57" ~ "57",
                    Day == "Day_70" ~ "70",
                    Day == "Day_77" ~ "77",
                    Day == "Day_81" ~ "81"),
    Switch = factor(Switch, levels = c("B", "A")),
    Product = factor(Product, levels = c("Cellulose", "Chito_oligosaccharides", "Fructans", "Galacturonan", "Hemicellulose", "Peptidoglycan", "Polysaccharides",       
                                         "Acetate", "Butyrate", "Ethanol", "Formate", "Lactate", "Propionate", "Succinate"))) %>%
  group_by(Type, Product, Switch, Day, Alt_Tax_Label, Relative_Abundance) %>%
  summarise(
    nMAGs = length(unique(Alt_Tax_Label)),
    nGenes = length(unique(Enzyme_description)),
    sumFreq = sum(Frequency)) 
head(carbz3)
# sum(carbz3$nMAGs)
# unique(carbz3$Product)
# range(carbz3$sumFreq)
#View(carbz3)

carbz3_A <- carbz3 %>% filter(Relative_Abundance > 0)
carbz3_B <- carbz3 %>% filter(Relative_Abundance > 3)
length(unique(filter(carbz3_A, Switch == "B")$Alt_Tax_Label))
length(unique(filter(carbz3_A, Switch == "A")$Alt_Tax_Label))
length(unique(filter(carbz3_B, Switch == "B")$Alt_Tax_Label))
length(unique(filter(carbz3_B, Switch == "A")$Alt_Tax_Label))

#carbz3_above5 <- carbz3 %>% filter(Relative_Abundance > 10)

colors <- c("tan", "goldenrod", "pink", "tomato", "tan4", "cyan4", "yellow2",       
            "purple4", "olivedrab", "grey35", "maroon4", "orange", "royalblue", "greenyellow")


p1 <- ggplot(carbz3_A, aes(x = Switch))+
  geom_boxplot(aes(y = sumFreq, fill = Product))+
  facet_nested(cols = vars(Product))+
  scale_fill_manual(values = colors) +
  ylim(0,16)+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.x = element_blank())+
  guides(fill = guide_legend(ncol = 7, byrow = TRUE))

p2 <- ggplot(carbz3_B, aes(x = Switch))+
  geom_boxplot(aes(y = sumFreq, fill = Product))+
  facet_nested(cols = vars(Product))+
  scale_fill_manual(values = colors) +
  ylim(0,16)+
  ylab(NULL)+
  xlab(NULL)+
  theme_minimal() +
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.5, hjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text.x = element_blank())+
  guides(fill = guide_legend(ncol = 7, byrow = TRUE))

# p3 <- ggplot(carbz3_above5, aes(x = Day))+
#   geom_boxplot(aes(y = sumFreq, fill = Product))+
#   facet_nested(cols = vars(Product), space = "free", shrink = TRUE, scales = "free")+
#   scale_fill_manual(values = colors) +
#   ylab(NULL)+
#   xlab(NULL)+
#   theme_minimal() +
#   theme(axis.text.x = element_blank(),#text(size = 16, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
#         axis.text.y = element_text(size = 16, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
#         legend.text = element_text(size = 16, color = "black"),
#         legend.position = "bottom",
#         legend.title = element_blank(),
#         strip.text.x = element_blank())+
#   guides(fill = guide_legend(ncol = 7, byrow = TRUE))

#http://127.0.0.1:26079/graphics/plot_zoom_png?width=1205&height=693
ggarrange(p1, p2, ncol = 1, common.legend = TRUE, legend = "bottom")

