# load libraries
library(csv)
library(tidyverse)
library(ggh4x)
library(forcats)
library(readr)
library(ggvegan2)
library(stringr)
#install.packages("remotes")
#remotes::install_github("jarioksa/ggvegan2")


# load data

taxQC <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/taxQC_info.csv") %>%
  mutate(user_genome = Genome, MAG = Genome)
head(taxQC)
#View(taxQC)

relAb <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/abundance_mapping/coverm_relative_abundance_HQ_MQ.csv")
colnames(relAb) <- c("MAG", "Y5_8_14", "Y5_8_19E", "Y5_8_28",  "Y5_8_8",   "Y5_9_2B",   "Y5_9_9B")
head(relAb)
#View(relAb)

clean_genes <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/metaErg_clean_genes_out_09192025.csv")
head(clean_genes)
dim(clean_genes)
#View(clean_genes)

clean_mag_ids <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/clean_mag_ids.csv")
head(clean_mag_ids)
dim(clean_mag_ids)

# genes_of_interest<- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/guilds-of-interest-08112025.csv") %>%
#   select(EC_number, KO_number, Gene_name, Metabolic_category, Warfare_compound, Function) %>%
#   mutate(Metabolic_category = ifelse(is.na(Metabolic_category), "", Metabolic_category),
#          Warfare_compound = ifelse(is.na(Warfare_compound), "", Warfare_compound)) %>%
#   filter(!is.na(EC_number)) %>%
#   group_by(EC_number) %>%
#   summarise(across(everything(), ~paste(unique(.), collapse = " ")), .groups = "drop")
# head(genes_of_interest)

##
#metacyc_pathways <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/metacyc_pathways.csv") %>%
# filter(!is.na(Metabolism))
#head(metacyc_pathways)

COG_test <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/merged_and_cleaned_eggnog_output_11192025.csv") %>%
  select(gene_id, Description, COG_category, GOs, EC) %>%
  filter(str_detect(gene_id, "##", negate = TRUE)) %>%
  separate(gene_id, into = c("Genome", "gene_identifier"), sep = ".faa_", remove = FALSE)
head(COG_test)

#head(bacteriocin_df$value)

bacteriocin_regex <- "bacteriocin|acidocin|cocin|durancin|enterocin|lactam|
lactococcin|lantibiotic|microcin|plantaricin|pediocin|garvicin|actifensin|actinomycesin|
amylocyclicin|aureocin|avidumicin|baceridin|bactofencin|brevicillin|buwchitin|cypemycin|
darobactin|doderlin|epidermicin|geobacillin|griselimycin|hominicin|huazacin|hyicin|ipomicin|
lacticin|lactocin|lactolisterin|laterosporulin|lichenin|medipeptin|mobilisin|
paenibacterin|pantocin|paralichenysindy|pentapeptide|pep27|peptide|perfrin|ruminococcin|sclerosin|
sonorensin|subtilosin|thurincin|tolworthcin|tryglysin|warnericin|delta.lysin|propionicin|
orniglycocin|sublancin|achromonodin-1|astexin-1|cloacaenodin|lariatin|propeptin|streptomonomicin|trilenodin|
ubonodin|cacaoidin|carnolysin|cerecidin|estercin|salivaricin|subtilomycin|thuricin|ticin|warnerin|plantazolicin|
glycocin|subtilin|siamycin|sviceucin|carnocin|duramycin|entianin|epicidin|epilancin|formicin|haloduracin|
lichenicidin|microbisporicin|mutacin|nisin|nukacin|paenibacillin|paenicidin|penisin|pseudocin|rodencin|
capistruin|chaxapeptin|citrocin|blauticin|cesin|epidermin|ericin|gallidermin|lan-ce02|lan-df|lan-p49.1|
pep5|planosporicin|rombocin|streptin|actagardine|bovicin|cinnamycin|cytolysin|mersacidin|michiganin|
roseocin|staphylococcin|streptococcin|suicin|thusin|variacin|carnobacteriocin|curvacin|curvaticin|
mesentericin|mundticin|penocin|piscicolin|leucocin|coagulin|divercin|divergicin|listeriocin|sakacin|abp-118|
cerecyclin|circularin|closticin|gassericin|lactocyclicin|leucocyclicin|plantaricyclin|pumilarin|uberolysin|
arnobacteriocin|sakacin|brochocin|gallocin|infantaricin|lactacin|thermophilin|
butyrivibriocin|carnocyclin|garvieacin"


head(clean_genes)

clean_filt <- clean_genes %>%
  mutate(Matches = str_extract_all(Enzyme_description, paste(bacteriocin_regex, collapse = ", ")),
         Matches = as.character(Matches))

pre_warfare_genes <- relAb %>%
  left_join(taxQC) %>%
  left_join(clean_genes) %>%
  left_join(clean_mag_ids) %>%
  left_join(COG_test) %>%
  filter(Y5_8_14>2|Y5_8_19E>2|Y5_8_28>2|Y5_8_8>2|Y5_9_2B>2|Y5_9_9B>2) %>%
  mutate(
    Enzyme_description = str_to_lower(Enzyme_description),
    Tax_Label = ifelse(Species == "Unclassified", Genus, Species)) %>%
  select(Genome, MAG, Completeness, Contamination, Quality, Family, Genus, Species, Alt_Tax_Label, Enzyme_description, KO, EC, Specific_pathway) %>%
  mutate(#Other_peptides = case_when(str_detect(Enzyme_description, bacteriocin_regex) ~ "Miscellaneous"),
    Pediocin = case_when(str_detect(Enzyme_description, "pediocin.*immunity") ~ "Resistance",
                         str_detect(Enzyme_description, "pediocin.*transport") ~ "Transport",
                         str_detect(Enzyme_description, "pediocin") ~ "Miscellaneous"),
    Helveticin = case_when(str_detect(Enzyme_description, "bacteriocin.*helveticin") ~ "Production",
                           str_detect(Enzyme_description, "helveticin") ~ "Miscellaneous"),
    Leucocin = case_when(str_detect(Enzyme_description, "leucocin.*immunity ") ~ "Resistance",
                         str_detect(Enzyme_description, "leucocin") ~ "Miscellaneous"),
    Lactococcin = case_when(str_detect(Enzyme_description, "lactococcin.*family") ~ "Production",
                            str_detect(Enzyme_description, "lactococcin.*trans") ~ "Transport",
                            str_detect(Enzyme_description, "lactococcin") ~ "Miscellaneous"),
    Bacteriocin = case_when(str_detect(KO, "K20344|K20345") ~ "Transport",
                            str_detect(Enzyme_description, "bacteriocin.*protection|peptidase.*c39|cysteine.*protease|peptidase.*M16|peptidase.*M20|peptide.*binding.*protein|AMP.*binding.*enzyme|multidrug.*resistance") ~ "Resistance",
                            str_detect(Enzyme_description, "bacteriocin.*biosynthesis|bacteriocin.*class|bacteriocin.*type|bacteriocin.*protein|bacteriocin.*sequence|peptide.*chain.*release")|str_detect(KO, "K20482") ~ "Production", 
                            str_detect(Enzyme_description, "peptide.*import|peptide.*transport|peptide.*transporter|multi.*ug.*export|multi.*ug.*transport|multi.*ug.*efflux") ~ "Transport",
                            str_detect(Enzyme_description, "bacteriocin") ~ "Miscellaneous"),
    Microcin = case_when(str_detect(KO, "K06159|K13893|K13894|K13895|K13896")|str_detect(Enzyme_description, "microcin.*secretion") ~ "Transport",
                         str_detect(Enzyme_description, "microcin.*immunity") ~ "Resistance",
                         str_detect(Enzyme_description, "microcin") ~ "Miscellaneous"),
    Enterocin = case_when(str_detect(Enzyme_description, "enterocin.*immunity") ~ "Resistance",
                          str_detect(Enzyme_description, "enterocin") ~ "Miscellaneous"),
    Subtilosin = case_when(str_detect(KO, "K27612|K27613|K27614|K27615|K27616|K27618") ~ "Production",
                           str_detect(Enzyme_description, "subtilosin") ~ "Miscellaneous"),
    Lantibiotic = case_when(str_detect(Enzyme_description, "lantibiotic transport")|
                              str_detect(KO, "K20459|K20460|K20461|K20490|K20491|K20492|K20494") ~ "Transport",
                            str_detect(Enzyme_description, "lantibiotic.*biosynthesis|lantibiotic.*dehydratase|lantibiotic.*dehydratase|lantibiotic.*alpha")|
                              str_detect(KO, "K20482|K20487|K20488")|
                              str_detect(EC, "2.7.13.3") ~ "Production",
                            str_detect(Enzyme_description, "lantibiotic.*protection")|
                              str_detect(KO, "K20486|K20489")|
                              str_detect(EC, "3.4.21.-") ~ "Resistance",
                            str_detect(Enzyme_description, "lantibiotic") ~ "Miscellaneous"),
    Subtilin = case_when(str_detect(Enzyme_description, "subtilase") ~ "Resistance",
                         str_detect(Enzyme_description, "subtilin.*transport") ~ "Transport",
                         str_detect(Enzyme_description, "subtilin") ~ "Miscellaneous"),
    Colicin = case_when(str_detect(KO, "K03558")|str_detect(Enzyme_description, "colicin.*production|colicin.*bacteriocin|phage.*colicin") ~ "Production",
                        str_detect(Enzyme_description, "colicin.*receptor|colicin.*uptake|colicin.*import|colicin.*secretion")|str_detect(KO, "K03646|K16089") ~ "Transport",
                        str_detect(Enzyme_description, "colicin.*immunity") ~ "Resistance",
                        str_detect(Enzyme_description, "colicin") ~ "Miscellaneous")) %>%
  filter(!is.na(Lactococcin)|!is.na(Leucocin)|!is.na(Enterocin)|!is.na(Subtilosin)|!is.na(Pediocin)|!is.na(Helveticin)|!is.na(Bacteriocin)|!is.na(Microcin)|!is.na(Subtilin)|!is.na(Colicin)|!is.na(Lantibiotic)) %>%
  mutate(Bacteriocin = ifelse(!is.na(Helveticin)|!is.na(Colicin)|!is.na(Lactococcin), NA, Bacteriocin)) 
head(pre_warfare_genes)
#View(pre_warfare_genes)
#write_csv(pre_warfare_genes, "/Users/lauraschaerer/Desktop/pre_warfare_genes_11262025.csv")

unique_proteins <- pre_warfare_genes %>%
  pivot_longer(cols = c("Enterocin", "Leucocin", "Subtilosin", "Pediocin", "Helveticin", "Bacteriocin", "Microcin", "Lantibiotic", "Lactococcin", "Subtilin", "Colicin"), names_to = "Peptide_name", values_to = "Mechanism_type", values_drop_na = TRUE) %>%
  filter(Mechanism_type != "Transport") %>%
  select(Enzyme_description, Mechanism_type, Peptide_name) %>%
  distinct() %>%
  left_join(warfare_genes3)
head(unique_proteins)
#write_csv(unique_proteins, "/Users/lauraschaerer/Desktop/feedstock_switch/uniQ.csv")

bacteriocin_info <- read_csv("/Users/lauraschaerer/Desktop/clean_bacteriocin_summary_12012025.csv")
head(bacteriocin_info)

uniQ <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/uniQ.csv") %>%
  select(Enzyme_description, Resistance_desc, Mechanism_type, Peptide_name) %>%
  distinct()
unique(uniQ$Resistance_desc)


warfare_genes2 <- pre_warfare_genes %>%
  pivot_longer(cols = c("Enterocin", "Leucocin", "Subtilosin", "Pediocin", "Helveticin", "Bacteriocin", "Microcin", "Lantibiotic", "Lactococcin", "Subtilin", "Colicin"), names_to = "Peptide_name", values_to = "Mechanism_type", values_drop_na = TRUE) %>%
  filter(Mechanism_type != "Transport") %>%
  filter(str_detect(Alt_Tax_Label, "JALCQ", negate = TRUE)) %>%
  mutate(Gram_stain = ifelse(str_detect(Alt_Tax_Label, "Xan|Prev|Dial", negate = FALSE), "G-", "G+")) %>%
  left_join(bacteriocin_info) %>%
  select(-Peptide_name) %>%
  left_join(uniQ) %>%
  mutate(Count = 1,
         Acts_against2 = ifelse(str_detect(Peptide_name, "ulti"), "M", Acts_against2),
         Mechanism_Type = factor(Mechanism_type, levels = c("Production", "Resistance")),
         Fill_Var1 = ifelse(Mechanism_type == "Production", paste0(Peptide_name, " (", Acts_against2, ")"), Acts_against2),
         Fill_Var1 = ifelse(Mechanism_type == "Resistance", paste0(Resistance_desc, " (", Acts_against2, ")"), Fill_Var1),
         Fill_Var = factor(Fill_Var1, levels = c("Lantibiotic (B)", "Colicin (G-)", "Helveticin (G+)", "Bacteriocin (NS)",
                                                "Degradation (M)", "Protection (M)",
                                                "Degradation (B)", 
                                                "Immunity (G+)",
                                                "Degradation (NS)", "Protection (NS)", "Binding (NS)", "Transport (NS)"))) %>%
  group_by(Alt_Tax_Label, Peptide_name, Mechanism_type, Acts_against2, Resistance_desc, Gram_stain, Fill_Var1, Fill_Var) %>%  
  summarise(Count = sum(Count))
head(warfare_genes2)
unique(warfare_genes2$Fill_Var1)
unique(warfare_genes2$Fill_Var)
#View(warfare_genes2)


colors <- c("darkgreen", "blue3", "gold", "grey44", 
            "maroon4", "pink3", 
            "olivedrab4", 
            "orange",
            "grey20", "grey30", "grey40", "grey50")


# colors <- c("grey50", "grey40", "grey30", "grey20",
#             "orange",
#             "olivedrab4",
#             "pink3", "maroon4",
#             "grey44", "gold", "blue3", "darkgreen")

#http://127.0.0.1:36716/graphics/plot_zoom_png?width=1425&height=740
ggplot(warfare_genes2, aes(x = Count, y = Alt_Tax_Label))+
  facet_nested(cols = vars(Mechanism_type), rows = vars(Gram_stain), space = "free", scales = "free")+
  geom_col(aes(fill = Fill_Var), color = "black", position = position_stack(reverse = TRUE))+
  theme_linedraw()+
  scale_fill_manual(values = colors)+
  theme(axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 1, color = "black"),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 18, angle = 0, hjust = 1, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 20, color = "black"),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        strip.text.x = element_text(size = 22, face = "bold", angle = 0, color = "black", hjust = 0),
        strip.text.y = element_text(size = 25, face = "bold", angle = 0, color = "black"))+
  guides(fill = guide_legend(ncol = 4, byrow = TRUE, reverse = FALSE))

