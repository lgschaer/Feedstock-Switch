
# load libraries

library(csv)
library(tidyverse)
library(ggh4x)
library(stringr)
library(forcats)


#### Figure 2.A VFA Production


vfa <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/fatty acids/VFA_data_RYAN.csv")
head(vfa)

filt_vfa <- vfa %>%
  filter(Reactor == "Y5") %>%
  filter(ExpDay >= 43 & ExpDay <= 99) %>%
  mutate(ExpDay = as.character(ExpDay),
         ExpDay = factor(ExpDay, levels = c("43",  "45", "51",  "56",  "57",  "59",  "65",
                                            "70",  "71",  "73",  "79",  "84", "91",  "95",  "99"))) %>%
  pivot_longer(cols = c("Lactic", "Acetic","Butyric", "Valeric","Caproic", "Propionic", "Isovaleric", "Heptanoic", "Ethanol"), 
               names_to = "VFA", values_to = "g_L")
head(filt_vfa)
unique(filt_vfa$ExpDay)
#write_csv(filt_p3, "/Users/lauraschaerer/Desktop/nrel_project3_experiment2_sample_summary.csv")


filt_vfa2 <- vfa %>%
  filter(Reactor == "Y5") %>%
  filter(ExpDay >= 43 & ExpDay <= 99) %>%
  pivot_longer(cols = c("Lactic", "Acetic","Butyric", "Valeric","Caproic", "Propionic", "Isovaleric", "Heptanoic", "Ethanol"), 
               names_to = "VFA", values_to = "g_L") %>%
  mutate(ExpDay = as.integer(ExpDay))
head(filt_vfa2)
#View(filt_vfa2)

vlines <- dplyr::bind_rows(
  expand.grid(Reactor = c("Y5"),
              ExpDay  = 56, 
              color   = "M+PPB"),
  expand.grid(Reactor = c("Y5"), ExpDay  = c(70,99), color   = "LFW")) %>%
  mutate(Reactor_Desc = "Grass+Rumen->LFW")
head(vlines)
unique(vlines$Reactor)

old_colors <- c("blue", "maroon3", "purple", "green4", "cyan", "gold4", "royalblue", "firebrick", "gold", "grey20", "grey")
colors <- c("purple4", "olivedrab", "cyan3", "grey35", "lightblue", "firebrick", "orange", "royalblue", "pink3")

#http://127.0.0.1:26079/graphics/plot_zoom_png?width=1090&height=550
ggplot()+
  geom_vline(data = vlines, aes(xintercept = ExpDay, color = color), linetype = "dashed", size = 1, show.legend = TRUE)+
  geom_line(data = filt_vfa2, aes(y = g_L, x = ExpDay, color = VFA), size = 2, na.rm = TRUE, show.legend = FALSE) + 
  #ylab("VFA Concentration (g/L)") +
  #xlab("Time (days)")+
  xlab(NULL)+
  ylab(NULL)+
  #scale_color_manual(values = colors) +
  scale_x_continuous(breaks = seq(45, 95, by = 10)) +
  theme_minimal()+
  theme(axis.text.y = element_text(size = 20, color = "black"),
        axis.title.y = element_text(size = 23, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 1, color = "black"),
        axis.title.x = element_text(size = 23, color = "black"),
        legend.text = element_text(size = 18),
        legend.position = "bottom",
        legend.spacing.x = unit(0.1, 'mm'),
        legend.spacing.y = unit(0.05, 'mm'),
        plot.margin=grid::unit(c(0.1,0.1,0.1,0.1), "mm"),
        strip.text.y = element_text(size = 18, face = "bold", angle = 0, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.title = element_blank())#+
#guides(color = guide_legend(override.aes = list(linetype = "dashed")))

#### Figure 2.B Relative Abundance (MAGs)



# load data

taxQC <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/taxQC_info.csv") %>%
  mutate(user_genome = Genome)
head(taxQC)
#View(taxQC)

relAb <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/data/abundance_mapping/coverm_relative_abundance_HQ_MQ.csv")
head(relAb)
#View(relAb)

clean_mag_ids <- read_csv("/Users/lauraschaerer/Desktop/feedstock_switch/clean_mag_ids.csv")
head(clean_mag_ids)

relAb_tidy <- relAb %>%
  select(Genome, Y5_8_14, Y5_8_19E, Y5_8_28, Y5_8_8, Y5_9_2B, Y5_9_9B) %>%
  pivot_longer(cols = c("Y5_8_14", "Y5_8_19E", "Y5_8_28", "Y5_8_8", "Y5_9_2B", "Y5_9_9B"), names_to = "Reactor", values_to = "RelAb") %>%
  left_join(taxQC) %>%
  left_join(clean_mag_ids) %>%
  filter(RelAb != 0) %>%
  group_by(Reactor, Genome, Phylum, Family, Genus, Species, Alt_Tax_Label) %>%
  summarise(
    RelAb = sum(RelAb)) %>%
  mutate(Alt_Tax_Label = ifelse(is.na(Alt_Tax_Label), "Unmapped", Alt_Tax_Label))
head(relAb_tidy)
#View(relAb_tidy)

relAb_tmp <- relAb_tidy %>% arrange(-RelAb)
top19 <- head(unique(relAb_tmp$Alt_Tax_Label), 19)
top19

relAb_tidy2 <- relAb_tidy%>%
  mutate(Reactor = factor(Reactor, levels = c("Y5_8_8", "Y5_8_14", "Y5_8_19E", "Y5_8_28", "Y5_9_2B", "Y5_9_9B")), 
         Alt_Tax_Label = ifelse(Alt_Tax_Label %in% top19, Alt_Tax_Label, "Other Low Abundance"),
         Day = case_when(
           Reactor == "Y5_8_8" ~ "45",
           Reactor == "Y5_8_14" ~ "51",
           Reactor == "Y5_8_19E" ~ "57",
           Reactor == "Y5_8_28" ~ "70",
           Reactor == "Y5_9_2B" ~ "77",
           Reactor == "Y5_9_9B" ~ "81"),
         Switch = case_when(Day %in% c("45", "51") ~ "A",
                            Day == "57" ~ "B",
                            Day %in% c("70", "77", "81") ~ "C")) %>%         
  group_by(Reactor, Switch, Day, Alt_Tax_Label) %>%
  summarize(RelAb = sum(RelAb)) %>%
  mutate(Alt_Tax_Label = forcats::fct_other(Alt_Tax_Label, keep = c(top19, "Unmapped", "Other Low Abundance")),
         Alt_Tax_Label = forcats::fct_relevel(Alt_Tax_Label, "Unmapped", "Other Low Abundance")) #%>%
#  group_by(Alt_Tax_Label, Day, Reactor) %>%
#  summarize(
#    RelAb = sum(RelAb)
#  ) 
head(relAb_tidy2)
unique(relAb_tidy2$Alt_Tax_Label)
levels(relAb_tidy2$Alt_Tax_Label)



#write_csv(relAb_tidy, "/Users/lauraschaerer/Desktop/feedstock_switch/metagenome_taxonomy_abundance_08152025.csv")


#highly_abundant_taxa <- unique(relAb_tidy$Tax_Label)

colors <- c("white", "grey70", "yellow", "darkgreen", "green2",
            "orange", "red3", "pink", "lightskyblue", "dodgerblue1", "dodgerblue3",
            "dodgerblue4",  "purple", "navajowhite1", "navajowhite2",
            "navajowhite3", "navajowhite4", "tan4", "cyan", "olivedrab2")

#http://127.0.0.1:26079/graphics/plot_zoom_png?width=1750&height=550
ggplot(relAb_tidy2)+
  geom_col(aes(x = Day, y = RelAb, fill = Alt_Tax_Label), color = "black", position = "stack", na.rm = TRUE)+
  facet_grid(cols = vars(Switch), scales = "free", space = "free" )+
  scale_fill_manual(values = colors)+
  theme_linedraw()+
  xlab(NULL)+
  ylab(NULL)+
  scale_y_continuous(breaks = seq(0, 95, by = 10)) +
  theme_void()+
  theme(axis.title.y = element_text(size = 23, color = "black"),
        axis.title.x = element_text(size = 23, color = "black"),
        axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5, color = "black"),
        legend.text = element_text(size = 20, angle = 0, hjust = 0, vjust = 0.5, color = "black", margin = margin(l = 6, r = 6, t = 6, b = 6)),
        legend.position = "right",
        legend.justification = c(1, 0.1),
        #legend.key.size = unit(0.8, "cm"),
        legend.key = element_rect(fill = "white", color = "black", size = 0.7),
        legend.key.size = unit(1, "cm"),
        legend.title = element_blank()) +
  guides(fill = guide_legend( ncol = 2, byrow = FALSE, keywidth = unit(1, "cm"), keyheight = unit(1, "cm"), default.unit = "cm"))

