######### Classified species Barplots, GAMs and Timelines selcted species#############

### Libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(plotly)
library(ggplot2)
library(gratia)
library(mgcv)
library(stringr)
library(Biostrings)
library(plotly)

#### Load Data
all_data <- read.xlsx("Supplementary Data S3 Relative Abundance Both Genes.xlsx")
species_info <- read.xlsx("Supplementary Data S4_SpStatus.xlsx")
metadata <- read.xlsx("metadata.xlsx")

# Filter for species-level entries
all_data <- all_data %>% filter(!is.na(species) & species != "")

## Long format
rel_abun_long <- all_data %>%
  pivot_longer(cols = matches("^HB|^HM"), names_to = "SampleID", values_to = "Reads")

## Join taxonomy info
rel_abun_long <- left_join(rel_abun_long, species_info %>% dplyr::select(species, NewStatus), by = "species") %>%
  filter(!is.na(NewStatus))

## Add metadata
rel_abun_long <- left_join(rel_abun_long, metadata, by = "SampleID")

### Relative abundance NIS both genes
phylum_colors <- c(
  "Annelida"          = "#66C2A5",  
  "Discosea"          = "#BFD3E6",  
  "Evosea"            = "#F6CED8",  
  "Arthropoda"        = "#FFD92F",  
  "Chlorophyta"       = "#B3E2CD",  
  "Nematoda"          = "#A6D854",  
  "Ascomycota"        = "#FDB462",  
  "Cnidaria"          = "#A1D99B",  
  "Rhodophyta"        = "#CAB2D6",  
  "Oomycota"          = "#FFED6F",  
  "Bacillariophyta"   = "#BC80BD",  
  "Heterolobosea"     = "#CCEBC5",  
  "Mucoromycota"      = "#D9D9D9",  
  "Loukozoa"          = "#BCB0D9",  
  "Mollusca"          = "#E78AC3", 
  "Prasinodermophyta"= "#F0E442",  
  "Gyrista"           = "#80B1D3",  
  "Chordata"          = "#FC8D62",  
  "Tubulinea"         = "#B3B3B3",  
  "Gnathostomulida"   = "#B8860B",  
  "Myzozoa"           = "#CCEBC5",  
  "Cercozoa"          = "#D95F02", 
  "Haptophyta"        = "#E5C494" , 
  "Streptophyta"     = "#A6CEE3",
  "Endomyxa"         = "#B2DF8A",
  "Platyhelminthes"  = "#FDBF6F", 
  "Aphelidiomycota"  = "#CAB2D6",
  "Apicomplexa"      = "#FB9A99",
  "Euglenozoa"       = "#E31A1C",
  "Bryozoa"          = "#33A02C",
  "Choanozoa"        = "#1F78B4",
  "Nibbleridia"      = "#FF7F00",
  "Foraminifera"     = "#6A3D9A",
  "Perkinsozoa"      = "#B15928",
  "Sulcozoa"         = "#999999",
  "Ochrophyta" ="#9e0c29ff"
)

species_rel <- rel_abun_long %>%
  group_by(Age, phylum, Env) %>%
  summarise(TotalRelAbund = sum(Reads), .groups = "drop")

ggplot(species_rel, aes(x = TotalRelAbund, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = phylum_colors) +  
  labs(
    title = "Relative abundance by phylum NIS",
    x = "Reads",
    y = "Age",
    fill = "Kingdom"
  ) +
  theme_minimal()

### Richness NIS Both genes

status_palette <- c(
  "CRY"     = "#E69F00",  
  "NAT"     = "#009E73",  
  "NIS"     = "#56B4E9",  
  "Unknown" = "#7851A9"   
)

species_richness <- rel_abun_long %>%
  filter(!is.na(species), !is.na(NewStatus), Reads > 0) %>%
  distinct(SampleID, species, NewStatus, Age, Env) %>%
  group_by(Age, NewStatus, Env) %>%
  summarise(Richness = n_distinct(species), .groups = "drop")

ggplot(species_richness, aes(x = Age, y = Richness, color = NewStatus)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "gam", se = TRUE, formula = y ~ s(x, k = 10), linewidth = 0.5, alpha = 0.2) +
  scale_x_continuous(breaks = pretty(rel_abun_long$Age, n = 30)) +
  facet_wrap(~ Env, scales = "free_y", ncol=1) +
  labs(title = "Species Richness Metazoa and Plants Detected", x = "Age (cal BP)", y = "Richness", color = "Status") +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold")
  )

ggplot(species_richness, aes(x = Age, y = Richness, color=NewStatus, fill = NewStatus)) +
  geom_line(linewidth = 1) +
  geom_smooth(method = "gam", se = TRUE, formula = y ~ s(x, k = 4), linewidth = 0.5, alpha = 0.3) +
  scale_x_continuous(breaks = pretty(rel_abun_long$Age, n = 30)) +
  scale_y_continuous(breaks = pretty(species_richness$Richness, n = 10))+
  facet_wrap(~ Env, scales = "free_y", ncol=1) +
  labs(title = "Species Richness Metazoa and Plants Detected", x = "Age (cal BP)", y = "Richness", color = "Status") +
  theme_minimal()+
  scale_color_manual(values = status_palette) +
  scale_fill_manual(values = status_palette) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=12),
    axis.text.y = element_text(size=12),
    strip.text = element_text(face = "bold", size=18),
    axis.title.x = element_text(size=12),
    axis.title.y = element_text(size=12),
    title = element_text(size=20)
  )


### GAM stadistics

## Paper: https://www.frontiersin.org/journals/ecology-and-evolution/articles/10.3389/fevo.2018.00149/full

library(tidyverse)
library(mgcv)
library(gratia)

## ALL
species_richness_all <- species_richness

#GAM Model with smooth separated by NewStatus and Environment
species_rel_all <- as.factor(species_richness_all$NewStatus)
gam_model <- gam(Richness ~ NewStatus + Env +
                   s(Age, by = interaction(NewStatus, Env), k = 4),
                   data = species_richness_all, method = "REML")
summary(gam_model)
appraise(gam_model) .

plot(gam_model)

se_data <- smooth_estimates(gam_model)
glimpse(se_data)


##### Timeline selected NIS

all_data <- read.xlsx("RelAbund_AllRep_BOTH_Clean_Kingdom_nonsingl_asv_pooled.xlsx")
species_info <- read.xlsx("Supplementary Data S2.xlsx")
metadata <- read.xlsx("metadata.xlsx")



species_in <- species_info %>%
  filter(`Luke.100%.marine.NIS`=="x")

selected <- all_data %>%
  filter(species %in% species_in$species)

## Long format
selected_long <- selected %>%
  pivot_longer(cols = matches("^HB|^HM"), names_to = "SampleID", values_to = "Reads")

## Join taxonomy info
selected_long <- left_join(selected_long, species_info %>% select(species, NewStatus), by = "species") %>%
  filter(!is.na(NewStatus))

## Add metadata
selected_long <- left_join(selected_long, metadata, by = "SampleID")

selected_long <- selected_long %>%
  filter(NewStatus=="NIS")


### Relative abundance NIS both genes

species_rel <- selected_long %>%
  group_by(Age, Env, species) %>%
  summarise(TotalRelAbund = sum(Reads), .groups = "drop")

palette <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(unique(species_rel$species)))

ggplot(species_rel, aes(x = Age, y = log10(TotalRelAbund), color = species)) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = pretty(species_rel$Age, n = 30)) +
  scale_color_manual(values = palette) +
  facet_wrap(~ Env, ncol = 1, scales = "free_y") +
  labs(
    title = "Interesting NIS species",
    x = "Age (cal BP)",
    y = "Relative Abundance of Reads",
    color = "Species"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold")
  )



##### Figure 5 with detection points and replicate support #####

library(openxlsx)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

nis_summary <- read.xlsx("NIS_summary.xlsx")
metadata <- read.xlsx("metadata.xlsx")

focal_species <- c(
  "Neoporphyra haitanensis",
  "Nematostella vectensis",
  "Styela plicata",
  "Botrylloides violaceus"
)

abundance_cols <- names(nis_summary)[
  str_detect(names(nis_summary), "^(HB|HM)") &
    !str_detect(names(nis_summary), "_PCR_REPS$")
]

rep_cols <- names(nis_summary)[
  str_detect(names(nis_summary), "_PCR_REPS$")
]

abund_long <- nis_summary %>%
  filter(species %in% focal_species) %>%
  pivot_longer(
    cols = all_of(abundance_cols),
    names_to = "SampleID",
    values_to = "RelAbund"
  )

rep_long <- nis_summary %>%
  filter(species %in% focal_species) %>%
  pivot_longer(
    cols = all_of(rep_cols),
    names_to = "SampleID_rep",
    values_to = "PCR_REPS"
  ) %>%
  mutate(
    SampleID = str_remove(SampleID_rep, "_PCR_REPS$")
  ) %>%
  select(ASV, species, gene, SampleID, PCR_REPS)

# Join abundance + replicate support + metadata
nis_long <- abund_long %>%
  left_join(
    rep_long,
    by = c("ASV", "species", "gene", "SampleID")
  ) %>%
  mutate(
    RelAbund = replace_na(RelAbund, 0),
    PCR_REPS = replace_na(PCR_REPS, 0),
    Detection = ifelse(RelAbund > 0, "Present", "Absent")
  ) %>%
  left_join(metadata, by = "SampleID")

fig5_data <- nis_long %>%
  filter(RelAbund > 0)

pseudo <- min(fig5_data$RelAbund[fig5_data$RelAbund > 0], na.rm = TRUE) / 2

ggplot(fig5_data, aes(x = Age, y = log10(RelAbund + pseudo), colour = species)) +
  geom_line(aes(group = species), linewidth = 0.7, alpha = 0.6) +
  geom_point(aes(size = PCR_REPS), alpha = 0.9) +
  facet_grid(Env ~ ., scales = "free") +
  scale_size_continuous(
    range = c(2, 7),
    breaks = c(1, 2,3, 4, 5,6,7, 8),
    name = "PCR replicates"
  ) +
  labs(
    x = "Age",
    y = "log10(relative read abundance)",
    colour = "Species"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    strip.text = element_text(face = "bold")
  )

heatmap_data <- nis_long %>%
  mutate(
    PCR_REPS = replace_na(PCR_REPS, 0),
    Age = as.numeric(Age),
    Env = factor(Env, levels = c("Hythe", "Bursledon"))
  )
ggplot(heatmap_data, aes(x = factor(Age), y = species, fill = PCR_REPS)) +
  
  geom_tile(colour = "grey80", linewidth = 0.5) +
  
  geom_text(
    aes(label = ifelse(PCR_REPS > 0, PCR_REPS, "")),
    size = 6,
    colour = "black",
    fontface = "bold"
  ) +
  
  facet_grid(scales = "free_y", space = "free_x") +
  
  scale_fill_gradient(
    low = "white",
    high = "#0000ff99",
    name = "PCR replicates",
    breaks = 0:8,
    limits = c(0, 8)
  ) +
  
  labs(
    x = "Age",
    y = "NIS"
  ) +
  
  theme_minimal(base_size = 22) +
  guides(
    fill = guide_colorbar(
      barwidth = 2,
      barheight = 20,
      title.position = "top",
      title.hjust = 0.5
    )
  )+
  
  theme(
    
    strip.text = element_text(
      face = "bold",
      size = 24
    ),
    panel.grid.major = element_blank(),
panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(
      size = 24,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 24,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 18
    ),
    
    axis.text.y = element_text(
      size = 18,
      face = "italic"
    ),
    
    # LEGEND TITLE
    legend.title = element_text(
      size = 22,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 18
    ),
    
    panel.spacing = unit(1.2, "cm")
  )

ggsave(
  filename = "NIS_heatmap_large.svg",
  width = 18,
  height = 10,
  units = "in",
  dpi = 600,
  bg = "white"
)


heatmap_data2 <- heatmap_data %>%
  mutate(
    xmin = Age - 2,
    xmax = Age + 2,
    ymin = as.numeric(factor(species)) - 0.45,
    ymax = as.numeric(factor(species)) + 0.45
  )

ggplot(heatmap_data2) +
  
  geom_rect(
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = PCR_REPS
    ),
    colour = "grey80",
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      x = Age,
      y = as.numeric(factor(species)),
      label = ifelse(PCR_REPS > 0, PCR_REPS, "")
    ),
    size = 6,
    colour = "black",
    fontface = "bold"
  ) +
  
  scale_y_continuous(
    breaks = seq_along(unique(heatmap_data2$species)),
    labels = unique(heatmap_data2$species)
  ) +
  
  scale_x_continuous(
    breaks = seq(
      floor(min(heatmap_data2$Age, na.rm = TRUE) / 10) * 10,
      ceiling(max(heatmap_data2$Age, na.rm = TRUE) / 10) * 10,
      by = 10
    ),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  
  scale_fill_gradient(
    low = "white",
    high = "#0000ff99",
    name = "PCR replicates",
    breaks = 0:8,
    limits = c(0, 8)
  ) +
  
  labs(
    x = "Age",
    y = "NIS"
  ) +
  
  theme_minimal(base_size = 22) +
  
  guides(
    fill = guide_colorbar(
      barwidth = 2,
      barheight = 20,
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    
    axis.title.x = element_text(
      size = 24,
      face = "bold"
    ),
    
    axis.title.y = element_text(
      size = 24,
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 18
    ),
    
    axis.text.y = element_text(
      size = 18,
      face = "italic"
    ),
    
    legend.title = element_text(
      size = 22,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 18
    )
  )
