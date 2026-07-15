########## Bar Plots all relative abundance

all_data <- read.xlsx("Supplementary Data S3 Relative Abundance Both Genes.xlsx")
metadata <- read.xlsx("metadata.xlsx")

all_COI <- all_data %>%
  filter(gene=="COI")
all_18S<- all_data %>%
  filter(gene=="18S")

## Long format
all_COI_long <- all_COI %>%
  pivot_longer(cols = matches("^HB|^HM"), names_to = "SampleID", values_to = "Reads") 

all_18S_long <- all_18S %>%
  pivot_longer(cols = matches("^HB|^HM"), names_to = "SampleID", values_to = "Reads") 


all_COI_long <- left_join(all_COI_long, metadata, by = "SampleID")
all_18S_long <- left_join(all_18S_long, metadata, by = "SampleID")

## Kingdom palette
king_colors <- c(
  "Chromista"          = "#66c2a3ff",  
  "Fungi"          = "#bb12afff",  
  "Metazoa"            = "#8c9ecaff", 
  "Unassigned"        = "#f7ca00ff", 
  "Protist"       = "#98a286ff",  
  "Plantae"          = "#59bc2dff")  

ggplot(all_COI_long, aes(x = Reads, y = factor(Age), fill= kingdom)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = king_colors) +
  labs(
    title = "Relative abundance by kingdom - COI",
    x = "Reads",
    y = "Age",
    fill = "Kingdom"
  ) +
  theme_minimal()


ggplot(all_18S_long, aes(x = Reads, y = factor(Age), fill= kingdom)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = king_colors) +
  labs(
    title = "Relative abundance by kingdom - 18S",
    x = "Reads",
    y = "Age",
    fill = "Kingdom"
  ) +
  theme_minimal()


#### Plot Metazoa

COI_M <- all_COI_long %>% filter(kingdom=="Metazoa")
S18_M <- all_18S_long %>% filter(kingdom=="Metazoa")

## Phylum Palette
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
  "Sulcozoa"         = "#999999"
)

ggplot(COI_M, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = phylum_colors) +  
  labs(
    title = "Metazoa Phylums - COI",
    x = "Relative Abundance of Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


ggplot(S18_M, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = phylum_colors) +  labs(
    title = "Metazoa Phylums - 18S",
    x = "Relative Abundance of Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


#### Plot Metazoa and Plantae

COI_MP <- all_COI_long %>% filter(kingdom==c("Metazoa", "Plantae"))
S18_MP <- all_18S_long %>% filter(kingdom==c("Metazoa", "Plantae"))

ggplot(COI_MP, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = phylum_colors) +  
  labs(
    title = "Metazoa and Plantae Phylums - COI",
    x = "Relative Abundance of Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


ggplot(S18_MP, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  scale_fill_manual(values = phylum_colors) +  labs(
    title = "Metazoa and Plantae Phylums - 18S",
    x = "Relative Abundance of Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


### Plot Plantae

COI_P <- all_COI_long %>% filter(kingdom=="Plantae")
S18_P <- all_18S_long %>% filter(kingdom=="Plantae")

ggplot(COI_P, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  labs(
    title = "Plantae Relative abundance by kingdom - COI",
    x = "Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


ggplot(S18_P, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  labs(
    title = "Plantae Relative abundance by kingdom - 18S",
    x = "Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


#### Plot Protist

COI_Pr <- all_COI_long %>% filter(kingdom=="Protist")
S18_Pr <- all_18S_long %>% filter(kingdom=="Protist")

ggplot(COI_Pr, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  labs(
    title = "Protist Relative abundance by kingdom - COI",
    x = "Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


ggplot(S18_Pr, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  labs(
    title = "Protists Relative abundance by kingdom - 18S",
    x = "Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


#### Plot Fungi

COI_F <- all_COI_long %>% filter(kingdom=="Fungi")
S18_F <- all_18S_long %>% filter(kingdom=="Fungi")

ggplot(COI_F, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  labs(
    title = "Relative abundance by kingdom - COI",
    x = "Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()


ggplot(S18_F, aes(x = Reads, y = factor(Age), fill= phylum)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Env)+
  labs(
    title = "Relative abundance by kingdom - COI",
    x = "Reads",
    y = "Age",
    fill = "Phylum"
  ) +
  theme_minimal()

